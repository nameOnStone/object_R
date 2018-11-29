
#---------------------------------本脚本坚持几个原则-----------------------------------------
# 1.接口：坚持接口数据类型通用化，坚持接口接入后立即做接口检查，有可能的话应尽量窄化接口；
# 2.变量：坚持起名字时，增加数据类型，坚持词之间用“_”连接；
# 3.数据与画图之间的关系：坚持数据的处理与画图之间有明显的分割，避免数据或变量交叉以至于影响到测试；
# 4.包：坚持“不轻易升级包”的原则，以免脚本出现问题；
# 5.坚持可维护性、可扩展性、可读性、低耦合性；
#若未安装下列的包，则自动安装；
if(!requireNamespace('R6')){
  install.packages('R6')
}
if(!requireNamespace('tidyverse')){
  install.packages('tidyverse')
}
if(!requireNamespace('magrittr')){
  install.packages('magrittr')
}
if(!requireNamespace('svglite')){
  install.packages('svglite')
}
library('R6')
library('tidyverse')
#为增加代码简洁、易于理解，遂增加此包；
library('magrittr')
#为了产出svg格式图，增加此包；
library('svglite')
#I decide to build a class called God and put all 
#behavior that we needed into God's attribute;
God <- R6Class(classname = "God",
               public = list(
                 f_handled_data=NULL,
                 f_filtered_data=NULL,
                 c_species=NULL,
                 c_extractGene=NULL,
                 box_plot = NULL,
                 f_name_data = NULL,
                 f_handled_box_plot=NULL,
                 initialize = function(f_unhandled_data,f_name_data,c_species) {
                   #先将几个不合法的符号（如：‘ ’，‘-’，‘.’）改成‘_’，防止由于很多读进来的数据因为特殊符号问题而不能跑脚本；
                   colnames(f_unhandled_data) %<>% str_replace_all(string =,pattern = ' +',replacement = '_' )
                   colnames(f_unhandled_data) %<>% str_replace_all(string =,pattern = '\\.+',replacement = '_' )
                   colnames(f_unhandled_data) %<>% str_replace_all(string =,pattern = '-+',replacement = '_' )
                   colnames(f_name_data) %<>% str_replace_all(string =,pattern = ' +',replacement = '_' )
                   colnames(f_name_data) %<>% str_replace_all(string =,pattern = '\\.+',replacement = '_' )
                   colnames(f_name_data) %<>% str_replace_all(string =,pattern = '-+',replacement = '_' )
                   
                   #------------------------------------接口检查----------------------------------------
                   f_unhandled_data%>%{
                     stopifnot(all(is.data.frame(.),all(tibble::has_name(x =.,name = c('Fold_Change','P_val','FDR_P_val','Gene_Symbol')))))
                   }
                   
                   f_name_data%>%{
                     stopifnot(all(is.data.frame(.),
                                   all(tibble::has_name(x = .,name=c('Chip_ID','Sample_ID','Group_Name')))))
                   }
                   stopifnot(is.character(c_species))
                   #检查完成
                   
                   #接下来对数据源做必要的处理
                   #先改下名字
                   for(i in 1:nrow(f_name_data)){
                     colnames(f_unhandled_data) <- str_replace(string = colnames(f_unhandled_data),
                                                               pattern = paste0('^',f_name_data$Chip_ID[i],'.*ignal$'),
                                                               replacement = paste0(f_name_data$Sample_ID[i],'_signal'))
                   }
                   
                   #增加UP_DOWN列
                   f_unhandled_data[f_unhandled_data$`Fold_Change`>=0,'UP_DOWN']<-'UP'
                   f_unhandled_data[f_unhandled_data$`Fold_Change`<0,'UP_DOWN']<-'DOWN'
                   #通过system.time函数测试，使用sapply函数的CPU时间是for循环的250倍以上，所以改用sapply
                   #增加log值列
                   #message('The Handing need some times,please wait patiently')
                   #for(j in 1:length(f_unhandled_data$`Fold_Change`)){
                     #若 Fold_Change为正值，则直接取log2值
                     #if(f_unhandled_data$`Fold_Change`[j]>=0){
                       #f_unhandled_data[f_unhandled_data$`Fold_Change`[j],'log2(Fold_Change)'] <- log2(f_unhandled_data$`Fold_Change`[j])
                     #}
                     #若 Fold_Change为负值，则需取负倒数后才可取log2值
                    # else if(f_unhandled_data$`Fold_Change`[j]<0){
                       #f_unhandled_data[f_unhandled_data$`Fold_Change`[j],'log2(Fold_Change)'] <- log2(-1/f_unhandled_data$`Fold_Change`[j])
                     #}
                   #}
                   #通过system.time函数测试，使用sapply函数的CPU时间是for循环的250倍以上，所以改用sapply
                   f_unhandled_data$Fold_Change %>% sapply(X = ,FUN = function(x){
                     if(x>=0){log2(x)}
                     else if(x<0){log2(-1/x)}
                   }) -> f_unhandled_data[,'log2(Fold_Change)']
                   
                   #以上完成了必要的处理；
                   #把处理完成的数据源赋值给实例属性，方便引用；
                   self$f_handled_data <- f_unhandled_data
                   self$f_name_data <- f_name_data
                   self$c_species <- c_species
                 },#初始化函数结束；
                 #我们得先创建文件夹，这样比较合乎逻辑
                 create_report_dir=function(){
                   #暂时留空，方便后续添加代码
                   
                 },#创造文件夹函数结束；
                 create_zhenzhong_report_dir=function(nameofreport='TEST',pathofmkdir='.',comp_name='test_compare_name'){
                   current_path <- getwd()#不应变换工作目录
                   setwd(pathofmkdir)#两个作用：1，检查路径是否存在（如若不存在，则会直接报错） 2，setwd本身的作用
                   dir.create(paste(pathofmkdir,nameofreport,sep='/'), showWarnings = FALSE)
                   main_path <- paste(pathofmkdir, nameofreport, sep="/")
                   dir.create(paste(main_path,"1_原始数据分析实验相关内容部分","1_样本信息",sep='/'), recursive = T,showWarnings=FALSE)
                   dir.create(paste(main_path,"1_原始数据分析实验相关内容部分","2_芯片质控",sep='/'), recursive = T,showWarnings=FALSE)
                   dir.create(paste(main_path,"1_原始数据分析实验相关内容部分","3_实验数据",sep='/'), recursive = T,showWarnings=FALSE)
                   dir.create(paste(main_path,"2_探索性分析","1_箱线图",sep='/'), recursive = T,showWarnings=FALSE)
                   dir.create(paste(main_path,"2_探索性分析","2_Peason相关系数图",sep='/'), recursive = T,showWarnings=FALSE)
                   dir.create(paste(main_path,"2_探索性分析","3_PCA图",sep='/'), recursive = T,showWarnings=FALSE)
                   dir.create(paste(main_path,"3_mRNA分析",sep='/'), showWarnings = FALSE)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,sep='/'),showWarnings=FALSE)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,"1_差异基因",sep='/'),showWarnings=FALSE)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Barplot',sep='/'),showWarnings=FALSE,recursive=T)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Enrishment',sep='/'),showWarnings=FALSE,recursive = T)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Result',sep='/'),showWarnings=FALSE,recursive = T)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG_Result",sep='/'),showWarnings=FALSE,recursive = T)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG_Barplot & Dotplot",sep='/'),showWarnings=FALSE,recursive = T)
                   dir.create(paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","Pathway",sep='/'),showWarnings=FALSE,recursive = T)
                   dir.create(paste(main_path,"4_注释信息及说明文档",sep='/'), showWarnings = FALSE)
                   setwd(current_path)#转回原来的工作目录
                   li_mkdir<<-list(#将li_mkdir声明成为全局变量，目的是方便后面的方法可以调用此路径
                     path_exploratory_BoxWhisker=paste(main_path,"2_探索性分析","1_箱线图",sep='/'),
                     path_exploratory_PCA=   paste(main_path,"2_探索性分析","3_PCA图",sep='/'),
                     path_GO_Barplot=      paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Barplot',sep='/'),
                     path_GO_Enrishment=   paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Enrishment',sep='/'),
                     path_GO_result=     paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Result',sep='/'),
                     path_KEGG_Result=     paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG Result",sep='/'),
                     path_KEGG_Barplot=    paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG Barplot & Dotplot",sep='/'),
                     path_KEGG_pathway=    paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","Pathway",sep='/'),
                     path_diff_gene=     paste(main_path,'3_mRNA分析',comp_name,"1_差异基因",sep='/')
                   )
                 },
                 #箱线图画图方法
                 draw_box_plot=function(data=self$f_handled_data,pathofdownload=li_mkdir$path_exploratory_BoxWhisker){
                   #画box_plot图，需要二组数据源，1：f_handled_data(用来获取每组样品的数值),2:f_name_data(用来获取分组信息)
                   #----------------------------------------接口检查-----------------------------------------
                   stopifnot(all(is.data.frame(data)))
                   if(!any(is.na(pathofdownload),is.null(pathofdownload))){port <- T}
                   else(port <- 'empty')
                   #magrittr包的代码规则，利于节省中间变量，类似于linux中的管道符，不过貌似不熟悉的人感觉上挺难看懂。鱼和熊掌难兼得
                   #-------------------------------------数据源开始处理--------------------------------------------------
                   c_extract_samplename <- self$f_handled_data %>% colnames() %>% grep(x=,pattern = '.*signal$',value = T) 
                   
                   f_reconstruct <- NULL
                   for (i in c_extract_samplename){
                     bb <- data_frame(Name=i,Value=(self$f_handled_data[[i]]))
                     f_reconstruct <- rbind(f_reconstruct,bb)
                     }#循环结束
                   #我需要根据f_reconstruct中的Name列再次增加一列Group，方便在画bar_plot时，可以区分把分组信息用颜色标记出来
                   #我是根据f_name_data和f_reconstruct此两组数据框中的信息作出增加f_reconstruct的第三列的Group组；
                   message('The handling data need some times,Please wait!')
                   str_remove(string = f_reconstruct$Name,pattern = '_signal')%>%
                   sapply(X = .,FUN = function(i){
                        self$f_name_data[[which(self$f_name_data$Sample_ID==i),'Group_Name']]
                           })->f_reconstruct[,'Group']
                   #保留小数点后两位
                   f_reconstruct$Value <- round(f_reconstruct$Value,2)
                   #要画图了，所以还是要把signal这个字符串尾巴去掉；
                   f_reconstruct$Name <- str_remove(f_reconstruct$Name,'_signal')
                   #-------------------------------------数据源处理完成----------------------------------------
                   
                   #以上代码已经完成了数据处理，应当将处理好的数据源赋值给实例属性，以便后续有人需要更改图片时可调用画图的数据源
                   #不过，这里有一点需要注意，画box_plot不应当是成对画的，而应该是整个项目的数据一起画才对，所以此处代码从逻辑上是不是有问题的呢？
                   #仍需要厘清。暂时先这样，留下Mark，方便提醒自己。
                   self$f_handled_box_plot <- f_reconstruct
                   #ggplot的接口只接受data.frame,所以,输入前不妨检查下接口的数据类型
                   #----------------------------------------接口检查----------------------------------------------
                   stopifnot(all(is.data.frame(f_reconstruct),
                                 is.numeric(f_reconstruct$Value),
                                 is.character(f_reconstruct$Name),
                                 is.character(f_reconstruct$Group)))
                   #--------------------------------------开始画图---------------------------------------------
                   box_plot <- ggplot(data = f_reconstruct,aes(x = Name,y = Value,color=Group))+geom_boxplot()+theme_light()
                   if(port==T){
                    ggsave(nameofpic,path = pathofdownload)
                     message(' has been completed!')
                   }
                   #--------------------------------------结束画图---------------------------------------------
                   #将画图的结果赋值给实例属性，方便查看调用
                   self$box_plot <- box_plot
                   return(box_plot)
                 },#draw_box_plot方法结束；
                 #主成分分析画图方法
                 draw_prcomp_plot=function(data=self$f_handled_data,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 #皮尔森图画图方法
                 draw_pearson_plot=function(data=self$f_handled_data,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 
                 #以上是步骤过滤之前的方法
                 
                 #过滤方法
                 filter_data=function(data=self$f_handled_data,FC=1,p_value=0.05){# mo ren zhi
                   #暂时留空，方便后续添加代码
                 },
                 
                 #以下是步骤过滤之后的方法
                 draw_cluter_plot=function(data=self$f_filtered_data,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                   unigene <- unique(data[which(data[,"Gene Symbol"]!= "NA" & data[,"Gene Symbol"]!= ""),"Gene Symbol"])
                 },
                 draw_volcano_plot=function(data=self$f_filtered_data,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 
                 
                 #因为本质上做GO和KEGG的数据源都是基因或探针ID，因此，以下方法的输入接口宜做 
                 #改变成字符串向量。
                 
                 #从f_filtered_data中提取基因或探针名称的字符串向量
                 fun_extract_Gene=function(data=f_filtered_data,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 fun_enrich_GO=function(data=c_extractGene,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 
                 #以下方法都是做完GO分析之后的方法，所以，一般情况在做GO之前此方法不宜调
                 #除非只是做测试。
                 draw_bar_plot_GO=function(data=GOresult,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 fun_plot_GO=function(data=GOresult,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 
                 
                 #富集KEGG
                 #一般富集试验之间宜做geneID和Gene Symbol之间的转换，本人特意留住此方法。
                 fun_biter=function(data=c_extractGene,pathofdownload='.'){
                   library('pathview')
                   library('DOSE')
                   library('clusterProfiler')
                   
                   
                 },
                 #富集
                 fun_enrich_KEGG=function(data=c_extractGene,pathofdownload='.'){
                   #暂时留空，方便后续添加代码
                 },
                 
                 #以下方法都是在做完KEGG之后宜调用，因此不建议在enrichKEGGfun调用之前做。除非是做测试。
                 draw_bar_plot_KEGG=function(data=NULL){
                   
                 },
                 draw_dot_plot_KEGG=function(data=NULL){
                   #暂时留空，方便后续添加代码
                 },
                 view_KEGG1_fun=function(data=NULL){
                   #暂时留空，方便后续添加代码
                 },
                 get_html_report=function(data=NULL){
                   #暂时留空，方便后续添加代码
                 }
               )
)


#以下是继承God类的类
RNA <- R6Class(classname = 'RNA',inherit = God)