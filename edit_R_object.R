
#---------------------------本脚本坚持几个原则—----------------------------------------
# 1.接口：坚持接口数据类型通用化，坚持接口接入后立即做接口检查，有可能的话应尽量窄化接口；
# 2.变量：坚持起名字时，增加数据类型，坚持词之间用“_”连接；
# 3.数据与画图之间的关系：坚持数据的处理与画图之间有明显的分割，避免数据或变量交叉以至于影响到测试；
# 4.包：应坚持“不轻易升级包”的原则，以免脚本出现问题；
#---------------------------本脚本坚持的几个原则---------------------------------------
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
library('R6')
library('tidyverse')
#为增加代码简洁、易于理解，遂增加此包
library('magrittr')
#I decide to build a class called God and put all 
#behavior that we needed into God's attribute
God <- R6Class(classname = "God",
               public = list(
                 f_handled_data=NULL,
                 f_filtered_data=NULL,
                 c_species=NULL,
                 c_extractGene=NULL,
                 box_plot = NULL,
                 initialize = function(f_unhandled_data,f_name_data,c_species) {
                   #先将几个不合法的符号（如：‘ ’，‘-’，‘.’）改成‘_’，防止由于很多读进来的数据因为特殊符号问题而不能跑脚本；
                   colnames(f_unhandled_data) %<>% str_replace_all(string =,pattern = ' +',replacement = '_' )
                   colnames(f_unhandled_data) %<>% str_replace_all(string =,pattern = '\\.+',replacement = '_' )
                   colnames(f_unhandled_data) %<>% str_replace_all(string =,pattern = '-+',replacement = '_' )
                   #接口检查
                   stopifnot(is.data.frame(f_unhandled_data)&is.data.frame(f_name_data)&has_name(f_unhandled_data,'Fold_Change')&all(has_name())&is.character(c_species))
                   f_unhandled_data%>%{
                     stopifnot(all(is.data.frame(.),
                                   all(has_name(x = .,which = c('Fold_Change','P_val','FDR_P_val','Gene_Symbol')))))
                     
                   }
                   f_name_data%>%{
                     stopifnot(all(is.data.frame(.),
                                   all(has_name(x=.,which=c('Chip_ID','Sample_ID','Group_Name')))))
                   }
                   #检查完成
                   
                   #接下来对数据源做必要的处理
                   #先改下名称
                   for(i in 1:nrow(f_name_data)){
                     colnames(f_unhandled_data) <- str_replace(string = colnames(f_unhandled_data),pattern = paste0('^',f_name_data$Chip_ID[i],'.*ignal$'),replacement = paste0(f_name_data$Sample_ID,'_signal'))
                   }
                   #增加UP&DOWN列
                   f_unhandled_data[f_unhandled_data$`Fold_Change`>=0,'UP_DOWN']<-'UP'
                   f_unhandled_data[f_unhandled_data$`Fold_Change`<0,'UP_DOWN']<-'DOWN'
                   #增加log值列
                   for(j in 1:length(f_unhandled_data$`Fold_Change`)){
                     #若 Fold_Change为正值，则直接取log2值
                     if(f_unhandled_data$`Fold_Change`[j]>=0){
                       f_unhandled_data[f_unhandled_data$`Fold_Change`[j],'log2(Fold_Change)'] <- log2(f_unhandled_data$`Fold_Change`[j])
                     }
                     #若 Fold_Change为负值，则需取负倒数后才可取log2值
                     else if(f_unhandled_data$`Fold_Change`[j]<0){
                       f_unhandled_data[f_unhandled_data$`Fold_Change`[j],'log2(Fold_Change)'] <- log2(-1/f_unhandled_data$`Fold_Change`[j])
                     }
                   }
                   #以上完成了必要的处理；
                   #把处理完成的数据源赋值给实例属性，方便引用；
                   self$f_handled_data <- f_unhandled_data
                   self$c_species <- c_species
                 },
                 #我们得先创建文件夹，这样比较合乎逻辑
                 create_report_dir=function(){
                   #暂时留空，方便后续添加代码
                   
                 },
                 create_zhenzhongtuijian_report_dir=function(nameofreport='TEST',pathofmkdir='.',comp_name='test_compare_name'){
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
                   #接口检查
                   stopifnot(all(is.data.frame(data)))
                   if(!any(is.na(pathofdownload),is.null(pathofdownload))){port <- T}
                   #magrittr包的代码规则，利于节省中间变量,类似于linux中的管道符，不过貌似不熟悉的人感觉上挺难看懂。鱼和熊掌难兼得
                   c_extract_samplename <- test$f_handled_data %>% colnames() %>% grep(x=,pattern = '.*signal$',value = T) 
                   f_reconstruct <- NULL
                   for (i in c_extract_samplename){
                     print(i)
                     bb <- data_frame(value=(test$f_handled_data[[i]]),name=i,group='a')
                     f_reconstruct <- rbind(bb,f_reconstruct)
                   }
                   #保留小数点后两位
                   f_reconstruct$value <- round(f_reconstruct$value,2)
                   f_reconstruct$name <- str_replace(f_reconstruct$name,'_signal','')
                   aa <- str_replace(string = aa,'_signal','')
                   f_reconstruct[name==a,'groupname']
                   #ggplot的接口只接受data.frame,所以,输入前不妨检查下接口的数据类型
                   stopifnot(all(is.data.frame(f_reconstruct),is.numeric(f_reconstruct$value),is.character(f_reconstruct$name)))
                   #画图
                   box_plot <- ggplot(data = f_reconstruct,aes(x = name,y = value,color=))+geom_boxplot()+theme_light()
                   if(port==T){ggsave('box_plot.svg',path = pathofdownload)}
                   #将画图的结果赋值给实例属性，方便查看调用
                   self$box_plot <- box_plot
                   
                 },
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
                 tidy_data=function(data=self$f_handled_data,FC=1,p_value=0.05){# mo ren zhi
                   
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