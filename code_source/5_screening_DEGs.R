#假如一个探针的ID对应多个基因,也没关系.因为我已加入if语句进行了条件判断
#此函数参数可接受data.frame和文件名称(character)
#多加了个参数outputworkdic,方便调用此函数,免得总是因为调取外面的变量而不能用此函数
screen_DEGs <- function(data 		=NULL,
						filename	='filtered_diff_result.csv',
						FC 			=NULL,
						p 			=NULL,
						outputworkdic=TRUE  ) {
#----------------------------------处理接口data来的数据---------------------------------------------
stopifnot(is.logical(outputworkdic))
if(class(data)=='character'){
 	if (file.exists('data')&file.info('data')['data','size']>1000){
	allData <- read.table(data, header = TRUE, sep = '\t', quote = '', dec = '.', na.strings = '---',row.names = NULL, skip = 0,
						   check.names = FALSE,blank.lines.skip = TRUE, strip.white = TRUE, comment.char = '#', fill = TRUE, as.is = TRUE)
	}
	else{
	stop('gene_data is not exist or the size of gene_data is wrong!')
			}
		}
else if(class(data)=='data.frame'){
	allData<-data
	}
#-----------------------------------完成了接口data数据处理---------------------------------------------------------
#-----------------------------------保留一个Gene Symbol---------------------------------------------
if(length(grep("; ", allData$`Gene Symbol`))!=0){#grep会返回个数字向量,所以自然它有长度
		#这两行很简单,把匹配多个基因的探针单独拿出来,
		 gsubed<-gsub(';.*','',allData$`Gene Symbol`)#此行有正则表达式,请注意理解
		 allData[,'Gene Symbol']<-gsubed
		unfilter_diff_result<-allData
												}
else if(length(grep("; ", allData$`Gene Symbol`))==0){
		unfilter_diff_result<-allData
		}
allData<<-allData#将allData设置成全局变量，目的是方便排查可能的错误
View(allData)
#-----------------------------------完成保留一个Gene Symbol---------------------------------------------
#-----------------------------------FC和P值筛选交互-----------------------------------------------------
# result<-0
# while(result!='y'){
			# cat("please set p value and FC value\nIf you don't want p or FC,please set 0\nIf the result of filtered is good,then please set y\n")
			# p<-as.numeric(readline(prompt='p value = '))
			# FC<-as.numeric(readline(prompt='FC value = '))
			# if(p==0&FC!=0){
				# filtered_diff_result<-subset(unfilter_diff_result,abs(unfilter_diff_result$`Fold Change`)>FC)
			# }
			# else if(p!=0&FC==0){
				# filtered_diff_result<-subset(unfilter_diff_result,unfilter_diff_result$`P-val`<p)
			# }
			# else if(FC!=0&p!=0){
				# filtered_diff_result<-subset(unfilter_diff_result,abs(unfilter_diff_result$`Fold Change`)>FC&unfilter_diff_result$`P-val`<p)#返回个data.frame
			# }
			# else if(FC==0&p==0){
				# filtered_diff_result<-unfilter_diff_result
			# }
			# View(filtered_diff_result)
			# result=readline(prompt='result of filtered is good? good set y bad set n\n\t')
		# }
filtered_diff_result<-subset(unfilter_diff_result,abs(unfilter_diff_result$`Fold Change`)>FC&unfilter_diff_result$`P-val`<p)#返回个data.frame

#-----------------------------------完成FC和P值筛选交互-----------------------------------------------------
#-----------------------------------加上UP&DOWN列,方便给人看------------------------------------------------
unfilter_diff_result[unfilter_diff_result$`Fold Change`>0,'UP&DOWN'] <- 'UP'
unfilter_diff_result[unfilter_diff_result$`Fold Change`<0,'UP&DOWN'] <- 'DOWN'
filtered_diff_result[filtered_diff_result$`Fold Change`>0, "UP&DOWN"] <- "UP"
filtered_diff_result[filtered_diff_result$`Fold Change`<0, "UP&DOWN"] <- "DOWN"
#-----------------------------------完成加上UP&DOWN列,方便给人看------------------------------------------------
#-----------------------------------给已经过滤的差异列表加log值------------------------------
data_FC1<-filtered_diff_result$"Fold Change"
FC_new1 <- cbind(data_FC1,rep(0,length(data_FC1)))
for(i in 1:nrow(FC_new1)){
	if (FC_new1[i,1]>0){
		FC_new1[i,2] <- FC_new1[i,1]
	}
	else if (FC_new1[i,1]<0){
		FC_new1[i,2] <-(-1/FC_new1[i,1])
		}
	}
FC_new1 <-FC_new1[,2]
log_FC_new1 <- log2(FC_new1)
#-----------------------------------完成给已经过滤的差异列表加log值------------------------------
#-----------------------------------给未过滤的差异列表加log值------------------------------------
data_FC2<-unfilter_diff_result$"Fold Change"
FC_new2 <- cbind(data_FC2,rep(0,length(data_FC2)))
for(i in 1:nrow(FC_new2)){
	if (FC_new2[i,1]>0){
		FC_new2[i,2] <- FC_new2[i,1]
	}
	else if (FC_new2[i,1]<0){
		FC_new2[i,2] <-(-1/FC_new2[i,1])
		}
	}
FC_new2 <-FC_new2[,2]
log_FC_new2 <- log2(FC_new2)
#-----------------------------------完成给未过滤的差异列表加log值------------------------------------
#注意:变量的名称虽是相同.但是内容已经不同了
filtered_diff_result<-cbind(filtered_diff_result,"log2(Fold Change)"=log_FC_new1)
unfilter_diff_result<-cbind(unfilter_diff_result,'log2(Fold Change)'=log_FC_new2)
# filtered_diff_result_demand_of_company
# unfilter_diff_result_demand_of_company
#-------------------------------------是否写入硬盘----------------------------------------------
if(outputworkdic==FALSE){
	# openxlsx::write.xlsx(unfilter_diff_result,file=paste(diff_gene_path,paste0(comp_name,'_diff_gene_unfilter.xlsx'),sep='\\'),asTable=T)
	# openxlsx::write.xlsx(filtered_diff_result,file=paste(diff_gene_path,paste0(comp_name,'_diff_gene_filtered.xlsx'),sep='\\'),asTable=T)
	write.csv(unfilter_diff_result,file=paste(diff_gene_path,paste0(comp_name,filename),sep='\\'),row.names=FALSE)
	write.csv(filtered_diff_result,file=paste(diff_gene_path,paste0(comp_name,filename),sep='\\'),row.names=FALSE)
	}
else if(outputworkdic==TRUE){
	# openxlsx::write.xlsx(filtered_diff_result,file=filename,asTable=TRUE)
	# openxlsx::write.xlsx(unfilter_diff_result,file=filename,asTable=TRUE)
	write.csv(filtered_diff_result,file=filename,row.names=FALSE)
	write.csv(unfilter_diff_result,file=filename,row.names=FALSE)
	}
else if(outputworkdic==NA){
	#如果outputworkdic==NA,则什么也不做
}
return(filtered_diff_result)
#-------------------------------------完成是否本地打印-------------------------------------------
}