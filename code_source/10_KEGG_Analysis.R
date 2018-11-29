#如果outputworkdic=FAlSE的话,则不必给filename这个参数
KEGG_Analysis <- function(	data     		= diff_result,
							filename 		= 'KEGGresult.csv',
							orgdb    		= "org.Hs.eg.db",
							species  		= "hsa",
							minGsize 		= 1,
							p        		= 1,
							q        		= 1,
							outputworkdic	=TRUE) {
	# library("RDAVIDWebService")
	library("DOSE")
	library("clusterProfiler")
	library("pathview")
	library('openxlsx')

	stopifnot(is.logical(outputworkdic))
	
	#convert the format of gene ID
	unigene <- unique(data[which(data[,"Gene Symbol"]!= "NA" & data[,"Gene Symbol"]!= ""),"Gene Symbol"])
	unigene<-as.character(unigene)
	loc1 <- match(unigene, data$`Gene Symbol`)
	diff_gene <- data[loc1, c("Gene Symbol", "log2(Fold Change)")]
	
	geneID <- bitr(			gene = diff_gene$`Gene Symbol`,
				 		fromType = "SYMBOL",
						  toType = "ENTREZID",
						   OrgDb = orgdb)
	
	tryresultKEGG<-tryCatch({ enrichKEGG(gene              = geneID$`ENTREZID`,
										keyType            = "kegg",
										organism           = species,
										pAdjustMethod      = "fdr",#KEGG若没有结果,可改此处为none ,annotation of wangduolin 
										minGSSize          = minGsize,
										pvalueCutoff       = p,
										qvalueCutoff       = q,
										use_internal_data  = FALSE)},error=function(e) {
		stop('KEGG analysis 模块中enrichKEGG函数出问题了~ ')
		})
	KEGG_result <- as.data.frame(tryresultKEGG)
	if (dim(KEGG_result)[1]==0){
		stop('please attention : The KEGG_result is empty,and you need to handle this problem~~~')#annotation of wangduolin 
	}
	diff_gene <- as.matrix(KEGG_result[,"geneID"])
	
	gene_symbol<-NULL
		for(i in 1:length(diff_gene)){
		first<-as.numeric(unlist(strsplit(diff_gene[i,],"/")))
		 geneID <- bitr(	gene = first,
						fromType = "ENTREZID",
						  toType = "SYMBOL",
						   OrgDb = orgdb)
		second<-geneID[,"SYMBOL"]
		gene_symbol<-rbind(gene_symbol,as.matrix(paste(second,collapse="/")))
		 }
		KEGG_result<-cbind(KEGG_result[,1:7],"Gene Symbol"=gene_symbol,"Count"=KEGG_result[,9])
	#write the output file
	a<-list("KEGG"=tryresultKEGG, "KEGG_result"=KEGG_result)
	if(outputworkdic==TRUE){
		#openxlsx::write.xlsx(x=a['KEGG_result'],file=filename,asTable=TRUE)
		write.csv(KEGG_result,file='KEGG_reuslt.csv',row.names=FALSE)
	}
	else if(outputworkdic==FALSE){
		write.csv(KEGG_result,file=paste(KEGG_Result_path, paste0(comp_name, '_',filename),sep  = "\\"),row.names=FALSE)
		# openxlsx::write.xlsx(x      = KEGGresult['KEGG_result'], 
							# file    = paste(KEGG_Result_path, paste0(comp_name, '_',filename),sep  = "\\"), 
					 		# asTable = TRUE)
	}
	else if(outputworkdic==NA){
		#如果outputworkdic==NA,则什么也不做
	}
	return(a)
}