GO_Analysis <- function(data     		= diff_result,
                        filename 		= "GOresult.xlsx",
                        orgdb    		= "org.Hs.eg.db",
                        minGsize 		= 1,
                        p   			= 1,
                        q      			= 1,
                        outputworkdic	= TRUE) {

	# library("RDAVIDWebService")
	library("DOSE")
	library("clusterProfiler")

	stopifnot(is.logical(outputworkdic))
	unigene <- unique(data[which(data[,"Gene Symbol"]!= "Ndata" & data[,"Gene Symbol"]!= ""),"Gene Symbol"])
	BP <- enrichGO(gene          = unigene,
		           OrgDb         = orgdb,
	               keyType       = "SYMBOL",
	               ont           = "BP",
	               minGSSize     = minGsize,
	               pAdjustMethod = "BH",
	               pvalueCutoff  = p,
	               qvalueCutoff  = q)
	CC <- enrichGO(gene          = unigene,
	               OrgDb         = orgdb,
	               keyType       = "SYMBOL",
	               ont           = "CC",
	               minGSSize     = minGsize,
	               pAdjustMethod = "BH",
	               pvalueCutoff  = p,
	               qvalueCutoff  = q)
	MF <- enrichGO(gene          = unigene,
	               OrgDb         = orgdb,
	               keyType       = "SYMBOL",
	               ont           = "MF",
	               minGSSize     = minGsize,
	               pAdjustMethod = "BH",
	               pvalueCutoff  = p,
	               qvalueCutoff  = q)
	#请注意:上面BP和CC和MF代码,王铎霖在2018.4.23进行过修改,把参数keytype改成keyType.
	#transfer format
	BP_result <- as.data.frame(BP)
	CC_result <- as.data.frame(CC)
	MF_result <- as.data.frame(MF)
	
	a<-list("BP"=BP, "CC"=CC, "MF"=MF, "BP_result"=BP_result,"CC_result"=CC_result,"MF_result"=MF_result)

	if(outputworkdic==TRUE){
		# write.csv(BP_result,file='GO_BP_result.csv',row.names=FALSE)
		# write.csv(CC_result,file='GO_CC_result.csv',row.names=FALSE)
		# write.csv(MF_result,file='GO_MF_result.csv',row.names=FALSE)
		openxlsx::write.xlsx(x		= a["BP_result"], 
							file 	= paste0(comp_name, filename),
							asTable = T)
		openxlsx::write.xlsx(x		= a["CC_result"], 
							file 	= paste0(comp_name, filename),
							asTable = T)
		openxlsx::write.xlsx(x		= a["MF_result"], 
							file 	= paste0(comp_name, filename),
							asTable = T)
		}
	else if(outputworkdic==FALSE){
		# write.csv(BP_result,file=paste(GO_result_path,paste0(comp_name, "_GO_BP_result.csv"),sep='\\'),row.names=FALSE,quote=FALSE)
		# write.csv(CC_result,file=paste(GO_result_path,paste0(comp_name, "_GO_CC_result.csv"),sep='\\'),row.names=FALSE)
		# write.csv(MF_result,file=paste(GO_result_path,paste0(comp_name, "_GO_MF_result.csv"),sep='\\'),row.names=FALSE)
		openxlsx::write.xlsx(x		= a["BP_result"], 
							file 	= paste(GO_result_path,paste0(comp_name, filename), 
							sep 	= "\\"), 
							asTable = T)
		openxlsx::write.xlsx(x		= a["CC_result"], 
							file 	= paste(GO_result_path,paste0(comp_name, filename), 
							sep 	= "\\"), 
							asTable = T)
		openxlsx::write.xlsx(x		= a["MF_result"], 
							file 	= paste(GO_result_path,paste0(comp_name, filename), 
							sep 	= "\\"), 
							asTable = T)
	}
	else if(outputworkdic==NA){
		#如果outputworkdic==NA,则什么也不做
		}
	return(a)

}