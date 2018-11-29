#data is a numeric data frame of the DEGs' expression.
#samplename is a character vector of sample names.
Cluster_Plot <- function(data         = NULL, 
						 plotname     = "Unsupervised cluster analysis",
						 samplename   = NULL,
						 outputworkdic=TRUE
						 ) {
library('gplots')

stopifnot(is.logical(outputworkdic))

	unigene <- unique(data[which(data[,"Gene Symbol"]!= "NA" & data[,"Gene Symbol"]!= ""),"Gene Symbol"])
	s_exp <- NULL
	for(i in 1:length(unigene)){
	probes <- data[unigene[i] == data$`Gene Symbol`, samplename]
	mean <- colMeans(probes)
	s_exp <- rbind(s_exp, mean)
	rownames(s_exp)[i] <- unigene[i]
	}
if(outputworkdic==TRUE){
	pdf(paste0(plotname,'(无基因名).pdf'), height = 8, width = 8)
	par(mar=c(12, 10, 7, 12),oma=c(3,3,3,5),no.readonly=T)
	heatmap.2(s_exp,col = colorRampPalette(c("green","black","red"))(100),dendrogram="both", labRow = FALSE,Colv="Rowv",scale = "row",
			  density.info = "none", trace = "none", cexCol = 0.8,cexRow = 0.8,
			  lhei = c(2,8), lwid = c(2,5))
	dev.off()
	pdf(paste0(plotname,'(有基因名).pdf'), height = 1.6+0.16*nrow(s_exp), width = 10)
	par(mar=c(12, 10, 7, 12),oma=c(2,10,1,10),pin=c(3,0.5),no.readonly=T,pty='m')
	heatmap.2(s_exp,col = colorRampPalette(c("green","black","red"))(100),dendrogram="both",Colv="Rowv",scale = "row",#Colv是用来画列线的
			  density.info = "none", trace = "none", cexCol = 0.8,cexRow = 0.8,margins = c(2,5),
			  lhei = c(10,nrow(s_exp)), lwid = c(2,5))
	dev.off()
	}

else if(outputworkdic==FALSE){
	pdf(paste(diff_gene_path, paste0(comp_name,
		plotname,'(无基因名).pdf'), sep = "\\"), height = 8, width = 8)
	par(mar=c(12, 10, 7, 12),oma=c(3,3,3,5),no.readonly=T)
	heatmap.2(s_exp,col = colorRampPalette(c("green","black","red"))(100),dendrogram="both", labRow = FALSE,Colv="Rowv",scale = "row",
			  density.info = "none", trace = "none", cexCol = 0.8,cexRow = 0.8,
			  lhei = c(2,8), lwid = c(2,5))
	dev.off()
	pdf(paste(diff_gene_path, paste0(comp_name,
		plotname,'(有基因名).pdf'), sep = "\\"), height = 1.6+0.16*nrow(s_exp), width = 10)
	par(mar=c(12, 10, 7, 12),oma=c(2,10,1,10),pin=c(3,0.5),no.readonly=T,pty='m')
	heatmap.2(s_exp,col = colorRampPalette(c("green","black","red"))(100),dendrogram="both",Colv="Rowv",scale = "row",#Colv是用来画列线的
			  density.info = "none", trace = "none", cexCol = 0.8,cexRow = 0.8,margins = c(2,5),
			  lhei = c(10,nrow(s_exp)), lwid = c(2,5))
	dev.off()
	}
}
