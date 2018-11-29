volcano_plot<-function(data 			= allData,
						filename		= "coding_Volcanoplot.pdf",
						outputworkdic	= FALSE,
						FC 				= 2,
						pvalue 			=0.05) {
	stopifnot(is.logical(outputworkdic))
	#-----------------------------------------数据处理--------------------------------------------------
	volcano_coding <- allData[,c("Fold Change","FDR P-val","P-val")]
	loc_coding_negFC <- which(volcano_coding[,"Fold Change"]<0)
	volcano_coding[loc_coding_negFC,"Fold Change"] <- -1/(volcano_coding[loc_coding_negFC,"Fold Change"])
	log2FoldChange <- log2(volcano_coding[,"Fold Change"])
	volcano_coding <- as.data.frame(cbind(log2FoldChange,volcano_coding[,c("FDR P-val","P-val")]))
	log2FC <- log2(FC)
	volcano_coding$`catalog`  <-  "no"
	x1 <- volcano_coding$`log2FoldChange` >= log2FC
	volcano_coding$catalog[x1] <- "red"
	y1 <- volcano_coding$`log2FoldChange` <= -log2FC
	volcano_coding$catalog[y1] <- "green"
	z1 <- volcano_coding$`P-val` >= pvalue
	volcano_coding$catalog[z1] <- "no"
	volcano_coding <- na.omit(volcano_coding)
	volcano_coding$`log10_pvalue` <- log10(volcano_coding$`P-val`)
	loc_no <- which(volcano_coding$`catalog`=="no")
	coding_grey <- volcano_coding[loc_no,]
	coding_volplot <- volcano_coding[-loc_no,]
	#-------------------------------------------完成数据处理----------------------------------------------
	#-------------------------------------------开始画图--------------------------------------------------
	library("ggplot2")
	if(outputworkdic==FALSE){
		pdf(paste(diff_gene_path,filename,sep="\\"))

	}
	else if(outputworkdic==TRUE){
		pdf(filename)
	}
	ggplot(coding_volplot)  +
	  	geom_point(aes(log2FoldChange, -log10_pvalue, color = catalog)) +
	  	scale_color_manual(values = c(green="green",red="red"),
	    labels = c(paste("Down:", table(volcano_coding$catalog == "green")[2]),
	      			paste("Up:", table(volcano_coding$catalog == "red")[2])),
	    guide = guide_legend( title = paste( "mRNA:", nrow(volcano_coding) ) ))+
	  	geom_point(data=coding_grey,aes(x=log2FoldChange,y=-log10_pvalue),colour="grey",size=1)
	dev.off()
	#-------------------------------------------完成画图---------------------------------------------------
}