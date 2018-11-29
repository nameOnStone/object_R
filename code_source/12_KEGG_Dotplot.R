#此函数的参数data只接受data.frame数据格式,请注意
#如果outputworkdic=FASLE,则不必给出plotname参数
KEGG_Dotplot <- function(	data     = KEGGresult$KEGG_result,
                    		plotname = 'Dotplot.pdf',
                    		outputworkdic=TRUE
                    ) {

	library(ggplot2)
	stopifnot(is.logical(outputworkdic))
	##data preprocessing
	GeneRatio <- data$`GeneRatio`
	a <- as.numeric(unlist(strsplit(GeneRatio, split = "/")))
	gene_numbers <- matrix(data = a, ncol = 2, byrow = TRUE)
	data$`GeneRatio` <- gene_numbers[,1]/gene_numbers[,2]
	KEGG_rank <- data[order(data$p.adjust), ]
	if(nrow(KEGG_rank)>20){
		PlotData <- KEGG_rank[1:20, ] 
		}
	else{PlotData <- KEGG_rank
		 }
	b <- PlotData[order(PlotData$p.adjust, decreasing = TRUE),]
	picture <- ggplot(b)+geom_point(aes(x=b$GeneRatio,y=b$Description,colour=b$p.adjust,size=b$Count))+
		scale_colour_gradientn(colours = c("red","green"),limit=c(0,max(b$p.adjust)))+
		labs(colour="p.adjust",size="Gene_number",x="GeneRatio", y="")+
		theme_bw()+
		ggtitle("Statistics of Pathway Enrichment")
	if(outputworkdic==TRUE){
	ggsave(filename=plotname,plot=picture,width = 10, height = 7,device=NULL)
	}
	else if(outputworkdic==FALSE){
	ggsave(filename=paste(KEGG_Barplot_path,paste0(comp_name, "_",plotname), sep = "\\"),plot=picture,width = 10, height = 7,device=NULL)
	}
}