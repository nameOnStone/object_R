
KEGG_Plot <- function(data    =diff_result,
					  KEGG    = KEGGresult$`KEGG`,
					  orgdb   = "org.Hs.eg.db",
					  species = "hsa",
					  outputworkdic=FALSE) {

	# library("RDAVIDWebService")
	library("DOSE")
	library("clusterProfiler")
	library("pathview")
	stopifnot(is.logical(outputworkdic))

	if(outputworkdic==FALSE){
		a<-getwd()
		setwd(KEGG_pathway_path)
	}else if(outputworkdic==TRUE){
		a<-getwd()
	}
	#convert the format of gene ID
	unigene <- unique(data[which(data[,"Gene Symbol"]!= "NA" & data[,"Gene Symbol"]!= ""),"Gene Symbol"])
	loc1 <- match(unigene, data$`Gene Symbol`)
	diff_gene <- data[loc1, c("Gene Symbol", "log2(Fold Change)")]
	
	geneID <- bitr(gene = diff_gene$`Gene Symbol`,
								 fromType = "SYMBOL",
								 toType = "ENTREZID",
								 OrgDb = orgdb)
	
	loc2 <- match(geneID$`SYMBOL`, diff_gene$`Gene Symbol`)
	FC <- as.matrix(diff_gene[loc2, "log2(Fold Change)"])
	rownames(FC) <- geneID$`ENTREZID`
	
	#plot
		allKEGG <- viewKEGG1(obj = KEGG,
							 pathwayID = "all",
							 foldChange = FC,
							 species = species,
							 color.low="green",
							 color.high="red",
							 kegg.native=TRUE,
							 out.suffix="pathway")
	
	plotnames <- dir()
	loc <- grep('.pathway.png',plotnames)
	removeplot <- plotnames[-loc]
	file.remove(removeplot)
	setwd(a)
}
	