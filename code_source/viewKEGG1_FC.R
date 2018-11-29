#所有KEGG的通路图的绘制，改自clusterProfiler包中viewKEGG
viewKEGG1 <- function(obj, pathwayID, foldChange,
					 color.low="green",
					 color.high="red",
					 species = NULL,
					 kegg.native=TRUE,
					 out.suffix="clusterProfiler") {
  
  if (class(obj) != "enrichResult")
	stop("only enrichResult object supported.")
  if (obj@ontology != "KEGG")
	stop("only KEGG supported.")
  
  print("viewKEGG is a wrapper function of pathview")
  citation("pathview")
  
  pkg <- "pathview"
  suppressMessages(require(pkg, character.only=TRUE))
  if (is.numeric(pathwayID)) {
	pathwayID <- as.data.frame(obj)[pathwayID, 1]
  }
  if (length(pathwayID) == 1 & pathwayID == "all") {
	pathwayID <- as.data.frame(obj)[, 1]
  }
	   
  m.fc <- max(abs(foldChange))
  bins <- ceiling(m.fc) * 2
  if (bins < 10)
		bins <- 10
		pathview <- eval(parse(text=pkg))
		res <- lapply(pathwayID, function(pid) {
		pathview(gene.data=foldChange,
				pathway.id = pid,
				species    = species,
				limit = list(gene=m.fc, cpd=1),
				bins = list(gene=bins, cpd=10),
				low = list(gene=color.low, cpd="blue"),
				high = list(gene=color.high, cpd="yellow"),
				kegg.native=kegg.native,
				plot.col.key=TRUE,
				key.pos='topright',
				out.suffix=out.suffix,
				new.signature=FALSE)
   })
  return (res)
 }

