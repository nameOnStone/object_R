GO_Plot <- function(data     = NULL,
					plotname = 'GO_Plot.pdf',
	numsofshowingsignificant = NULL,
				outputworkdic=TRUE) {
	stopifnot(is.logical(outputworkdic))
	#loading database
	# library("RDAVIDWebService")
	library("DOSE")
	library("clusterProfiler")
	#plot
	tryCatch({
		if(outputworkdic==TRUE){
			pdf(plotname, family = "Helvetica", height = 10, width = 12)
				}
		else if(outputworkdic==FALSE){
			pdf(paste(GO_Enrishment_path, paste0(comp_name,'_',plotname), sep = "\\"), family = "Helvetica", height = 10, width = 12)
		}
				plotGOgraph(data,firstSigNodes=numsofshowingsignificant)
				dev.off()},error = function(e){cat("ERROR:", conditionMessage(e),"\n")})
	if(length(dev.list()) != 0){
		graphics.off()
	}
}