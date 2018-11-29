#data参数接受KEGGresult
KEGG_Barplot <- function(data     = KEGGresult$KEGG_result,
                         plotname = 'KEGG_Barplot.pdf',
                         outputworkdic=TRUE
                         ) {
  stopifnot(is.logical(outputworkdic))
  #data preprocessing
	KEGG_rank <- data[order(data[, "p.adjust"]), ]
	if(nrow(KEGG_rank)>20){
		KEGG_PlotData <- KEGG_rank[1:20, c("Description", "Count")] }else{
		KEGG_PlotData <- KEGG_rank[, c("Description", "Count")] }
		
	PlotData <- as.numeric(KEGG_PlotData[, "Count"])
	name <- as.character(KEGG_PlotData[, "Description"])
	#plot
	if(outputworkdic==TRUE){
		if(nrow(KEGG_PlotData) >= 7){
	  pdf(plotname, width = 0.444*nrow(KEGG_PlotData)+1.75, height = 0.0437*max(nchar(name))+3.726)
	  	if((0.0367*nchar(name[1])+0.11)>=1.25){
	    par(mai=c(0.0437*max(nchar(name))+0.226,0.0367*nchar(name[1]),1,0.5))
	  	}else{
	    par(mai=c(0.0437*max(nchar(name))+0.226,1.15,1,0.6))
	  }
	  xbar <- barplot(PlotData, width = rep(0.7,nrow(KEGG_PlotData)),space =0.7, las = 2, col =rep("firebrick1", times = nrow(KEGG_PlotData)),
	                  border = NA, main = 'Gene Function Classification (KEGG)', ylab = 'Numbers of genes',xpd = T, axisnames = T,
	                  cex.main=1, cex.axis=0.7, ylim=c(0,max(PlotData)))
	  text(x=xbar, y=-0.038*max(PlotData), labels = name, srt = 50, adj = 1, cex = 0.7, xpd = T)
	  dev.off()
		}else{
	  pdf(plotname, width = 5, height = 0.0437*max(nchar(name))+3.726)
	  par(mai=c(0.0437*max(nchar(name))+0.226,2.5-0.222*nrow(KEGG_PlotData),1,2.5-0.222*nrow(KEGG_PlotData)))
	  xbar <- barplot(PlotData, width = rep(0.7,nrow(KEGG_PlotData)),space =0.7, las = 2, col =rep("firebrick1", times = nrow(KEGG_PlotData)),
	                  border = NA, main = 'Gene Function Classification (KEGG)', ylab = 'Numbers of genes',xpd = T, axisnames = T,
	                  cex.main=1, cex.axis=0.7, ylim=c(0,max(PlotData)))
	  text(x=xbar, y=-0.038*max(PlotData), labels = name, srt = 50, adj = 1, cex = 0.7, xpd = T)
	  dev.off()
		}
	}
	else if(outputworkdic==FALSE){
		if(nrow(KEGG_PlotData) >= 7){
	  pdf(paste(KEGG_Barplot_path,paste0(comp_name, '_',plotname), sep = "\\"), width = 0.444*nrow(KEGG_PlotData)+1.75, height = 0.0437*max(nchar(name))+3.726)
	  	if((0.0367*nchar(name[1])+0.11)>=1.25){
	    par(mai=c(0.0437*max(nchar(name))+0.226,0.0367*nchar(name[1]),1,0.5))
	  	}else{
	    par(mai=c(0.0437*max(nchar(name))+0.226,1.15,1,0.6))
	  }
	  xbar <- barplot(PlotData, width = rep(0.7,nrow(KEGG_PlotData)),space =0.7, las = 2, col =rep("firebrick1", times = nrow(KEGG_PlotData)),
	                  border = NA, main = 'Gene Function Classification (KEGG)', ylab = 'Numbers of genes',xpd = T, axisnames = T,
	                  cex.main=1, cex.axis=0.7, ylim=c(0,max(PlotData)))
	  text(x=xbar, y=-0.038*max(PlotData), labels = name, srt = 50, adj = 1, cex = 0.7, xpd = T)
	  dev.off()
		}else{
	  pdf(paste(KEGG_Barplot_path,paste0(comp_name, '_',plotname), sep = "\\"), width = 5, height = 0.0437*max(nchar(name))+3.726)
	  par(mai=c(0.0437*max(nchar(name))+0.226,2.5-0.222*nrow(KEGG_PlotData),1,2.5-0.222*nrow(KEGG_PlotData)))
	  xbar <- barplot(PlotData, width = rep(0.7,nrow(KEGG_PlotData)),space =0.7, las = 2, col =rep("firebrick1", times = nrow(KEGG_PlotData)),
	                  border = NA, main = 'Gene Function Classification (KEGG)', ylab = 'Numbers of genes',xpd = T, axisnames = T,
	                  cex.main=1, cex.axis=0.7, ylim=c(0,max(PlotData)))
	  text(x=xbar, y=-0.038*max(PlotData), labels = name, srt = 50, adj = 1, cex = 0.7, xpd = T)
	  dev.off()
		}
	
	}
	
	}