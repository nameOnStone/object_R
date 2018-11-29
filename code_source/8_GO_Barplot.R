GO_Barplot <- function(data     = GOresult,
					   plotname = "Barplot(GO).pdf",
					   outputworkdic=TRUE
					   ) {
  
	stopifnot(is.logical(outputworkdic))
  #data preprocessing
	BP_data <- data$`BP_result`
	CC_data <- data$`CC_result`
	MF_data <- data$`MF_result`
	BP_rank <- BP_data[order(BP_data[,"p.adjust"]),]
	CC_rank <- CC_data[order(CC_data[,"p.adjust"]),]
	MF_rank <- MF_data[order(MF_data[,"p.adjust"]),]
	
	#extract the top 20 data
	if(nrow(BP_rank)>20){
		BP_PlotData <- BP_rank[1:20,c("Description", "Count")] }else{
		BP_PlotData <- BP_rank[,c("Description", "Count")] }
	if(nrow(CC_rank)>20){
		CC_PlotData <- CC_rank[1:20,c("Description", "Count")] }else{
		CC_PlotData <- CC_rank[,c("Description", "Count")] }
	if(nrow(MF_rank)>20){
		MF_PlotData <- MF_rank[1:20,c("Description", "Count")] }else{
		MF_PlotData <- MF_rank[,c("Description", "Count")] }
	
	#rbind	
	PlotData <- rbind(BP_PlotData, CC_PlotData, MF_PlotData)
	Data <- as.numeric(PlotData[,"Count"])
	name <- as.character(PlotData[,"Description"])
	
	#plot
	if(outputworkdic==TRUE){
		pdf(plotname, width = 10.8+0.1365*nrow(PlotData), height = 5.28+0.045*max(nchar(name)))
	}
	else if(outputworkdic==FALSE){
		pdf(paste(GO_Barplot_path, paste0(comp_name,'_',plotname), sep = "\\" ), width = 10.8+0.1365*nrow(PlotData), height = 5.28+0.045*max(nchar(name)))
	}
		if(0.0378*nchar(name[1]) >= 1.2){
		par(mai=c(0.045*max(nchar(name))+0.28,0.0378*nchar(name[1]),1,0.6))
		}else{
		par(mai = c(0.045*max(nchar(name))+0.28,1.2,1,0.6))
		}
		xbar=barplot(Data, space =0.5, las = 2, col =c(rep("firebrick1",nrow(BP_PlotData)),rep("goldenrod1",nrow(CC_PlotData)),rep("deepskyblue2",nrow(MF_PlotData))),
			border = NA, main = 'Gene Function Classification (GO)', ylab = 'Numbers of genes', xpd = T, axisnames = T,
			cex.main=1, cex.axis = 0.7, ylim=c(0, max(Data)+3))
		legend("topright", legend = c("Biological Process","Cellular Component","Molecular Function"), fill=c("firebrick1","goldenrod1","deepskyblue2"),cex = 0.7, border = NA)
		text(x=xbar, y=-(0.02*max(Data)+0.12), labels = name, srt = 50, adj = 1, xpd = T, cex = 0.7)
		dev.off()
}