Pearson_plot<-function(	data 			= allData,
						filename 		= 'Pearson.pdf',
						outputworkdic	= FALSE){
	#数据处理
	stopifnot(is.logical(outputworkdic))
	loc <- grep("Signal", colnames(data))
	QCplot_matrix <- allData[, loc]
	Pearson <- cor(QCplot_matrix,method = "pearson",use="complete.obs")
	#处理完成
	if(outputworkdic==TRUE){
		pdf(filename,family="GB1")
		par(mar=c(5, 5, 0.5, 1))
		corrplot(corr = Pearson, method = 'color',order = "AOE",tl.col = "black", 
			cl.ratio=0.1, cl.lim=c(0,1),cl.length = 5,cl.align.text = "c",
			cl.cex = 0.8,addCoef.col="grey")
		dev.off()
		write.csv(Pearson, "cor-pearson.csv")
	}
	else if(outputworkdic==FALSE){
		pdf(paste(exploratory_Peason_path,'Pearson.pdf',sep='\\'),family="GB1")
		par(mar=c(5, 5, 0.5, 1))
		corrplot(corr = Pearson, method = 'color',order = "AOE",tl.col = "black", 
			cl.ratio=0.1, cl.lim=c(0,1),cl.length = 5,cl.align.text = "c",
			cl.cex = 0.8,addCoef.col="grey")
		dev.off()
		write.csv(Pearson, paste(exploratory_Peason_path,"cor-pearson.csv",sep='\\'))
	}
}