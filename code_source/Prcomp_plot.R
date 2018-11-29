Prcomp_plot<-function(	data 			= allData,
						filename 		= 'Procomp.pdf',
						outputworkdic	= TRUE){
	stopifnot(is.logical(outputworkdic))
	#数据处理
	loc <- grep("Signal", colnames(data))
	QCplot_matrix <- data[, loc]
	tSignal <-t(QCplot_matrix)
	PCA_exp <- prcomp(tSignal,scale=F)
	SummaryPCA <- summary(PCA_exp)
	PC <- round(PCA_exp$x,1 )
	group <- as.matrix(ID2ID1[ID2ID1[,2]==Order_handled,3])
	PCA_plot <- cbind(PC[,1:2], group, Order_handled) 
	#PCA_plot <- cbind(PC[,1:2], Group)
	colnames(PCA_plot)[3] <- "Group"
	rownames(PCA_plot) <- NULL
	PCA_plot <- PCA_plot[order(PCA_plot[,"Order_handled"]),]
	PCA_plot <- data.frame(PCA_plot)
	percentVar <- round(100*SummaryPCA$importance[2,],digit=1)
	#开始画图
	library(ggplot2)
	ggplot(PCA_plot, aes(x =PC1, y =PC2,color =Group, label =Order_handled))+
		geom_point(size = 2)+
		xlab(paste0("PC1: ", percentVar[1], "% variance")) +
		ylab(paste0("PC2: ", percentVar[2], "% variance"))+ 
		coord_fixed(ratio = 1/2)
	#p + geom_text(hjust = 0, nudge_x = 1)#添加样本名
	if(outputworkdic==TRUE){
		ggsave(filename)
	}
	else if(outputworkdic==FALSE){
		ggsave(paste(exploratory_PCA_path,"PCA.pdf",sep='\\'))
	}
}