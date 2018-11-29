BoxWhisker<-function(	data 		 = allData,
						ID 			 = ID2ID,
						filename	 ='BoxWhisker.pdf',
						outputworkdic=TRUE){
	stopifnot(is.logical(outputworkdic))
	ID2ID1<-ID2ID[order(ID2ID[,3]),]
	GroupNUM<-length(unique(ID2ID[,3]))#group numbers
	freq<-as.data.frame(table(ID2ID[,3]))
	NUM<-freq$'Freq'#sample numbers for each group
	Order<-as.character(ID2ID1[,2])

	QCplot_matrix <- allData[, loc]
	QCplot_matrix <- QCplot_matrix[,Order]
	N<-sum(NUM)
	color <- c("firebrick1","dodger blue","gold","lightskyblue","yellow","deepskyblue",
        "cornflowerblue","lightpink","steelblue","#FD8338","#DED847","#E37CF8",
        "#FC6672")
	Color <- rep(color[1:GroupNUM],each=NUM)
	if(outputworkdic==TRUE){
	pdf(filename, family="GB1")
	}
	else if(outputworkdic==FALSE){
	pdf(paste(exploratory_BoxWhisker_path,filename,sep='\\'), family="GB1")#family="Helvetica"
	}
  	par(mar=c(6, 5, 3, 3))
  	boxplot(QCplot_matrix, col=Color, font.lab=2, cex.lab=1.5,pars = list(boxwex = 0.6, 
         	staplewex = 0.5, outwex = 0.5), ylab="Signal", las='2', cex.axis=0.8)
	dev.off()
}

	