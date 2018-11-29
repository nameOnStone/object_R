
Diffdata <- read.table()
ID2ID <- read.table()

for (i in 1:dim(ID2ID)[1]){
	 colnames(Diffdata) <- sub(ID2ID[i,1], ID2ID[i,2], colnames(Diffdata))
	} 


colnames(Diffdata)<-gsub("\\_.*", "", colnames(Diffdata))
Diffdata <- Diffdata[,-1]

Pearson <- cor(Diffdata, method = "pearson", use="complete.obs")
  write.table(Pearson, "cor-pearson.txt", quote = F)
  #plot
  library(corrplot)
  pdf('Pearson.pdf',family="GB1" ,width = 8, height = 8)
  par(mar=c(10, 10, 2, 1))
  corrplot(corr = Pearson, method = 'color',order = "original",tl.col = "black", 
		   cl.ratio=0.1, cl.lim=c(0,1),cl.length = 5,cl.align.text = "c",cl.cex = 0.8,addCoef.col="grey")
   dev.off()