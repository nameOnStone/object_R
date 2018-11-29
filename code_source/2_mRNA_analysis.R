mRNA_analysis <- function(  filename1 = file_Name_Gene, 
							filename2 = file_Name_ID2ID,
							comp_name = NULL, 
							organism  = NULL,
							FC 		  = NULL,
							pvalue 	  = NULL) {
# 此函数可接受参数filename1=charcter或data.frame
	library("openxlsx")
	library("gplots")
	library("pathview")
	library("clusterProfiler")
	library('corrplot')
	library('tidyverse')

	org_info <- organismMapper(organism)
	# ======================================判断filename1和filename2的数据类型====================================
	if (class(filename1) == "character") {
	# =====================================判断是否存在以及大小============================
		if (file.exists(filename1) & file.info(filename1)[filename1, "size"] >1000) {
			allData <- read.table(filename1, header = TRUE, sep = "\t", quote = "",
				dec = ".", na.strings = "---", row.names = NULL, skip = 0, check.names = FALSE,
				blank.lines.skip = TRUE, strip.white = TRUE, comment.char = "#",
				fill = TRUE, as.is = TRUE)
		} else {
		stop('gene_data is not exist or size of gene_data is wrong!')
		}
	} 
	else if (class(filename1) == "data.frame") {
	allData <- filename1
	} 
	if (file.exists(filename2) & file.info(filename2)[filename2, "size"] < 1000) {
		ID2ID <- read.table(filename2, header = F, sep = "\t",stringsAsFactors = F)
	} else {
	stop("ChipIDSampleID.txt is not exist or size of ChipIDSampleID.txt is wrong")
	}
	# ====================================换名========================================
	for (i in 1:dim(ID2ID)[1]) {
	colnames(allData) <- sub(ID2ID[i, 1], ID2ID[i, 2], colnames(allData))
	}
	loc <- grep("Signal", colnames(allData))
	colnames(allData)[loc] <- gsub("\\_.*", "_Signal", colnames(allData)[loc])
	# ====================================完成换名========================================
#-------------------------------------以下代码是画箱线图---------------------------------------
	ID2ID1<-ID2ID[order(ID2ID[,3]),]
	GroupNUM<-length(unique(ID2ID[,3]))#group numbers
	freq<-as.data.frame(table(ID2ID[,3]))
	NUM<-freq$'Freq'#sample numbers for each group
	Order<-as.character(ID2ID1[,2])
	#注意上一行Order是经过排序后的样品名称,而且是全部的。
	QCplot_matrix <- allData[, loc]
	#
	aaa<-c()
	for (i in Order){
		bbb<-str_subset(string = colnames(QCplot_matrix),pattern = i)
		aaa<-c(aaa,bbb)

	}
	QCplot_matrix<-QCplot_matrix[,aaa]
	#
	fff<-c()
	for (i in Order){
		ggg<-str_detect(string=colnames(QCplot_matrix),pattern = i)
		if(TRUE%in%ggg){
			fff<-c(fff,i)
		}
	}
	Order_handled<-fff

	color <- c("firebrick1","dodger blue","gold","lightskyblue","yellow","deepskyblue",
			"cornflowerblue","lightpink","steelblue","#FD8338","#DED847","#E37CF8",
			"#FC6672")
	Color <- rep(color[1:GroupNUM],each=NUM)
	pdf(paste(exploratory_BoxWhisker_path,'BoxWhisker.pdf',sep='\\'),family="GB1")
		par(mar=c(6, 5, 3, 3))
		boxplot(x=QCplot_matrix, col=Color, font.lab=2, cex.lab=1.5,pars = list(boxwex = 0.6, 
				staplewex = 0.5, outwex = 0.5), ylab="Signal", las='2', cex.axis=0.8)
	dev.off()
#-------------------------------------完成画箱线图---------------------------------------
#-------------------------------------以下代码是画相关系数图---------------------------------------
	Pearson <- cor(QCplot_matrix,method = "pearson",use="complete.obs")
	pdf(paste(exploratory_Peason_path,'Pearson.pdf',sep='\\'),family="GB1")
	par(mar=c(5, 5, 0.5, 1))
	corrplot(corr = Pearson, method = 'color',order = "AOE",tl.col = "black", 
			cl.ratio=0.1, cl.lim=c(0,1),cl.length = 5,cl.align.text = "c",
			cl.cex = 0.8,addCoef.col="grey")
	dev.off()
	write.csv(Pearson, paste(exploratory_Peason_path,"cor-pearson.csv",sep='\\'))
#-------------------------------------完成画相关系数图---------------------------------------
#-------------------------------------以下是画主成分分析图---------------------------------------------
	#计算主成分PCA
	tSignal <-t(QCplot_matrix)
	PCA_exp <- prcomp(tSignal,scale=F)
	SummaryPCA <- summary(PCA_exp)
	
	#画图数据整理
	PC <- round(PCA_exp$x,1 )

	group <- as.matrix(ID2ID1[match(Order_handled,ID2ID1[,2]),3])
	PCA_plot <- cbind(PC[,1:2], group, Order_handled) 
	#PCA_plot <- cbind(PC[,1:2], Group)
	colnames(PCA_plot)[3] <- "Group"
	rownames(PCA_plot) <- NULL
	PCA_plot <- PCA_plot[order(PCA_plot[,"Order_handled"]),]
	PCA_plot <- data.frame(PCA_plot)
	percentVar <- round(100*SummaryPCA$importance[2,],digit=1)
	
	#library(ggplot2)
	ggplot(data = PCA_plot, aes(x =PC1, y =PC2,color =Group))+
    geom_point(size = 1.8)+
    xlab(paste0("PC1 : ", percentVar[1], "% variance")) +
    ylab(paste0("PC2 : ", percentVar[2], "% variance"))+ 
    coord_fixed(ratio = 1/2)+
    theme(panel.grid.major = element_line(colour = 'grey50'))+
    scale_color_manual(values = c('purple ','orange'))
	#p + geom_text(hjust = 0, nudge_x = 1)#添加样本名
	ggsave(paste(exploratory_PCA_path,"PCA.pdf",sep='\\'))
#-------------------------------------完成主成分分析图---------------------------------------------
#-------------------------------------以下是画火山图----------------------------------------------------------
	#数据处理
	#取值log2(Fold Change)和-log10(FDR)用于画图
	volcano_coding <- allData[,c("Fold Change","FDR P-val","P-val")]
	aa1 <- which(volcano_coding[,"Fold Change"]<0)
	volcano_coding[aa1,"Fold Change"] <- -1/(volcano_coding[aa1,"Fold Change"])
	log2FoldChange <- log2(volcano_coding[,"Fold Change"])
	volcano_coding <- as.data.frame(cbind(log2FoldChange,volcano_coding[,c("FDR P-val","P-val")]))
	log2FC <- log2(FC)
	volcano_coding$catalog  <-  'Nonsignificant'
	x1 <- volcano_coding$log2FoldChange > log2FC
	volcano_coding$catalog[x1] <- "UP"
	y1 <- volcano_coding$log2FoldChange < -log2FC
	volcano_coding$catalog[y1] <- "Down"
	z1 <- volcano_coding$`P-val` > pvalue
	volcano_coding$catalog[z1] <- "Nonsignificant"
	volcano_coding <- na.omit(volcano_coding)
	volcano_coding$`log10_p_value` <- log10(volcano_coding$`P-val`)
	#library("ggplot2")
	#开始画图
	# pdf(paste(,"coding_volcanoplot.pdf",sep='\\'),family="GB1")
	# ggplot(data=volcano_coding)  +
		# geom_point(aes(log2FoldChange, -log10_p_value, color = catalog)) +
		 # scale_color_manual(values = c(green="green",red="red",grey='grey'),
		# labels = c(	paste("Down:", table(volcano_coding$catalog == "green")[2]),
					# paste('Nonsignificant:',table(volcano_coding$catalog == 'grey')[2]),
					# paste("Up:", table(volcano_coding$catalog == "red")[2])),
		# guide = guide_legend( title = paste( "mRNA:", nrow(volcano_coding))))
	# ggsave(paste(diff_gene_path,"coding_volcanoplot.pdf",sep='\\'))

	ggplot(data=volcano_coding,aes(x=log2FoldChange,y=-log10_p_value,colour=catalog))+
	geom_point(size=0.6)+
	scale_color_manual(values = c(Down='darkgreen',Nonsignificant='grey',UP='darkred'),
		labels = c(	paste("Down :", table(volcano_coding$catalog == "Down")[2]),
					paste('Nonsignificant :',table(volcano_coding$catalog == 'Nonsignificant')[2]),
					paste("Up :", table(volcano_coding$catalog == "UP")[2])),
		guide = guide_legend( title = paste( "Total :", nrow(volcano_coding))))+
	labs(x='log2(FoldChange)',y='-log10(pvalue)')+
	# theme(axis.title.y=element_text(size=10,vjust=0.5))+
	geom_vline(data=volcano_coding,mapping=aes(xintercept=log2FC))+
	geom_vline(data=volcano_coding,mapping=aes(xintercept=-log2FC))+
	geom_hline(data = volcano_coding,mapping = aes(yintercept=-log10(pvalue)))+
	theme_light(base_size=8)
	
	ggsave(paste(diff_gene_path,"coding_volcanoplot.pdf",sep='\\'))
#----------------------------------------完成火山图---------------------------------------------------
	# samplename等于下面的,请思考一下.
	samplename <- colnames(allData[loc])

	diff_result <<- screen_DEGs(	data 	= allData,
								filename 	='_diff_gene_unfilter.csv',
								FC 			=FC,
								p 			=pvalue,
							outputworkdic 	= FALSE)

	Cluster_Plot(			data = diff_result, 
						plotname = "_Unsupervised cluster analysis", 
					samplename 	 = samplename,
				outputworkdic	 = FALSE)  

	GOresult <<- GO_Analysis(    data  	= diff_result, 
							filename 	= "_GO_result.xlsx", 
								orgdb 	= org_info$orgdb,
							minGsize 	= 1, 
								p 		= 1, 
								q 		= 1,
						outputworkdic 	= FALSE) 
	
	total <- rbind(GOresult$BP_result, GOresult$CC_result, GOresult$MF_result)
	if (nrow(total) != 0) {
		GO_Barplot( data 		= GOresult, 
					plotname 	= "_Barplot(GO).pdf",
				outputworkdic	= FALSE)

		if(nrow(GOresult$BP_result)!=0){
			if(nrow(GOresult$BP_result)>=10){
			GO_Plot(	data 	= GOresult$BP, 
					plotname	= "_BP.pdf",
	numsofshowingsignificant	= 10,
				outputworkdic	= FALSE)
			}
			else{
			GO_Plot(	data 	= GOresult$BP, 
					plotname	= "_BP.pdf",
	numsofshowingsignificant	= nrow(GOresult$BP_result),
				outputworkdic	= FALSE)
			}
		}

		if(nrow(GOresult$CC_result)!=0){
			if(nrow(GOresult$CC_result)>=10){
			GO_Plot(	data 	= GOresult$CC, 
					plotname	= "_CC.pdf",
	numsofshowingsignificant	= 10,
				outputworkdic	= FALSE)
			}
			else{
			GO_Plot(	data 	= GOresult$CC, 
					plotname	= "_CC.pdf",
	numsofshowingsignificant	= nrow(GOresult$CC_result),
				outputworkdic	= FALSE)
			}
		}

		if(nrow(GOresult$MF_result)!=0){
			if(nrow(GOresult$MF_result)>=10){
			GO_Plot(	data 	= GOresult$MF, 
					plotname	= "_MF.pdf",
	numsofshowingsignificant	= 10,
				outputworkdic	= FALSE)
			}
			else{
			GO_Plot(	data 	= GOresult$MF, 
					plotname	= "_MF.pdf",
	numsofshowingsignificant	= nrow(GOresult$MF_result),
				outputworkdic	= FALSE)
			}
		}
	} else if (total == 0) {
		stop("GO_Analysis is no result,find the reason~~")
	}

	KEGGresult <<- KEGG_Analysis(data 	= diff_result, 
							filename 	= 'KEGGresult.csv', 
								orgdb 	= org_info$orgdb, 
							species 	= org_info$species,
							minGsize 	= 1, 
								p 		= 1, 
								q 		= 1, 
						outputworkdic 	= FALSE)

	if (nrow(KEGGresult$KEGG_result) != 0) {
		# KEGG Barplot  
		KEGG_Barplot(data 		= KEGGresult$KEGG_result, 
					plotname 	= 'KEGG_Barplot.pdf',
				outputworkdic 	= FALSE)
		# KEGG Dotplot
		KEGG_Dotplot(data 		= KEGGresult$KEGG_result, 
					plotname 	= 'Dotplot.pdf',
				outputworkdic	= FALSE)
		# KEGG plot 
		KEGG_Plot(	data 		= diff_result, 
					KEGG 		= KEGGresult$KEGG, 
					orgdb 		= org_info$orgdb,
					species 	= org_info$species, 
				outputworkdic 	= FALSE)
	} else if (KEGGresult$KEGGresult == 0) {
		stop("KEGG_Analysis is no result~~")
	}
}