signals <- read.table("sample_signals.txt", header = TRUE,
  comment.char = '#', sep = "\t", check.names = FALSE)

colnames(signals) <- gsub("_\\(.*Signal", "", colnames(signals))
rownames(signals) <- signals$ID
signals <- signals[, -1]

controls <- read.table("tac_unique.txt")
controls <- controls$V1

drop <- rownames(signals) %in% controls
signals <- signals[!drop, ]

# copy and modified from miRNA 104 code
library(corrplot)

Data <- signals
Pearson = cor(Data, method = "pearson",use="complete.obs")
svg('Pearson.svg',family="GB1" ,width = 8, height = 8)
par(mar=c(5, 5, 0.5, 1))
   corrplot(corr = Pearson, method = 'color',tl.col = "black", 
           cl.ratio=0.1, cl.lim=c(0,1),cl.length = 10,cl.align.text = "c",cl.cex = 0.8,
		   addCoef.col="grey", order = "alphabet")
dev.off()

# pick up DEG
deg <- read.table("DEG_p0.05.txt", head = TRUE,
  skip = 4, sep = "\t")
deg_id <- deg$ID

deg_idx <- rownames(signals) %in% deg_id
deg_signals <- signals[deg_idx, ]

# copy and modified from miRNA 203 code
library(gplots)

s_exp <- as.matrix(deg_signals)


svg("heatmap.svg", family="GB1", height = 8, width = 8)
par(mar=c(10, 8, 5, 10),oma=c(1,1,1,4.5))
heatmap.2(s_exp,
  col = colorRampPalette(c("green", "black", "red")), 
  scale = "row",
  labRow = FALSE,
  density.info = "none",
  trace = "none", 
  cexCol = 1.2,
  cexRow = 1.1,
  lwid = c(2,5),
  Colv = FALSE,
  dendrogram = "row")
dev.off()

# DEPRECATED
# copy and modified from de novo 12.R code
library(pheatmap)

# HEIGHT
DElogcpm <- signals

nrow(DElogcpm)/68

# pheatmap(b, cellheight = 5, scale = "row", show_rownames = FALSE)
svg(paste("et_1_", round(nrow(DElogcpm)/68), ".svg", sep = ""), 
    height = nrow(DElogcpm)/68, onefile = FALSE)
pheatmap(DElogcpm, cellheight = 1, scale = "row", show_rownames = FALSE)
dev.off()

pdf(paste("et_1_", round(nrow(DElogcpm)/68), "_annotate.pdf", sep = ""),
    height = nrow(DElogcpm)/68, onefile = FALSE)
pheatmap(DElogcpm, cellheight = 1, scale = "row", show_rownames = TRUE, fontsize = 1)
dev.off()

# 10520 ~ 22000
# 17995 ~ 40000
# 15511 ~ 32000

# HEIGHT
2.14*nrow(DElogcpm)

png(paste("et_2_", round(2.14*nrow(DElogcpm)), ".png", sep = ""), height = 2.14*nrow(DElogcpm))
pheatmap(DElogcpm, cellheight = 2, scale = "row", show_rownames = FALSE)
dev.off()

# DEPRECATED ENDS HERE

# copy and modified from miRNA4.0 204 code
MatureData <- read.table("mouse_mature_result.txt", header = TRUE, sep = "\t",
  check.names = FALSE)
BaseData <- read.table("mouse_targets_result.txt", header = TRUE, sep = "\t",
  check.names = FALSE)


miRNA_data <- deg[, c("ID", "Sequence")]
#   miRNA_data <- unique(Diffdata[,c("Transcript ID(Array Design)", "Sequence")])
   
  #---------------------???ҳ???miRNA------------------------------------
#   loc1 <- match(miRNA_data$`Sequence`, MatureData$`Sequence`)
#   Diff_mature <- na.omit(MatureData[loc1,])
   
   mature_idx <- MatureData$Sequence %in% miRNA_data$Sequence
   Diff_mature <- MatureData[mature_idx, ]


#   loc2 <- match(BaseData$`MiRBase Accession`, Diff_mature$`miRBase Accession`)
#   loc3 <- !is.na(loc2)
#   targetgene <- BaseData[loc3,]

   db_idx <- BaseData$`MiRBase Accession` %in% Diff_mature$`miRBase Accession`
   targetgene <- BaseData[db_idx, ]


#   setwd(path)
  #--------------------?ֱ???ȡ??ͬ???ݿ??ал???-------------------------------------				 
  if(species=="hsa"){
    A=targetgene[targetgene$`miRecords`==1, c("MiRBase ID", "Gene Symbol")]
    A1=unique(paste(A$`MiRBase ID`, A$`Gene Symbol`, sep="_"))

    B=targetgene[targetgene$`miRTarbase`==1, c("MiRBase ID", "Gene Symbol")]
    B1=unique(paste(B$`MiRBase ID`, B$`Gene Symbol`, sep="_"))

    C=targetgene[targetgene$`TargetScan`==1, c("MiRBase ID", "Gene Symbol")]
    C1=unique(paste(C$`MiRBase ID`, C$`Gene Symbol`, sep="_"))

    D=targetgene[targetgene$`miRanda`==1, c("MiRBase ID", "Gene Symbol")]
    D1=unique(paste(D$`MiRBase ID`, D$`Gene Symbol`, sep="_"))

    E=targetgene[targetgene$`miRDB`==1, c("MiRBase ID", "Gene Symbol")]
    E1=unique(paste(E$`MiRBase ID`, E$`Gene Symbol`, sep="_"))
    
     #----------------------------------??Τ??ͼ----------------------------------------------------------------------------------------
   library(VennDiagram)

  VennName = "venn.tiff"

    venn.plot <- venn.diagram(x=list(A=A1, B=B1,C=C1,D=D1,E=E1),filename = VennName, resolution =500, imagetype = "tiff",col="black",
	                 fill = c("red","blue3","green","yellow","purple"), alpha = 0.5,fontfamily = "serif",fontface = "bold",cat.default.pos = "text",
				     cex = 1.25,cat.col =c("darkred", "darkblue", "darkgreen","olivedrab","purple3"),cat.cex = 1.25,main = "A=miRecords  B=miRTarbase  C=TargetScan D=miRanda  E=miRDB",
				     main.pos= c(0.5, 0.23), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 1.1, 
				     main.just = c(0.5, 0.8),cat.fontfamily = "serif", cat.dist = c(0.05,0.05,-0.05,-0.05,0.05),cat.pos = 0,cat.fontface = "bold",margin = 0.47)

	 x <- dir()
	 y <- x[grep('\\.log.*', x)]
	 file.remove(y)
	 
	}else if(species=="mmu"){
	   A=targetgene[targetgene$`miRanda`==1, c("MiRBase ID", "Gene Symbol")]
	   A1=unique(paste(A$`MiRBase ID`, A$`Gene Symbol`, sep="_"))
	   B=targetgene[targetgene$`miRDB`==1, c("MiRBase ID", "Gene Symbol")]
	   B1=unique(paste(B$`MiRBase ID`, B$`Gene Symbol`, sep="_"))
	   C=targetgene[targetgene$`TargetScan`==1, c("MiRBase ID", "Gene Symbol")]
	   C1=unique(paste(C$`MiRBase ID`, C$`Gene Symbol`, sep="_"))	
	   
	  library(VennDiagram)
	  venn.plot <- venn.diagram(x=list(A=A1, B=B1,C=C1),filename = VennName, resolution =500, imagetype = "tiff",col="black",fill = c("red","blue3","green"), 
	                  alpha = 0.5,fontfamily = "serif",fontface = "bold",cat.default.pos = "text",cex = 1.25,cat.col =c("darkred", "darkblue", "darkgreen"),
					  cat.cex = 1.25,main = "A=miRanda   B=miRDB   C=TargetScan",main.pos = c(0.5, 0.2), main.fontface = "plain", main.fontfamily = "serif", 
					  main.col = "black", main.cex = 1.2, main.just = c(0.5, 0.5),cat.fontfamily = "serif", cat.dist = c(0.06, 0.06, 0.03),cat.pos = 0,cat.fontface = "bold",margin = 0.2)
# use following to give svg format

svg('venn.svg',family="GB1" ,width = 8, height = 8)
venn.diagram(x=list(A=A1, B=B1,C=C1),filename = VennName, resolution =500, imagetype = "tiff",col="black",fill = c("red","blue3","green"), 
	                  alpha = 0.5,fontfamily = "serif",fontface = "bold",cat.default.pos = "text",cex = 1.25,cat.col =c("darkred", "darkblue", "darkgreen"),
					  cat.cex = 1.25,main = "A=miRanda   B=miRDB   C=TargetScan",main.pos = c(0.5, 0.2), main.fontface = "plain", main.fontfamily = "serif", 
					  main.col = "black", main.cex = 1.2, main.just = c(0.5, 0.5),cat.fontfamily = "serif", cat.dist = c(0.06, 0.06, 0.03),cat.pos = 0,cat.fontface = "bold",margin = 0.2)
dev.off()



	  x <- dir()
	  y <- x[grep('\\.log.*', x)]
	  file.remove(y)
	  
	}else if(species=="rno"){
	   A=targetgene[targetgene$`miRanda`==1, c("MiRBase ID", "Gene Symbol")]
	   A1=unique(paste(A$`MiRBase ID`, A$`Gene Symbol`, sep="_"))
	   B=targetgene[targetgene$`miRDB`==1, c("MiRBase ID", "Gene Symbol")]
	   B1=unique(paste(B$`MiRBase ID`, B$`Gene Symbol`, sep="_"))
	   C=targetgene[targetgene$`miRTarbase`==1, c("MiRBase ID", "Gene Symbol")]
	   C1=unique(paste(C$`MiRBase ID`, C$`Gene Symbol`, sep="_"))	
	   
	  library(VennDiagram)
	  venn.plot <- venn.diagram(x=list(A=A1, B=B1,C=C1),
  filename = VennName, resolution =500, imagetype = "tiff",
  col="black",fill = c("red","blue3","green"), 
  alpha = 0.5,fontfamily = "serif",fontface = "bold",
  cat.default.pos = "text",cex = 1.25,
  cat.col =c("darkred", "darkblue", "darkgreen"),
  cat.cex = 1.25,
  main = "A=miRanda   B=miRDB   C=miRTarbase",
  main.pos = c(0.5, 0.2), main.fontface = "plain", main.fontfamily = "serif", 
  main.col = "black", main.cex = 1.2,
  main.just = c(0.5, 0.5),cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),cat.pos = 0,
  cat.fontface = "bold",margin = 0.2)

	  x <- dir()
	  y <- x[grep('\\.log.*', x)]
	  file.remove(y)
	}
  #--------------------ɸѡ??��?????ݿ??о??ܱ?Ԥ?⵽?İл???---------------------
   Diff_targetgene <- targetgene[targetgene$`Sum`>=2,]
   Diff_targetgene_unique <- unique(Diff_targetgene$`Gene Symbol`)
   
  #---------------------------------д?ļ?------------------------------------------------------------------------------------------

   write.table(Diff_targetgene, "target_gene_table.txt", row.names = F,col.names = T,sep ='\t',quote=F)
   write.table(Diff_targetgene_unique, "target_gene_list.txt", row.names = F,col.names = F,sep ='\t',quote=F)




