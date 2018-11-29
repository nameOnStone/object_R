# copy and modified from latest ClariomD_code.R

##GO function annotation
#setwd(path_coding_comp)
#dir.create("GO")
#path_coding_GO <- paste(path_coding_comp,"GO",sep="/")

#use SYMBOL to GO Analysis

GeneID <- as.character("Diff_targetgene_unique")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db", version = "3.8")

library(org.Mm.eg.db)
library(clusterProfiler)


#GeneID <- unique(diff_coding[,"Gene Symbol"])

GeneID <- as.character(Diff_targetgene_unique)

 
#Get functional enrichment results 
go_BP <- enrichGO(
    gene = GeneID,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 1,
    pAdjustMethod = "BH",  #P value adjust method
    pvalueCutoff  = 1,  #p threhold
    qvalueCutoff  = 1
  )
go_CC <- enrichGO(
    gene = GeneID,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "CC",
    minGSSize = 1,
    pAdjustMethod = "BH",  
    pvalueCutoff  = 1,  
    qvalueCutoff  = 1
  )
go_MF <- enrichGO(
    gene = GeneID,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "MF",
    minGSSize = 1,
    pAdjustMethod = "BH",  
    pvalueCutoff  = 1,  
    qvalueCutoff  = 1
  )
BP_result <- as.data.frame(go_BP)
CC_result <- as.data.frame(go_CC)
MF_result <- as.data.frame(go_MF)
 
setwd(path_coding_GO)
if(nrow(BP_result)!=0){
  write.csv(BP_result, "GO_BP.csv", row.names = FALSE)
}else{
  print_error1 <- paste(CompName,"have no enrichGO for unfilter BP ")
  error_data <- rbind(error_data,print_error1)
}
if(nrow(CC_result)!=0){
  write.csv(CC_result, "GO_CC.csv", row.names = FALSE)
}else{
  print_error2 <- paste(CompName,"have no enrichGO for unfilter CC ")
  error_data <- rbind(error_data,print_error2)
}
if(nrow(MF_result)!=0){
  write.csv(MF_result, "GO_MF.csv", row.names = FALSE)
}else{
  print_error3 <- paste(CompName,"have no enrichGO for unfilter MF ")
  error_data <- rbind(error_data,print_error3)
}

##Draw directed acyclic graph for GO result
go_plot_BP <- enrichGO(
    gene = GeneID,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 1,
    pAdjustMethod = "BH", 
    pvalueCutoff  = 0.05,  
    qvalueCutoff  = 1
  )
go_plot_CC <- enrichGO(
    gene = GeneID,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "CC",
    minGSSize = 1,
    pAdjustMethod = "BH",  
    pvalueCutoff  = 0.05,  
    qvalueCutoff  = 1
  )
go_plot_MF <- enrichGO(
    gene = GeneID,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "MF",
    minGSSize = 1,
    pAdjustMethod = "BH",  
    pvalueCutoff  = 0.05,  
    qvalueCutoff  = 1
  )
BP_plot_result <- as.data.frame(go_plot_BP)
CC_plot_result <- as.data.frame(go_plot_CC)
MF_plot_result <- as.data.frame(go_plot_MF)



BP_n <- nrow(BP_plot_result)
if(BP_n!=0){
#  pdf("GO_BP.pdf", family="GB1", height = 10, width = 12)
  svg("GO_BP.svg", family="GB1", height = 10, width = 12)

  if(BP_n>10){
    plotGOgraph(go_plot_BP, firstSigNodes = 10)
  }else{
    plotGOgraph(go_plot_BP,firstSigNodes = BP_n)
  }
  dev.off()
}else{
  print_error4 <- paste(CompName," have no enrichGO for filter BP ")
  error_data <- rbind(error_data,print_error4)
}

CC_n <- nrow(CC_plot_result)
if(CC_n!=0){
  pdf("GO_CC.pdf", family="GB1", height = 10, width = 12)
  if(CC_n >10){
    plotGOgraph(go_plot_CC,firstSigNodes = 10)
  }else{
    plotGOgraph(go_plot_CC,firstSigNodes = CC_n)
  }
  dev.off()
}else{
  print_error5 <- paste(CompName,"have no enrichGO for filter CC ")
  error_data <- rbind(error_data,print_error5)
}

MF_n <- nrow(MF_plot_result)
if(MF_n!=0){
  pdf("GO_MF.pdf", family="GB1", height = 10, width = 12)
  if(MF_n >10){
    plotGOgraph(go_plot_MF,firstSigNodes =10)
  }else{
    plotGOgraph(go_plot_MF,firstSigNodes = MF_n)
  }
  dev.off()
}else{
  print_error6 <- paste(CompName,"have no enrichGO for filter MF")
  error_data <- rbind(error_data,print_error6)
}

##Draw barplot for GO result
BP_rank <- BP_plot_result[order(BP_plot_result[,"p.adjust"]),]
CC_rank <- CC_plot_result[order(CC_plot_result[,"p.adjust"]),]
MF_rank <- MF_plot_result[order(MF_plot_result[,"p.adjust"]),]

if(nrow(BP_rank)>20){
  BP_plot <- BP_rank[1:20,] 
}else{
  BP_plot <- BP_rank
}

if(nrow(CC_rank)>20){
  CC_plot <- CC_rank[1:20,] 
}else{
  CC_plot <- CC_rank
}

if(nrow(MF_rank)>20){
  MF_plot <- MF_rank[1:20,] 
}else{
  MF_plot <- MF_rank 
}
  
GO_plot <- rbind(BP_plot, CC_plot, MF_plot)
GO_bar <- as.numeric(GO_plot[,"Count"])
GO_name <- as.character(GO_plot[,"Description"])

GeneRatio_GO <- GO_plot$`GeneRatio`
gene_rate <- as.numeric(unlist(strsplit(GeneRatio_GO, split = "/")))
gene_numbers <- matrix(data = gene_rate, ncol = 2, byrow = TRUE)
GO_plot$`GeneRatio` <- gene_numbers[,1]/gene_numbers[,2]


setwd(path_coding_GO)
if(length(GO_bar)!=0){  

  #barplot not working
  pdf(paste(CompName,"GO_barplot.pdf",sep="_"), width = 5.8+0.1365*nrow(GO_plot), 
     height = 5.28+0.045*max(nchar(GO_name)))
    #width = 1.8+0.1365*nrow(GO_bar)
    if(0.0378*nchar(GO_name[1]) >= 1.2){
      par(mai=c(0.045*max(nchar(GO_name))+1.28,0.0378*nchar(GO_name[1])+1,1,0.6))
    }else{
      par(mai = c(0.045*max(nchar(GO_name))+1.28,1.2,1,0.6))
    }
    xbar=barplot(GO_bar, space =0.5, las = 2, col =c(rep("firebrick1",nrow(BP_plot)),
             rep("goldenrod1",nrow(CC_plot)),rep("deepskyblue2",nrow(MF_plot))),
             border = NA, main = 'Gene Function Classification (GO)', cex.main=1, 
             ylab = 'Numbers of genes', xpd = T, axisnames = T,
	       cex.axis = 0.7, ylim=c(0, max(GO_bar)+3))
    legend("topright", legend = c("Biological Process","Cellular Component","Molecular Function"), 
           fill=c("firebrick1","goldenrod1","deepskyblue2"),cex = 0.7, border = NA)
    text(x=xbar, y=-(0.02*max(GO_bar)+0.12), labels = GO_name, srt = 50, adj = 1, xpd = T, cex = 0.7)
  dev.off()

### by de novo develop14.R

mf.table <- BP_plot
bp.table <- CC_plot
cc.table <- MF_plot

table <- rbind(mf.table, bp.table, cc.table)

color = c( rep( "deepskyblue2", nrow(mf.table) ),
           rep( "springgreen4", nrow(bp.table) ),
           rep( "goldenrod2", nrow(cc.table) )
)


term <- table$Description

need_short_idx <- nchar(term) > 40

term[need_short_idx] <- paste(
  substr(term[need_short_idx], 1, 37), "...", sep = "")

svg("GO_Boxplot.svg", family="GB1", height = 10, width = 12)

par(mar = c(12,8,4,2) + 0.2)
barplot(table$Count, col = color, main = "GO Term summary",
        ylab = "Number of Gene", )

num <- nrow(mf.table) + nrow(bp.table) + nrow(cc.table)

text(seq(1, 1.2*num, 1.2), par("usr")[3]-0.1, labels = term,
     adj = 1.05, cex = 0.8, srt = 45, xpd = TRUE)

# borrow from RNA-seq
legend("top", inset=.05,
  c("Biological Process","Cellular Component","Molecular Function"),
  fill = unique(color), cex = 0.7, horiz=F,
  border = NA)


# extend before you save                
savePlot("sigTerm", type = "pdf")
savePlot("sigTerm", type = "png")



dev.off()


###


  #dotplot	
  GO_dot <- ggplot(GO_plot)+
            geom_point(aes(x=GO_plot$GeneRatio,y=GO_plot$Description,
                       colour=GO_plot$p.adjust,size=GO_plot$Count))+
	  	scale_colour_gradientn(colours = c("red","green"),
                       limit=c(0,max(GO_plot$p.adjust)))+
  		labs(colour="p.adjust",size="Gene_number",x="GeneRatio", y="")+
		theme_bw()+
		ggtitle("Statistics of GO Enrichment")
    ggsave("GO_dotplot.svg", GO_dot, width = 10, height = 7)



}else{}

#运行KEGG1小函数
viewKEGG1 <- function(obj, pathwayID, foldChange,
                      color.low="green",
                      color.high="red",
                      species = "hsa",
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
    pathwayID <- summary(obj)[pathwayID, 1]
  }
  if (length(pathwayID) == 1 & pathwayID == "all") {
    pathwayID <- summary(obj)[, 1]
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
             out.suffix=out.suffix,
             new.signature=FALSE)
  })
  return (res)
}

#KEGG通路注释
#setwd(path_coding_comp)
#dir.create("KEGG")
#path_coding_KEGG <- paste(path_coding_comp,"KEGG",sep="/")
#setwd(path_coding_KEGG) 

unigene_entrez <- bitr(gene = GeneID,
	               fromType = "SYMBOL",
	               toType = "ENTREZID",
	               OrgDb = org.Mm.eg.db)
  
KEGG <- enrichKEGG(gene          = unigene_entrez$ENTREZID,
                   keyType           = "kegg",
                   organism          = "mmu",
                   pAdjustMethod     = "BH",
                   minGSSize         = 1,
                   pvalueCutoff      = 1,
                   qvalueCutoff      = 1,
                   use_internal_data = FALSE)

KEGG_result <- as.data.frame(KEGG)
if(nrow(KEGG_result)!=0){
  KEGG_gene <- as.matrix(KEGG_result[,'geneID'])
  KEGG_symbol <- NULL
  for(i in 1:length(KEGG_gene)){
    entrezID <- as.numeric(unlist(strsplit(KEGG_gene[i,],"/")))
    ChangeID <- bitr(gene = entrezID,
  	           fromType = "ENTREZID",
	           toType = "SYMBOL",
	           OrgDb = "org.Mm.eg.db")
    Symbol <- ChangeID[,"SYMBOL"]
    KEGG_symbol <- rbind(KEGG_symbol,as.matrix(paste(Symbol,collapse="/")))
  }

  KEGG_result <- cbind(KEGG_result[,1:7],
  "Gene Symbol"=KEGG_symbol,"Count"=KEGG_result[,9])

  write.table(KEGG_result,"KEGG_result.csv", sep=",", row.names=FALSE)
  
  KEGG2 <- enrichKEGG(gene          = unigene_entrez$ENTREZID,
                   keyType           = "kegg",
                   organism          = "mmu",
                   pAdjustMethod     = "BH",
                   minGSSize         = 1,
                   pvalueCutoff      = 0.05,
                   qvalueCutoff      = 1,
                   use_internal_data = FALSE)
  KEGG_filter <- as.data.frame(KEGG2)

  ##KEGG结果柱状图
  KEGG_bar_name <- "KEGG_barplot.svg"

  if(nrow(KEGG_filter)!=0){
    KEGG_rank <- KEGG_filter[order(KEGG_filter[,"p.adjust"]), ]
    if(nrow(KEGG_rank)>20){
      KEGG_plot <- KEGG_rank[1:20, c("Description", "Count")] 
    }else{
    KEGG_plot <- KEGG_rank[, c("Description", "Count")] 
    }
		
    KEGG_count <- as.numeric(KEGG_plot[, "Count"])
    KEGG_name <- as.character(KEGG_plot[, "Description"])
	
    if(nrow(KEGG_plot) >= 7){
	  svg(KEGG_bar_name, width = 0.444*nrow(KEGG_plot)+1.75, 
           height = 0.0437*max(nchar(KEGG_name))+3.726)
	    if((0.0367*nchar(KEGG_name[1])+0.11)>=1.25){
	      par(mai=c(0.0437*max(nchar(KEGG_name))+0.226,0.0367*nchar(KEGG_name[1]),1,0.5))
	    }else{
	      par(mai=c(0.0437*max(nchar(KEGG_name))+0.226,1.15,1,0.6))
	    }
	    xbar <- barplot(KEGG_count, width = rep(0.7,nrow(KEGG_plot)),space =0.7, las = 2, 
	                 col =rep("firebrick1", times = nrow(KEGG_plot)),border = NA, 
	                  main = 'Gene Function Classification (KEGG)', ylab = 'Numbers of genes',
					  xpd = T, axisnames = T, cex.main=1, cex.axis=0.7, ylim=c(0,max(KEGG_count)))
	    text(x=xbar, y=-0.038*max(KEGG_count), labels = KEGG_name, srt = 50, adj = 1, cex = 0.7, xpd = T)
     dev.off()
   }else{
     pdf(KEGG_bar_name, width = 5, height = 0.0437*max(nchar(KEGG_name))+3.726)
	   par(mai=c(0.0437*max(nchar(KEGG_name))+0.226,2.5-0.222*nrow(KEGG_plot),1,2.5-0.222*nrow(KEGG_plot)))
	   xbar <- barplot(KEGG_count, width = rep(0.7,nrow(KEGG_plot)),space =0.7, las = 2, 
	                 col =rep("firebrick1", times = nrow(KEGG_plot)),border = NA,
	                 main = 'Gene Function Classification (KEGG)', ylab = 'Numbers of genes',xpd = T, 
	                 axisnames = T, cex.main=1, cex.axis=0.7, ylim=c(0,max(KEGG_count)))
	   text(x=xbar, y=-0.038*max(KEGG_count), labels = KEGG_name, srt = 50, adj = 1, cex = 0.7, xpd = T)
     dev.off()
     }
  }else{
    print_error7 <- paste(CompName,"have no enrichKEGG for filter KEGG ")
    error_data <- rbind(error_data,print_error7)
  }

  #dotplot
  GeneRatio <- KEGG_filter$`GeneRatio`
  b <- as.numeric(unlist(strsplit(GeneRatio, split = "/")))
  gene_numbers <- matrix(data = b, ncol = 2, byrow = TRUE)
  KEGG_filter$`GeneRatio` <- gene_numbers[,1]/gene_numbers[,2]
	
  KEGG_rank2 <- KEGG_filter[order(KEGG_filter$p.adjust), ]
  if(nrow(KEGG_rank2)!=0){
    if(nrow(KEGG_rank2)>20){
      KEGG_count <- KEGG_rank2[1:20, ] 
    }else{
      KEGG_count <- KEGG_rank2
    }
    c <- KEGG_count[order(KEGG_count$p.adjust, decreasing = TRUE),]
	
    KEGG_dot <- ggplot(c)+geom_point(aes(x=c$GeneRatio,y=c$Description,colour=c$p.adjust,size=c$Count))+
	  	scale_colour_gradientn(colours = c("red","green"),limit=c(0,max(c$p.adjust)))+
  		labs(colour="p.adjust",size="Gene_number",x="GeneRatio", y="")+
		theme_bw()+
		ggtitle("Statistics of Pathway Enrichment")
    ggsave("dotplot.svg", KEGG_dot, width = 10, height = 7)
  }else{}

#  dir.create("pathview")
#  setwd("pathview")





# 以下内容不清楚功能
  #整理KEGG表格
  #差异基因去NA和无基因名
  noN_gene <- diff_coding[which(diff_coding[,"Gene Symbol"]!= "NA" 
                                & diff_coding[,"Gene Symbol"]!= ""),c("Fold Change","Gene Symbol")]
  noN_gene <- as.data.frame(noN_gene)
  rownames(noN_gene)<-NULL
  repeat_times <- as.data.frame(table(as.character(noN_gene$'Gene Symbol')))
  #得到重复次数大于1次的基因名和其FC值
  repeat_n <- repeat_times[which(repeat_times$Freq>1),]
  #取重复区唯一基因的名称和FC值
  repeat_unique <- NULL 
  if(nrow(repeat_n)!=0){
    for (i in 1:nrow(repeat_n)) {
       repeat_gene <- noN_gene[noN_gene[,2]==repeat_n[i,1],]
       repeat_FC <- as.numeric(repeat_gene[,"Fold Change"])
       correct_FC <- NULL
       if (all(repeat_FC>0,na.rm=FALSE)==T|all(repeat_FC<0,na.rm=FALSE)==T) {
         if (repeat_FC[1]>0) {
            correct_FC <- max(repeat_FC)
         }else{
            correct_FC <- min(repeat_FC)
         }
       }else{
         correct_FC<-NULL
       } 
      correct_gene <- unique(repeat_gene[repeat_gene[,"Fold Change"]==correct_FC,])
      repeat_unique <- rbind(repeat_unique,correct_gene)
    }
   #与去差异基因NA和无基因名的数据框合并
    no_repeat_gene <- noN_gene[!duplicated(noN_gene$"Gene Symbol"),]
    gene_FC <- rbind(no_repeat_gene[-match(repeat_unique[,"Gene Symbol"],no_repeat_gene[,"Gene Symbol"]),],
                        repeat_unique)
  }else{
    gene_FC <- noN_gene
  }
  #取其FC值，若为负，则求其负倒数
  KEGG_all_FC <- as.numeric(gene_FC[,"Fold Change"])
  neg <- KEGG_all_FC < 0
  KEGG_all_FC[neg] <- -1/KEGG_all_FC[neg]
  loc_mappinggene <- match(unigene_entrez$`SYMBOL`,gene_FC[,"Gene Symbol"])
  KEGG_FC <- KEGG_all_FC[loc_mappinggene]
  names(KEGG_FC) <- unigene_entrez$ENTREZID
  KEGG_log2FC <- log2(KEGG_FC)

  if(nrow(as.data.frame(KEGG))>=20){
    n_plot <- c(1:20)
  }else{
    n_plot <- "all"
  }
  kegg_pathway <- tryCatch(pathway_coding <- viewKEGG1(KEGG,
                                                    pathwayID = n_plot,
                                                    foldChange = KEGG_log2FC,
                                                    species = "hsa",
                                                    color.low="green",
                                                    color.high="red",
                                                    kegg.native=TRUE,
                                                    out.suffix="pathway"),
	                  error=function(e){print(paste(CompName,"pathwayplot-error",
                                             conditionMessage(e),"\n",sep="_"))})
   x.inv6 <- try(pathway_coding <- viewKEGG1(KEGG,
                                         pathwayID = n_plot,
                                         foldChange = KEGG_log2FC,
                                         species = "hsa",
                                         color.low="green",
                                         color.high="red",
                                         kegg.native=TRUE,
                                         out.suffix="pathway"),											   
		   silent = TRUE)
   if ('try-error' %in% class(x.inv6)){
   }else {
     pathway_lnc <- x.inv6
   }
   if(length(dev.list())>0){
	replicate(length(dev.list()),dev.off())
   }
   plotnames <- list.files()
   loc <- grep('.pathway.png',plotnames)
   removeplot <- plotnames[-loc]
   file.remove(removeplot)

   pathway_error <- kegg_pathway[grep("-error_",kegg_pathway)]
   pathway_error <- as.matrix(pathway_error)
   error_data <- rbind(error_data,pathway_error)
   }else{
    print_error7 <- paste(CompName,"have no enrichKEGG for unfilter KEGG ")
    error_data <- rbind(error_data,print_error7)
}
