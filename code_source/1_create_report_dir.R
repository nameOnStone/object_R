#Generate a microarray data analysis report
create_report_dir <- function(	path = report_path,
								name = report_name,
								doc_path  = NULL,
							 	comp_name='default'
							) {
  
	dir.create(paste(path,name,sep='\\'), showWarnings = FALSE)
	dir.create(paste(path,paste0(comp_name,'_result_files'),sep='\\'), showWarnings = FALSE)
	main_path <- paste(path, name, sep="\\")
	# tryresult<<-tryCatch({
	#   file.copy(paste(doc, "数据分析报告.docx", sep = "\\"), main_path)
	#   },error=function(e){
	#     cat('error: ',conditonMessage(e))
	#     })
	dir.create(paste(main_path,"1_原始数据分析实验相关内容部分","1_样本信息",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"1_原始数据分析实验相关内容部分","2_芯片质控",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"1_原始数据分析实验相关内容部分","3_实验数据",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"2_探索性分析","1_箱线图",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"2_探索性分析","2_Peason相关系数图",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"2_探索性分析","3_PCA图",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"3_mRNA分析",sep='\\'), showWarnings = FALSE)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,sep='\\'),showWarnings=FALSE)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,"1_差异基因",sep='\\'),showWarnings=FALSE)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Barplot',sep='\\'),showWarnings=FALSE,recursive=T)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Enrishment',sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Result',sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG Result",sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG Barplot & Dotplot",sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","Pathway",sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,"4 注释信息及说明文档",sep='\\'), showWarnings = FALSE)
	exploratory_BoxWhisker_path<<-paste(main_path,"2_探索性分析","1_箱线图",sep='\\')
	exploratory_Peason_path<<-	paste(main_path,"2_探索性分析","2_Peason相关系数图",sep='\\')
	exploratory_PCA_path<<-		paste(main_path,"2_探索性分析","3_PCA图",sep='\\')
	GO_Barplot_path<<-			paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Barplot',sep='\\')
	GO_Enrishment_path<<-		paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Enrishment',sep='\\')
	GO_result_path<<-			paste(main_path,'3_mRNA分析',comp_name,"2_GO分析",'GO_Result',sep='\\')
	KEGG_Result_path<<-			paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG Result",sep='\\')
	KEGG_Barplot_path<<-		paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","KEGG Barplot & Dotplot",sep='\\')
	KEGG_pathway_path<<-		paste(main_path,'3_mRNA分析',comp_name,"3_Pathway分析","Pathway",sep='\\')
	diff_gene_path<<-			paste(main_path,'3_mRNA分析',comp_name,"1_差异基因",sep='\\')
	result_files_path<<-		paste(path,paste0(comp_name,'_result_files'),sep='\\')
}