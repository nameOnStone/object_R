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
	#   file.copy(paste(doc, "���ݷ�������.docx", sep = "\\"), main_path)
	#   },error=function(e){
	#     cat('error: ',conditonMessage(e))
	#     })
	dir.create(paste(main_path,"1_ԭʼ���ݷ���ʵ��������ݲ���","1_������Ϣ",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"1_ԭʼ���ݷ���ʵ��������ݲ���","2_оƬ�ʿ�",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"1_ԭʼ���ݷ���ʵ��������ݲ���","3_ʵ������",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"2_̽���Է���","1_����ͼ",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"2_̽���Է���","2_Peason���ϵ��ͼ",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"2_̽���Է���","3_PCAͼ",sep='\\'), recursive = T,showWarnings=FALSE)
	dir.create(paste(main_path,"3_mRNA����",sep='\\'), showWarnings = FALSE)
	dir.create(paste(main_path,'3_mRNA����',comp_name,sep='\\'),showWarnings=FALSE)
	dir.create(paste(main_path,'3_mRNA����',comp_name,"1_�������",sep='\\'),showWarnings=FALSE)
	dir.create(paste(main_path,'3_mRNA����',comp_name,"2_GO����",'GO_Barplot',sep='\\'),showWarnings=FALSE,recursive=T)
	dir.create(paste(main_path,'3_mRNA����',comp_name,"2_GO����",'GO_Enrishment',sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA����',comp_name,"2_GO����",'GO_Result',sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA����',comp_name,"3_Pathway����","KEGG Result",sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA����',comp_name,"3_Pathway����","KEGG Barplot & Dotplot",sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,'3_mRNA����',comp_name,"3_Pathway����","Pathway",sep='\\'),showWarnings=FALSE,recursive = T)
	dir.create(paste(main_path,"4 ע����Ϣ��˵���ĵ�",sep='\\'), showWarnings = FALSE)
	exploratory_BoxWhisker_path<<-paste(main_path,"2_̽���Է���","1_����ͼ",sep='\\')
	exploratory_Peason_path<<-	paste(main_path,"2_̽���Է���","2_Peason���ϵ��ͼ",sep='\\')
	exploratory_PCA_path<<-		paste(main_path,"2_̽���Է���","3_PCAͼ",sep='\\')
	GO_Barplot_path<<-			paste(main_path,'3_mRNA����',comp_name,"2_GO����",'GO_Barplot',sep='\\')
	GO_Enrishment_path<<-		paste(main_path,'3_mRNA����',comp_name,"2_GO����",'GO_Enrishment',sep='\\')
	GO_result_path<<-			paste(main_path,'3_mRNA����',comp_name,"2_GO����",'GO_Result',sep='\\')
	KEGG_Result_path<<-			paste(main_path,'3_mRNA����',comp_name,"3_Pathway����","KEGG Result",sep='\\')
	KEGG_Barplot_path<<-		paste(main_path,'3_mRNA����',comp_name,"3_Pathway����","KEGG Barplot & Dotplot",sep='\\')
	KEGG_pathway_path<<-		paste(main_path,'3_mRNA����',comp_name,"3_Pathway����","Pathway",sep='\\')
	diff_gene_path<<-			paste(main_path,'3_mRNA����',comp_name,"1_�������",sep='\\')
	result_files_path<<-		paste(path,paste0(comp_name,'_result_files'),sep='\\')
}