#设置必要路径
TAC_path <- "F:\\王铎霖做过的数据备份\\KPS082017044_首都儿科研究所AFFY HTA2.0表达谱芯片技术服务合同\\Project_File\\Require_of_Gene.2.0.R\\Coding\\test"
#TAC数据存放路径，请据实修改
report_path <- "F:\\王铎霖做过的数据备份\\KPS082017044_首都儿科研究所AFFY HTA2.0表达谱芯片技术服务合同\\Project_File\\Require_of_Gene.2.0.R\\Coding"
#数据报告生成路径，请据实修改
report_name <-"wangduolintest"
organism <- "human"    #物种，人("human")、小鼠("mouse")、大鼠("rat")
R_path <- "G:\\Desktop\\About_Framework\\Gene_2.0_Code_and_instruction\\Father_R"
file_Name_ID2ID <- "F:\\王铎霖做过的数据备份\\KPS082017044_首都儿科研究所AFFY HTA2.0表达谱芯片技术服务合同\\Project_File\\Require_of_Gene.2.0.R\\Coding\\ChipIDSampleID.txt"
instruction_files_path <- ""
FC <- 2
pvalue <- 0.05
#说名文档所在路径，请据实修改

library('openxlsx')
library('gplots')
library('tidyverse')

FunList <- dir(R_path)
sapply(FunList, function(x) source(paste(R_path,x,sep='\\')))

for(i in dir(TAC_path)){
	comp_name <-gsub('\\..*','',i)
	create_report_dir(	path 		= report_path,
						name		= report_name,
						doc_path	= NULL,
						comp_name	= comp_name)

	file_Name_Gene <- paste(TAC_path,i,sep='\\')
	mRNA_analysis(	filename1		= file_Name_Gene,
					filename2		= file_Name_ID2ID,
					comp_name 		= comp_name,
					organism 		= organism,
					FC 				= 2,
					pvalue 			= 0.05)
	#以下的变量已经成为全局变量，为避免每次循环时可能会使用这些全局变量，特将下列全局变量删除！
	rm(allData)
	rm(diffresult)
	rm(GOresult)
	rm(KEGGresult)
}
