---
title: "lncRNA"
author: "wangduolin"
date: "2018年11月21日"
output: html_document
---
  # lncRNA分析总结

---
今天，我将从以下几个方面来说：


1. 基本信息
2. 实验流程
3. 作用机制
4. 数据库介绍
5. 如何分析（重点）

# 基本信息：

大于200nt的非转录RNA为 long-non RNA；分为Sence(相向链重叠) Antisence(相斥链重叠) Bidirectional(相斥不重叠) Intronic(lncRNA来自一个transcript内部) Intergenic(相向链不重叠)
有部分lncRNA具有polyA尾巴，所以不能用磁珠方法调取

特点：有强烈的组织、时空特异性；自身变化多，因此可与转录因子、酶类结合；

# 它的作用：
参与转录调控（分为转录前、中、后调控）
转录前调控机制：转录前调控就是对转录事件发生之前的准备阶段进行调控，主要指染色质开放或关闭状态的调控（比如，想让不转录的基因开启转录，需要使染色质由关闭状态转化为开启状态，反之亦然），染色质的状态的改变过程：状态是由表观修饰因子决定（富含H3K4me3、H3K36me3及组蛋白乙酰化等激活型组蛋白修饰的为开放状态；富含H3K9me3、H3K27me3、H4K20me3及DNA甲基化等抑制型组蛋白修饰的为关闭状态），但是问题来了，表观修饰因子本质是蛋白。结构固定，怎么去和千变万化的染色质结合呢？答案是由lncRNA的中间商做牵线搭桥（比如是由lncRNA中的一部分结构和DNA结合，这样会导致折叠成高级结构，然后再和钢铁直男修饰因子结合）

目前，先进的测序技术表明：哺乳动物中lncRNA中有几万条，虽然目前认为大多数lncRNA是有功能的，但是，目前只要极少数被注释是有功能的（不到200个）；不过，目前也有很多数据表明，先前被证明是lncRNA，后来发现它有翻译蛋白的功能；
RNA转录在真核生物中是一个受到严密调控的过程。非编码RNA可以靶向该进程的多个方面，包括靶向转录激活因子或转录抑制因子、如RNA聚合酶（RNAP）Ⅱ等转录反应中的各组分、甚至是DNA双螺旋结构，以达到调控基因转录及表达的目的（Goodrich 2006）。非编码RNA将这些机制结合在一起可以组成为一个包括转录因子在内的调控网络，可以精细地调控复杂真核生物的基因表达；
lncRNA可调控transcription factor activity，作为 co-regulators行使功能

不仅如此：

In post-transcriptional regulation
与小调控RNAs，例如微小RNAs和小核仁RNAs，类似，ncRNAs 的功能包括与目标mRNA进行互补碱基配对。互补ncRNA和mRNA形成的RNA双链可能为需要结合反式作用因子的mRNA募集关键因子，可能影响转录后水平基因表达的每一步，包括前体mRNA加工，剪接，运输，翻译以及降解。

# 实验流程：

![](E://lncRNA_liucheng.png)

#作用机制：转录前、转录中、转录后

转录前：主要是对染色质开放或关闭状态的调控

转录中：(1) 通过结合增强子控制TFs与增强子的互作从而调控增强子活性；(2) 通过转录事件抑制RNA polII复合体与启动子结合从而实现转录干扰；(3) 因其与DNA很像，可作为诱饵捕获TFs，从而控制转录因子的活性（其实，只要涉及到蛋白质和DNA互作，lncRNA就有机可乘）。

转录后：一般是说对transcripts的后期加工和转运过程以及稳定性方面进行调控，前面我们说了，lncRNA的结构变化比较容易，所以：比如：lncRNA结合到RNA的剪切位点，则此处不能被切割了。再比如：lncRNA可以和mRNA结构互补形成更加容易降解或容易稳定的结构。再比如：lncRNA可以和miRNA结构互补，从而间接影响miRNA对mRNA的抑制或降解作用。

#数据库介绍：

主要有三大数据库：NCBI的Entrez。Ensemble（欧洲，自动注释，准确度一般，根据不同物种设置的前缀 + 数据类型（基因/蛋白质）+ 一系列的数字（例如版本号就以小数点体现））。Uniprot（之所以分三个是因为由日本 欧洲 美国三国数据库合并起来导致的，因此也可能导致东西一样，但是命名却不一样）；再由这三大数据库做主框架数据库，分支出其他的主流数据库，如Refseq（冗余度少，采用HGNC命名，可信度高,两个大写字母+一个下划线+大于6个数字）。GenBank（冗余较高，准确度一般）。dbSNP。GO（根据基因产物功能进行注释的，而不是基因功能）。KEGG数据库。以上数据库知道下就行；

1.NONCODE:
主要关注到lncRNA，共包含17个物种（human, mouse, cow, rat, chicken, fruitfly, zebrafish, celegans, yeast, Arabidopsis, chimpanzee, gorilla, orangutan, rhesus macaque, opossum platypus and pig），它也通过pubmed,Ensembl , RefSeq, lncRNAdb and GENCODE中用关键词在Supplementary Material中去收集ncRNA信息；
数据库中的ID是根据自己网站规则起的（但是，别担心，它自己有ID转换工具）；

lncRNA的命名规则：
人的序列是由唯一官方授权机构雨果命名委员会（HGNC）命名，对其中小一万个（约8500）个非编码进行了命名，其中，microRNA功能保守，比较好命名，但是，lncRNA则由于多变的功能及没有保守的序列而导致不好命名，一般情况，都是由lncRNA的产物功能来命名，但如果不知道产物功能的lncRNA则采用它位置信息（比如挨着哪些基因）来命名


#昂飞芯片：HTA芯片
它是WT芯片（whole Transcriptome），覆盖全，而且真对不同的剪切后的转录本设计了不同的探针，而且可以检测到lncRNA；
康普森公司采用的是一种叫“寡核苷酸探针”，长度在20-25bp之间，不过，由于较短，所以存在寡核苷酸的非特异性结合风险，对此，昂飞采用的对策是：对于同一根探针制作两根探针，一根完全匹配（MM），一根允许一个碱基错配（PM，这根探针是在第13个碱基被替换了碱基），然后将他们之间的表达差值作为信号强度；同时，一个序列设计多根独立探针，保证信号强度准确；

其中，CEL files 被当成“raw”data；但是，CEL文件并不是all the information；一般还需要CDF文件（Chip Description File)，它应该就是注释文件；它有Bioconductor CDF packages；





#如何分析：

	芯片的大致步骤：
	大致分为以下几步以及各个步骤的作用：
	1、ReadData #读入数据
    2、QualityControl #数据的质量控制
	3、backgroud adjustment
	#这里面包括信号背景处理和噪声处理（注意：噪声和背景是两种概念，要区分好），如果按照简单的处理方式则是：PM-MM得出的值则为强度,			但是，实际情况是：有多达30%的MM型探针的信号强度比			相应的PM型的信号强度还强。因此，研究者发明了两种方法（MAS和RMA），两种方法的细节没有探究
	4、DataNormalization #将数据归一化以利于芯片之间进行比较；
	5、DifferentiallyExpressionGenes #清洗数据，寻找有意义的Genes
	6、GO analysis
	7、KEGG analysis
		
	RNA_seq的大致步骤：
		1.map:首先将测序读段映射至参考基因组以lncRNA研究为目标的RNA-Seq数据	处理流程，主要采用能够识别剪接读段的映射方法。  
		可使用Tophat，此软件通过调用Bowtie或者Bowtie2将read比对到参考基因组上去，然后分析比对结果，从而找出外显子之间的结合位点。
		2.align：常见的有cufflinks软件，此软件可主要根据Tophat的比对结果，依托或不依托参考基因组GTF注释文件，计算出各个transcript的FPKM值，并给出transcripts.gtf的文件。然后我们需要对RNA进行分类使用Cuffcompare软件，
		3.如何鉴别是lncRNA呢？首先提取转录本长度超过200nt的，其次，过滤掉是编码蛋白质的转录本

参考：

[网址：https://www.jianshu.com/p/c3df672d9794](网址：https://www.jianshu.com/p/c3df672d9794 "数据库说明")

	
	

