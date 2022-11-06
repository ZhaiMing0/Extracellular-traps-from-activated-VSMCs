
data_dir<-"~/Desktop/ZM/PR"
list.files(data_dir)
target<- Read10X(data.dir=data_dir) #读取数据
dim(target)#查看维度，基因数以及细胞数，也可以在environment中查看
target[1:10,1:6]#查看表达矩阵（1-10行，1-6列）
#对于GEO处理后提供的一些原始数据的处理
target<-readRDS("GSM3823940_control.s2.dgecounts")
#RDS文件包括所有的测序数据，需要将其中的的外显子测序数据提取出来#
data = scRNActr2[["readcount"]][["exon"]][["all"]] #直接View数据包，然后在目录上单击即可提取#
#创建seurat对象和数据过滤
#数据集中测到的少于200个基因的细胞（min.features = 200）和少于3个细胞覆盖的基因（min.cells = 3）被过滤掉
target <- CreateSeuratObject(counts = target, min.cells = 3, min.features = 200)  
rm(MI) #删除矩阵
#两个不同数据的合并#注意合并需要进行不同的分组，之后如果分组进行细胞群落划分有依据#
seurat_object <- merge(Case,Control)
#计算每个细胞的线粒体基因转录本数的百分比（%），使用[[ ]] 操作符存放到 metadata 中； 
target[["percent.mt"]] <- PercentageFeatureSet(target, pattern = "^MT-")
#nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和，percent.mt代表测到的线粒体基因的比例。
VlnPlot(target, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.pdf", width = 28, height = 25, units = "cm")

#过滤细胞：保留 gene 数大于 200 小于 2500 的细胞；目的是去掉空 GEMs 和 1 个 GEMs 包 含 2 个以上细胞的数据；
#而保留线粒体基因的转录本数低于 5%的细胞，为了过滤掉死细胞 等低质量的细胞数据。 
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat_object<-target
## 对过滤后的 QC metrics 进行可视化（绘制散点图）；
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#表达量数据标准化：LogNormalize 的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 ) 
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
#鉴定细胞间表达量高变的基因（feature selection），用于下游分析，PCA
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# 提取表达量变变化最高的 10 个基因； 
top10 <- head(VariableFeatures(seurat_object), 10)
top10
# 绘制带有和不带有标签的变量特征的散点图
plot3 <- VariableFeaturePlot(seurat_object)+NoLegend()
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge=0, ynudge=0)
plot3+plot4

#PCA分析数据准备，使用ScaleData()进行数据归一化；默认只是标准化高变基因（2000 个），速度更快，不影响 PCA 和分群，但影响热图的绘制。 
scale.genes <-  VariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object, features = scale.genes)
#而对所有基因进行标准化的方法如下： 
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes, vars.to.regress = "percent.mt") ##耗时2min
#线性降维（PCA）,默认用高变基因集，但也可通过 features 参数自己指定； 
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object)) 
# 检查 PCA 分群结果， 这里只展示前 12 个 PC,每个 PC 只显示 3 个基因； 
print(seurat_object[["pca"]], dims = 1:12, nfeatures = 3) 
#绘制 pca 散点图； 去除图例
DimPlot(seurat_object, reduction = "pca")+ NoLegend() 
#画前 2 个主成分的热图； 
DimHeatmap(seurat_object, dims = 1:2, cells = 500, balanced = TRUE) 
ggsave("plot2.png", width = 28, height = 25, units = "cm")

#确定数据集的分群个数 
##方法 1：Jackstraw 置换检验算法；重复取样（原数据的 1%），重跑PCA,鉴定p-value较小的PC；计算‘null distribution’(即零假设成立时)时的基因 scores; 
seurat_object <- JackStraw(seurat_object, num.replicate = 100)  ##耗时3min
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20) 
JackStrawPlot(seurat_object, dims = 1:15)
#方法 2：肘部图（碎石图），基于每个主成分对方差解释率的排名； 
ElbowPlot(seurat_object)

#基于PCA空间中的欧氏距离计算 nearest neighbor graph，优化任意两个细胞间的距离权重（输入上一步得到的 PC 维数）； 
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
#接着优化模型，resolution 参数决定下游聚类分析得到的分群数，对于3K左右的细胞，设为0.4-1.2能得到较好的结果(官方说明)；如果数据量增大，该参数也应该适当增大； 
seurat_object <- FindClusters(seurat_object, resolution = 0.5) 
#使用 Idents（）函数可查看不同细胞的分群； 查看前8个细胞的分群ID
head(Idents(seurat_object), 5) 

#tsne非线性降维
seurat_object <- RunTSNE(seurat_object, dims = 1:10) 
#用TSNEPlot函数绘制tsne图
tsneplot<-TSNEPlot(seurat_object,label = TRUE, pt.size = 1.5)+ NoLegend()
tsneplot
#用DimPlot函数绘制tsne图
tsneplot1 <- DimPlot(seurat_object, reduction = "tsne", pt.size = 1.5)
tsneplot1
#注意：在这里，如果需要按照merge前的两组进行分组，需要使用group.by函数# #采用View，去观察meta_data下面的数据结构，并找寻之前命名project的term依据，并以此分群#
tsneplot1 <- DimPlot(seurat_object, group.by = "orig.ident", reduction = "tsne", pt.size = 0.05)
tsneplot1
#绘制 Marker 基因的 tsne 图； 
P2<-FeaturePlot(seurat_object, features = c("MYH6"),pt.size = 0.2)
P2
#UMAP非线性降维
seurat_object <- RunUMAP(seurat_object, dims = 1:10, label = T)
head(seurat_object@reductions$umap@cell.embeddings) # 提取UMAP坐标值。
#用DimPlot函数绘制UMAP图
umapplot <- DimPlot(seurat_object, reduction = "umap",pt.size = 0.8)
umapplot
#保存目前分析的数据#
save(seurat_object,file = "ctr+case merge.Rda")
#或者保存为RDS格式#
saveRDS(scRNA1, file="scRNA1.rds")

#接下来是对细胞亚群的结果进行鉴定，有三种鉴定的办法#
####细胞类型的注释一般有三种方法.1、利用marker基因查找网站进行注释  2、使用singler进行注释 3、根据已有的生物学知识或者文献，按照dotplot来注释。
##现在使用方法一寻找marker基因使用网站注释 找marker基因有以下方法三选一，建议第一种
#默认wilcox方法
diff.wilcox = FindAllMarkers(seurat_object)
##3 %>% 这是通道函数 起传递左右  可以自己百度深入理解一下，我这里只告诉你他是起传递作用的
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
###将marker基因保存一下
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
###把top10marker基因保存一下
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)


new.cluster.ids <- c("0", "1", "Padi4 high macrophage", "Padi4 high macrophage", "Stemlike-SMC", "5", "6", "Macrophage", "SMC","9") #自定义名称
names(new.cluster.ids) 
scRNA1<-seurat_object

levels(SS)
#将seurat_object的水平属性赋值给new.cluster.ids的names属性； 
names(new.cluster.ids) <- levels(scRNA1)
names(new.cluster.ids) 
scRNA1 <- RenameIdents(scRNA1, new.cluster.ids) 


library(monocle)
####抽样####
sc.sub=sample(colnames(T1),1500)

scRNA.Osteoclastic=scRNA1[,sc.sub]

data <- as(as.matrix(Good@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Good@meta.data)    #new函数是创建一个S4对象#
?new
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) #根据不同数据类型进行微调#
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())  


monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)  

######选择Seurat差异基因来做
deg.cluster  <- read.csv("diff_genes_wilcox.csv", header = FALSE)
colnames(deg.cluster)<-deg.cluster[1,]
deg.cluster<-deg.cluster[-1,]

diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
HSMM <- setOrderingFilter(HSMM, diff.genes)
plot_ordering_genes(HSMM)

setwd("~/Desktop")




###选择高变基因
HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

#####降低维度
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)

###降维
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
###空间维度 max_components = 2
####  按照轨迹排序细胞，这一步就产生state
HSMM <- orderCells(HSMM)
###可视化排序结果
plot_cell_trajectory(HSMM, color_by = "State")

plot_cell_trajectory(HSMM,color_by = "cellType",cell_size=0.6)

?brewer.pal

cell_type_cols<-c("#FC8D62","#8DA0CB","#E78AC3","#E5C494","#00CD00")  ####good####
cell_type_cols<-c("#66C2A5","#8DA0CB","#A6D854","#FFD92F","#90EE90")  ####bad###

plot_cell_trajectory(HSMM, color_by = "cellType",cell_size = 0.65) + scale_color_manual(values=cell_type_cols) 

?plot_cell_trajectory

plot_cell_trajectory(HSMM, color_by = "Pseudotime")
#5.2.比较细胞分化轨迹进程中功能基因的表达差异

#“%in%”:匹配符号，判断前面向量中的元素是否在“%in%”后面向量中存在，是match()函数等价方式。
to_be_tested <- row.names(subset(fData(HSMM),gene_short_name %in% c("Trem2",)))
to_be_tested
cds_subset <- HSMM[to_be_tested,]
dim(cds_subset)

#主要用到sm.ns()函数根据表达量拟合曲线；
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by = 'seurat_clusters')







genes <- c("Myh11")
p1 <- plot_genes_jitter(HSMM[genes,], grouping = "Pseudotime", color_by = "cellType")
p2 <- plot_genes_violin(HSMM[genes,], grouping = "seurat_clusters", color_by = "seurat_clusters")
p3 <- plot_genes_in_pseudotime(HSMM[genes,], color_by = "cellType",cell_size = 0.65)+scale_color_manual(values=cell_type_cols) 
p3
p1
p2

?plot_genes_jitter
plotc <- p1|p2|p3
plotc
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
###根据伪时间模型寻找差异基因(因为根据时间顺序来找差异基因)
to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% f))
cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(HSMM,cores = 15,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by ="State")


###根据伪时间表达pattern聚类基因

marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% f))

diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters =2,
                        cores = 1,
                        show_rownames = F)

### num_clusters =  可以自己调整


BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 15)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM[1:50,],
                            branch_point = 1,
                            num_clusters = 5,
                            cores = 10,
                            use_gene_short_name = T,
                            show_rownames = T)



save(HSMM,file = "good-HSMM.Rdata")

























