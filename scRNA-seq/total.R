
library(cowplot)
library(dplyr)
library(reshape2)
library(loomR)
library(patchwork)
library(garnett)
library(Seurat)
library("data.table")


data_dir_1<-"~/Desktop/New Zm/SR"
data_dir_2<-"~/Desktop/New Zm/SG"
data_dir_3<-"~/Desktop/New Zm/PR"
data_dir_4<-"~/Desktop/New Zm/PG"


SR<- Read10X(data.dir=data_dir_1)
SG<- Read10X(data.dir=data_dir_2)
PR<- Read10X(data.dir=data_dir_3)
PG<- Read10X(data.dir=data_dir_4)

scRNAlist<-list()  #创建一个可以包含数个S4对象的list（scRNAlist)
scRNAlist[[1]] <- CreateSeuratObject(SR, min.cells = 3, min.features =200)
scRNAlist[[2]] <- CreateSeuratObject(SG, min.cells = 3, min.features =200)
scRNAlist[[3]] <- CreateSeuratObject(PR, min.cells = 3, min.features =200)
scRNAlist[[4]] <- CreateSeuratObject(PG, min.cells = 3, min.features =200)

for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst",nfeatures = 3000)
}
#可以为list当中的seurat对象命名一下#
scRNAlist[[1]]@meta.data$orig.ident="S_B6_Td"
scRNAlist[[2]]@meta.data$orig.ident="S_B6_Zs"
scRNAlist[[3]]@meta.data$orig.ident="S_B6/P_Td"
scRNAlist[[4]]@meta.data$orig.ident="S_B6/P_Zs"

diff.wilcox = FindAllMarkers(scRNA1)
##3 %>% 这是通道函数 起传递左右  可以自己百度深入理解一下，我这里只告诉你他是起传递作用的
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
###将marker基因保存一下
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
###把top10marker基因保存一下
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)





#开始锚定，锚定点越多相当于每个细胞取寻找一定数量的锚定点，也就是花费的时间越多，一般默认2000个#
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,anchor.features = 2000)

scRNA1 <- IntegrateData(anchorset = scRNA.anchors)
save(scRNA,file = "scRNA1.Rda")
#值得注意的是，这里采用DefaultAssay函数来进行Assay的转化#
DefaultAssay(scRNA1) <- "integrated"
scRNA1=ScaleData(scRNA1)   #找好了锚定点后，对锚定点数据进行scaledata#
scRNA1 <- RunPCA(scRNA1, npcs = 30, verbose = T)

#基于PCA空间中的欧氏距离计算 nearest neighbor graph，优化任意两个细胞间的距离权重（输入上一步得到的 PC 维数）； 
SS <- FindNeighbors(SS, dims = 1:12)
#接着优化模型，resolution 参数决定下游聚类分析得到的分群数，对于3K左右的细胞，设为0.4-1.2能得到较好的结果(官方说明)；如果数据量增大，该参数也应该适当增大； 
SS <- FindClusters(SS, resolution = 1) 
#使用 Idents（）函数可查看不同细胞的分群； 查看前8个细胞的分群ID
head(Idents(scRNA1), 5) 




scRNA1 <- FindNeighbors(scRNA1, reduction = "pca", dims = 1:5)
scRNA1 <- FindClusters(scRNA1, resolution = 0.4)
head(Idents(scRNA1), 5) 
scRNA1 <- RunUMAP(scRNA1, reduction = "pca", dims = 1:8)
colnames(scRNA1@meta.data)
DimPlot(scRNA1, reduction = "umap", group.by = "cellType")
DimPlot(scRNA1, reduction = "umap", label = FALSE,pt.size = 0.1,split.by = "orig.ident")

rm(scRNA)
DefaultAssay(scRNA1) <- "RNA"        #后续看表达情况的话，比如小提琴、气泡图是用RNA数据来做#，前面的分cluster是按照之前分的来#
scRNA <- ScaleData(scRNA1)
colnames(scRNA)
DotPlot(object = scRNA, features = "Cdh5")
P1
DotPlot(object = scRNA, features = genes1)
genes1<-c("Acta2","Cnn1","Myh11")
genes2<-c("Col3a1","Fn1","Fsp1","Tnc","Vim")
genes3<-c("Vim","Lgals3")
genes4<-c("Col1a1","Col1a2")
genes5<-c("Hist1h2ap","Hist1h1b","Mki67")
genes6<-c("Cd68","Lgals3","Tnf")

genes<-c("Myh11","Acta2","Tagln","Dcn","Slc6a12","Hspa4l","Manf","Col1a2","Col1a1","Col3a1","Col5a2","Fn1","Tm4sf1","Ly6a","Tnf","Cxcl2","Vim","Hist1h2ap","Hist1h1b","Mki67","Spp1","Lyz2","Tgfb1","Cd68","Ctsk","Ccl3","Fcer1g","Tnfrsf11b","Ctsd","Ccl5","Cox17")
genes<-c("Cd68","Lgals3","Fcer1g","Tgfb1","Lyz2","Spp1","Ctss","Slc6a12","Hspa4l","Manf","Col1a1","Col3a1","Tnc","Vim","Cdk1","Pclaf","Cxcl2","Fabp4","Ccl3","Mki67","Tnf","Klf4","Tm4sf1","Hist1h1b","Col1a2","Fn1","Myh11","Acta2","Tagln","Dcn","Pecam1","Cdh5","Ly6a")

genes7<-c("Tnf")

DotPlot(scRNA, features = genes6)+coord_flip()+theme_bw()+#去除背景，旋转图片  
  theme(panel.grid = element_blank(),  
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
  scale_color_gradientn(values = seq(0,1,0.2),colours = c(cols = c("Blue","Orangered")))+#颜色渐变设置  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 









DotPlot(scRNA, features = "Cxcl2",cols = c("Blue","Orangered")) 
p

mycolor <- c('lightgrey', 'orange','red')#设置颜色  
FeaturePlot(seurat_object, features = "Bcl2", pt.size = 0.65,cols = mycolor)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 


FeaturePlot(seurat_object, features = 'Sox9', pt.size = 1.5,cols = mycolor)
VlnPlot(scRNA,features = "Stmn1")


?DotPlot
new.cluster.ids <- c("SMC", "Padi4(hi) Macrophage", "Sensecent SMCs", "Stem-Like SMC", "Sensecent SMCs", "Stem-Like SMC", "Fibroblast II", "Fibrochondrocyte", "Fibromyocyte","Stem-Like SMC","Macrophage I","SMC","Endothelial") #自定义名称
names(new.cluster.ids) 

levels(SS)
#将seurat_object的水平属性赋值给new.cluster.ids的names属性； 
names(new.cluster.ids) <- levels(scRNA1)
names(new.cluster.ids) 
scRNA1 <- RenameIdents(scRNA1, new.cluster.ids) 
scRNA1$cellType=Idents(scRNA1)
Idents(seurat_object)=seurat_object$class
seurat_object <- RenameIdents(seurat_object, new.cluster.ids) 
#绘制 tsne 图(修改标签后的)； 
sample<-scRNA1[scRNA1$orig.ident(S_B6_Td),]

#######自定义细胞分群颜色
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set2"), "#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")


DimPlot(scRNA1, reduction = "umap",label = F,   
                cols= cell_type_cols, 
                pt.size = 0.65,
                repel = T,split.by="orig.ident")

save(Bad,file = "Bad.Rda")

tsneplot2<-UMAPPlot(scRNA1,label = TRUE, pt.size = 0.8) 
tsneplot2

DotPlot(seurat_object,features = "Pecam1")


T1<-subset(scRNA1,idents=c("Padi4(hi) Macrophage","Macrophage I","Macrophage II", "SMC","Stem-Like SMC")) 


P1<-VlnPlot(object = T1, features = "Padi4")
P1

genes<-c("Padi2","Cd68","Myh11","Tm4sf1","Spp1","Timp1")
########丝滑版小提琴图展示######
install.packages("remotes")  
remotes::install_github("lyc-1995/MySeuratWrappers")#通过链接安装包  
library(MySeuratWrappers) 
library(ggplot2)
#需要展示的基因  
markers <- c('CD3D', 'S100A8', 'S100A9', 'CD79A', 'CCL5', 'NKG7', 'GZMA', 'IL32', 'CD4', 'CD8A', 'LTB', 'FCN1', 'MS4A1', 'SPON2','FCER1A','SERPINF1', 'TMEM40', 'CD3E')  
my36colors <-c("#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")#颜色设置  


mycolor<-c("#F18D00","#4D4398")
c<-duplicated(genes)
c
VlnPlot(scRNA, features = genes,  
stacked=T,pt.size=0,  
cols = cell_type_cols,#颜色  
direction = "horizontal", #水平作图  
x.lab = '', y.lab = '')+#横纵轴不标记任何东西  
theme(axis.text.x = element_blank(),   
      axis.ticks.x = element_blank())#不显示坐标刻度 
VlnPlot(seurat_object,features = "Acta2",group.by = orig.ident)

VlnPlot(seurat_object, features = "Caspas3",  
        stacked=T,pt.size=0.3,  
        cols = mycolor,#颜色  
         #水平作图  
        x.lab = '', y.lab = '',group.by = "orig.ident")+#横纵轴不标记任何东西  
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度 

#######RidgePlot#######
RidgePlot(scRNA,features = "C3",ncol = 1)
?RidgePlot
######### Cell Ratio作图#####
library("plyr")
library(ggplot2)
library(gplots)
samples.cellType.stat=as.data.frame(table(seurat_object$cellType,seurat_object$orig.ident))
samples.cellType.stat.prop=ddply(samples.cellType.stat,"Var2",transform,Ratio=Freq/sum(Freq))
head(samples.cellType.stat.prop)
samples.cellType.stat.prop
ggplot(samples.cellType.stat.prop,aes(x=Var2,y=Ratio,fill=Var1))+geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = cell_type_cols)+theme_bw()



########Dotplot2 作图（改进版）######

markers<-c("Padi4","Padi3","Padi2","Padi1")
DotPlot(scRNA, features = markers)+coord_flip()+theme_bw()+#去除背景，旋转图片  
   theme(panel.grid = element_blank(),  
                    axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
    scale_color_gradientn(values = seq(0,1,0.2),colours = c(cols = c("Blue","Orangered")))+#颜色渐变设置  
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

?DotPlot

FeaturePlot(T1,features = "Ly6a")

?FeaturePlot
c<-scRNA1$orig.ident
c
stem<-subset(seurat_object,idents=c("Stem-Like SMC","Fibrochondrocyte")) 


M2<-subset(T1,idents=c("Macrophage II")) 
Macrophage<-subset(seurat_object,idents=c("Padi4(hi) Macrophage","Macrophage"))
Sttem<-subset(seurat_object,idents = c("Stem-Like SMC"))
macrophage<-merge(M1,M2)
macrophage<-merge(macrophage,Mp)
Fibromyocyte<-subset(T1,idents = c("Fibromocyte"))
Fibrochondrocyte<-subset(T1,idents = c("Fibrochondrocyte"))
T1$cellType=Idents(T1)
VlnPlot(M2, features = genes,  
        stacked=T,pt.size=0,  
        cols = cell_type_cols,#颜色  
        direction = "horizontal",
        x.lab = '', y.lab = '',split.by = "orig.ident")+  
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度 

genes<-c("Tnf","Spp1","Tgfb1","Cd68","Mki67","Vim","Padi2","Gsdmd")
genes2<-c("Fn1","Tnfrsf11b","Col1a2")
genes3<-c("Mmp12","Fabp5","Lyve1","Mafb")
cell_type_cols2<-c(brewer.pal(2,"Set3"),"#4D4398","#F18D00")
cell_type_cols2<-c("#F18D00","#4D4398")
VlnPlot(subset(seurat_object,Tnf>0), features = "Tnf",  
        pt.size=1,  
        #颜色  
        cols = cell_type_cols2,
        group.by = "orig.ident")
VlnPlot(seurat_object,features = "Acta2",group.by = "orig.ident")
?VlnPlot
rm(scRNA)

####亚群提取####
Sample<-scRNA1[,scRNA1@meta.data$orig.ident %in% c("S_B6_Td","S_B6/P_Td")]
seurat_object<-Sample
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
seurat_object <- ScaleData(seurat_object, features = all.genes) ##耗时2min
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
seurat_object <- FindClusters(seurat_object, resolution = 0.7) 
#使用 Idents（）函数可查看不同细胞的分群； 查看前8个细胞的分群ID
head(Idents(seurat_object), 5) 


#UMAP非线性降维
seurat_object <- RunUMAP(seurat_object, dims = 1:10, label = T)
head(seurat_object@reductions$umap@cell.embeddings)

umapplot <- DimPlot(seurat_object, reduction = "umap",pt.size = 1.5,label = T,split.by = "orig.ident")
umapplot


diff.wilcox = FindAllMarkers(seurat_object)
##3 %>% 这是通道函数 起传递左右  可以自己百度深入理解一下，我这里只告诉你他是起传递作用的
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
###将marker基因保存一下
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
###把top10marker基因保存一下
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)


DotPlot(seurat_object, features = "Spp1",cols = c("Blue","Orangered")) 

VlnPlot(seurat_object,features = "Acta2")




#####巨噬细胞样转分化######

T1<-subset(Sample,idents=c("Padi4(hi) Macrophage","Macrophage I","Macrophage II","Stem-Like SMC")) 

new.cluster.ids <- c("Macrophage", "Fibroblast II", "Stem-Like SMC", "Fibromyocyte", "Fibroblast I", "Padi4(hi) Macrophage", "Stem-Like SMC", "Fibrochondrocyte", "Padi4(hi) Macrophage","SMC","Fibroblast I","Senescent SMC","Endothelial Cell","Fibrochondrocyte","Fibroblast II","Stem-Like SMC","Macrophage")

names(new.cluster.ids) <- levels(seurat_object)
names(new.cluster.ids) 
seurat_object <- RenameIdents(seurat_object, new.cluster.ids) 
seurat_object$cellType=Idents(seurat_object)
Idents(scRNA1)=scRNA1$integrated_snn_res.0.4
seurat_object <- RenameIdents(seurat_object, new.cluster.ids) 

DimPlot(seurat_object, reduction = "umap",label = F,   
        cols= cell_type_cols, 
        pt.size = 0.65,
        repel = T,split.by = "orig.ident")

DimPlot(scRNA1,group.by = "orig.ident",pt.size = 0.1)
save(seurat_object,file = "td.rds")

######### Cell Ratio作图#####
library("plyr")
library(ggplot2)
library(gplots)
samples.cellType.stat=as.data.frame(table(scRNA1$cellType,scRNA1$orig.ident))
samples.cellType.stat.prop=ddply(samples.cellType.stat,"Var2",transform,Ratio=Freq/sum(Freq))
head(samples.cellType.stat.prop)
samples.cellType.stat.prop
ggplot(samples.cellType.stat.prop,aes(x=Var2,y=Ratio,fill=Var1))+geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = cell_type_cols)+theme_bw()


########拟时序分析

Bad<-subset(seurat_object,idents=c("Macrophage","Padi4(hi) Macrophage","Stem-Like SMC","Senescent SMC","Fibroblast I")) 



#####good样细胞转分化######


Good<-subset(seurat_object,idents=c("Stem-Like SMC","Fibromyocyte","Fibrochondrocyte","Fibroblast II","Endothelial Cell")) 



######Control与Case组的差异分析：


deg[,6]<-rownames(deg)
plot(deg$avg_log2FC,(deg$pct.1-deg$pct.2))
library(EnhancedVolcano)
install.packages("EnhancedVolcano")
BiocManager::install("prog4")

install.packages("prog4")

library("prog4")

R.version


rm(scRNA1,seurat_object)


