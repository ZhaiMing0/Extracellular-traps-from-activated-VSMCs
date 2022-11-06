Sys.setenv(LANGUAGE = "en")
BiocManager::install('DESeq2')
library("DESeq2")
counts<-read.table("counts.txt",sep = "\t",header = T)
row.names(counts)=counts[,1]
counts=counts[,-1]
# 预处理，过滤低丰度的数据
countData <- counts[apply(counts, 1, sum) > 0 , ]
# 读取样本分组信息
colData <- read.table(
  "sample_group.txt",
  header=T,
  sep="\t")
# 构建DESeq2中的对象
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData, ~Group)
# 指定哪一组作为control
dds$Group <- relevel(dds$Group, ref = "Ctr")
#####计算每个样本的归一化系数####
dds <- estimateSizeFactors(dds)
#####通过下列代码估算每个基因的a值######
dds <- estimateDispersions(dds)
#######差异分析
dds <- nbinomWaldTest(dds)
res <- results(dds)

dds <- DESeq(dds)
res <- results(dds)
write.table(
  res,
  "DES12.27.txt",
  sep="\t",
  quote=F,
  col.names = NA)

library(EnhancedVolcano)

vals<-c('socs1','Tlr4','Myd88','Tbk1','Tmem173',"Col1a1","Col1a2","Tagln","Acta2","Akt","Tgfb")
vals<-c("socs1")
group<-ifelse(
  data$log2FoldChange<(-0.5)&data$pvalue<0.05,'#4D4398',
  ifelse(data$log2FoldChange>(0.5)&data$pvalue<0.05,'#F18D00',
         '#b5b5b5'))
group[is.na(group)]<-'#b5b5b5'
names(group)[group=='#F18D00']<-'Up'
names(group)[group=='#b5b5b5']<-'Nodiff'
names(group)[group=='#4D4398']<-'Down'

data<-read.table("DES12.27.txt",sep = "\t",header = T)
row.names(data)=data[,1]
data=data[,-1]

data[,7]<-rownames(data)

EnhancedVolcano(data,
                x="log2FoldChange",
                y="pvalue",
                lab=data$V7,
                pCutoff=10e-1/20,#y轴阈值线(水平)
                FCcutoff=0.5,#x轴阈值线（垂直）
                labSize=3.5,#标签大小
                xlim=c(-6, 6),#限制X轴范围
                ylim=c(0,10),#限制Y轴范围
                selectLab=c(vals),#使用selectLab参数选定所关注的标签
                xlab=bquote(~Log[2]~'fold change'),#将内容传递给xlab
                labCol='black',#标签颜色
                labFace='bold',#标签字体
                boxedLabels=TRUE,#是否在框中绘制标签
                drawConnectors=TRUE,#是否通过连线将标签连接到对应的点上
                widthConnectors=0.8,#连线的宽度
                endsConnectors="last",#连线绘制箭头的方向，可选first、both、last
                colConnectors='black',#连线的颜色
                colCustom=group,#用group覆盖默认配色方案
                colAlpha=0.6,#调整透明度
                cutoffLineType='longdash',#阈值线类型，可选“blank”、“solid”、“dashed”、“dotted”、“dotdash”、“longdash”和“twodash”
                cutoffLineCol='pink',#阈值线颜色
                cutoffLineWidth=0.88,#阈值线粗细
                title="Volcano Plot Exp",#主标题
                subtitle="Differential expression",#副标题
                caption=bquote(~Log[2]~"fold change cutoff,1;p-value cutoff,0.05"),#注释说明
                legendPosition='right',#图例位置
                legendLabSize=12,#图例文字大小
                legendIconSize=6,#图例符号大小
                ####↓新加入↓####
                #通过ifelse条件判断函数指定散点大小


pointSize=c(ifelse(data$V7 %in% vals,5,3)))+coord_flip()













