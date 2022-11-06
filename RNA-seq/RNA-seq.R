library(pheatmap)
library(ggplot2)
Sys.setenv(LANGUAGE = "en")
data2 <- read.table("gene name.txt",  # 读取的数据文件名称，这里文件是放在工作目录下
                    header=T, # 数据集第一行为变量名
                    row.names=1, # 第一列为行名
                    sep="\t") # 指定分隔符号
test <- fpkm[apply(fpkm, 1, function(x) sd(x)!=0),]
pheatmap(test,scale = "row")

data3 <- read.table("DEGs_heat.txt",  # 读取的数据文件名称，这里文件是放在工作目录下
                    header=T, # 数据集第一行为变量名
                    row.names=1, # 第一列为行名
                    sep="\t") # 指定分隔符号
############在对行进行数据放缩的时候，要先通过下列一步把有0/NA值的部分祛除#####
test1 <- data2[apply(data2, 1, function(x) sd(x)!=0),]
pheatmap(test1,scale = "row",show_rownames = T,cluster_rows = F,cluster_cols = F,cellwidth = 30,cellheight = 21,annotation_row = group)

group <- read.table("group.txt",  # 读取的数据文件名称，这里文件是放在工作目录下
                    header=T, # 数据集第一行为变量名
                    row.names=1, # 第一列为行名
                    sep="\t") # 指定分隔符号
