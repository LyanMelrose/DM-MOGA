# install.packages("ggpubr")
# install.packages("ggthemes")

library(ggpubr)
library(ggthemes)

# 读取差异表达基因数据
deg.data <- read.table(
 "D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19804/ENid_GSE19804_limma_genes_index.csv",
 header = T,
 sep = ","
)


# 对差异表达基因的 adj.P.Val 进行 -log10 转换
deg.data$logP <- -log10(deg.data$adj.P.Val)

# (绘制火山图)

deg.data$Group = "not-significant"
# 将adj.P.Val<0.05，logFC>1的基因设为显著上调基因
# 将adj.P.Val<0.05，logFC<-1的基因设为显著上调基因
deg.data$Group[which((deg.data$adj.P.Val < 0.05)&(deg.data$logFC > 1))] = "up-regulated"
deg.data$Group[which((deg.data$adj.P.Val < 0.05)&(deg.data$logFC < -1))] = "down-regulated"
# 继续绘制
ggscatter(
  deg.data,
  x = "logFC",
  y = "logP",
  color = "Group",
  palette = c("#2f5688","#BBBBBB","#CC0000"),
  size = 2) + theme_base()