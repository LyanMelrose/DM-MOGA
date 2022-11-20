if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
BiocManager::install("preprocessCore")

library(WGCNA)
library(R.matlab)
#library(rhdf5)
library( igraph )
if(packageVersion("igraph") < "1.0.0") {
  stop("Need to install igraph version 1.0.0 or above")
}

Mat<-readMat("Network.mat")
corMat<-Mat$Network
rm(Mat)

par(mfrow = c(2, 2))

##软阈值筛选以及画图过程##
powers = c(seq(1,10,by = 1), seq(12, 20, by = 2))
threshold = pickSoftThreshold(
  corMat,
  powerVector = powers,
  verbose = 0)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(threshold$fitIndices[,1], -sign(threshold$fitIndices[,3])*threshold$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(threshold$fitIndices[,1], -sign(threshold$fitIndices[,3])*threshold$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")
plot(threshold$fitIndices[,1], threshold$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(threshold$fitIndices[,1], threshold$fitIndices[,5], labels=powers, cex=cex1,col="red")

#选择6作为软阈值
softPower = 6;
#构建邻接矩阵
adjacency = adjacency(corMat, power = softPower);



##硬阈值选取及画图，corMat需位于[-1,1]
cuts = seq(0, 0.95, by = 0.05)
hard_threshold<-pickHardThreshold.fromSimilarity(
  corMat,
  RsquaredCut = 0.85,
  cutVector = cuts,
)
hard_get_Network<-signumAdjacencyFunction(corMat,hard_threshold$cutEstimate)
plot(hard_threshold$fitIndices[,1], -sign(hard_threshold$fitIndices[,4])*hard_threshold$fitIndices[,3],
     xlab="Hard Threshold",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(hard_threshold$fitIndices[,1], -sign(hard_threshold$fitIndices[,4])*hard_threshold$fitIndices[,3],
     labels=cuts,cex=cex1,col="red");

abline(h=0.90,col="red")
plot(hard_threshold$fitIndices[,1], hard_threshold$fitIndices[,6],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(hard_threshold$fitIndices[,1], hard_threshold$fitIndices[,6], labels=cuts, cex=cex1,col="red")

writeMat(
  "E:/个人数据/zhuxuhui/GSE18842/ENid_BP_GSE18842_HardThreshold_Net.mat",
  Network=hard_get_Network,
  fixNames=TRUE,
  matVersion="5",
  onWrite=NULL,
  verbose=FALSE
)

##############################################################################
#尝试使用clusterprofiler进行GO和KEGG分析
# library(BiocManager)
# BiocManager::install("clusterProfiler")
library(clusterProfiler)

C_ENid <- read.table(
  "D:/_Documents/_Papers/Paper_3/DATA/GEO/clusterProfiler_GO/ENid_BP_GSE18842_Final_gene.txt",
  header = FALSE,
  sep = "\t"
)

#C_ENid <- read.table(
#  "D:/_Documents/_Papers/Paper_3/DATA/GEO/clusterProfiler_GO/ENid_BP_GSE19804_Final_gene_100_50.txt",
#  header = FALSE,
#  sep = "\t"
#)

names(C_ENid)<-c("Entrez","group")
# GO分析
enrichGO_res <- compareCluster(Entrez~group,
                     data=C_ENid,
                     fun="enrichGO",
                     OrgDb="org.Hs.eg.db",
                     ont= "BP",
                     pvalueCutoff=0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
summary(enrichGO_res)
dotplot(enrichGO_res,showCategory=10,includeAll=TRUE)
# 在S3对象中，一般我使用$来访问一个对象的属性
# 但在S4对象中，我们只能使用@来访问一个对象的属性
write.table(
  enrichGO_res@compareClusterResult,
  file = "D:/_Documents/_Papers/Paper_3/DATA/GEO/clusterProfiler_GO/GSE18842_GO.txt",
  append = FALSE,
  quote = TRUE,
  sep = "\t",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = TRUE,
  qmethod = c("escape", "double")
)
# KEGG分析
enrichKEGG_res <- compareCluster(Entrez~group,
                               data=C_ENid,
                               fun="enrichKEGG",
                               pvalueCutoff=0.05,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05)
summary(enrichKEGG_res)
dotplot(enrichKEGG_res,showCategory=10,includeAll=TRUE)
# 在S3对象中，一般我使用$来访问一个对象的属性
# 但在S4对象中，我们只能使用@来访问一个对象的属性
write.table(
   enrichKEGG_res@compareClusterResult,
   file = "E:/_Experiments/MATLAB/CMDPSO/ENid_BP_31773_KEGG.txt",
   append = FALSE,
   quote = TRUE,
   sep = "\t",
   eol = "\n",
   na = "NA",
   dec = ".",
   row.names = FALSE,
   col.names = TRUE,
   qmethod = c("escape", "double")
)

# 使用igraph的模块识别功能处理所得网络
Mat<-readMat("E:/documents/_Papers/Paper_3/DATA/GEO/GSE43696_D_no_0.mat")
Network<-Mat$Network
XXX <- graph_from_adjacency_matrix( Network, weighted=TRUE, mode="undirected")
# 可以更换社团检测算法
Community <- cluster_label_prop(XXX, weights = NA, initial = NULL, fixed = NULL)
Index_of_C <- communities(Community)
Max_C_ID <- which.max(lengths(Index_of_C))
Out_C<-unlist(Index_of_C[Max_C_ID])
writeMat(
  "E:/documents/_Papers/Paper_3/DATA/GEO/igraph_Cres_43696.mat",
  Community=Out_C,
  fixNames=TRUE,
  matVersion="5",
  onWrite=NULL,
  verbose=FALSE
)

#Ctrl+Shift+C注释多行
# out<-data.frame( Network, stringsAsFactors=FALSE )
# write.table(
#   out,
#   file = "E:/documents/_Papers/Paper_3/DATA/GEO/Net_31773.txt",
#   append = FALSE,
#   quote = TRUE,
#   sep = " ",
#   eol = "\n",
#   na = "NA",
#   dec = ".",
#   row.names = FALSE,
#   col.names = FALSE,
#   qmethod = c("escape", "double")
# )

#corMat: a matrix of correlations or other measures of similarity.
#threshold: threshold for connecting nodes: all nodes whose corMat is above the threshold will be connected in the resulting network.