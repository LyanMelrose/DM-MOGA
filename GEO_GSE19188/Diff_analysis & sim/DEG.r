##limma Differential expression

##GSE14407   GSE14407_Gene_standardized.txt  12control,12case   
##GSE54388   GSE54388_Gene_standardized.txt  6control,16case
##GSE66957   GSE66957_Gene_standardized.txt  12control,57case

data<-read.table("D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/GSE19188_Sanger_ENTREZ_ID",sep="\t",header=TRUE,row.names = 1)

data<-na.omit(data)

# 去掉非control和case的样本数据
#对于GSE43696
#data<-data[,-c(21:70),]
#对于GSE31773
#data<-data[,-c(17:20),]
#data<-data[,-c(25:28),]

len<-length(data[1,])


#####limma
library(limma)

exp<-data[,1:len]

exp<-as.matrix(exp)

exp<-scale(exp, center = TRUE, scale = TRUE)

#修改control和case的个数
#GSE31773
#samps<-factor(c(rep("control",16),rep("case",16)))
#GSE43696
#samps<-factor(c(rep("control",20),rep("case",38)))
#GSE18842
#samps<-factor(c(rep("control",45),rep("case",46)))
#GSE19804
#samps<-factor(c(rep("case",60),rep("control",60)))
#GSE19188
samps<-factor(c(rep("case",91),rep("control",65)))

design <- model.matrix(~0+samps)
colnames(design) <- c("case","control") 
fit <- lmFit(exp, design) 
cont.matrix<-makeContrasts(case-control,levels=design) 
fit2 <- contrasts.fit(fit, cont.matrix) 
fit2 <- eBayes(fit2) 
final<-topTable(fit2, coef=1, number=dim(exp)[1], adjust.method="BH", sort.by="B", resort.by="M")


###筛选调整后的p值文件
##logFC>1,logFC<-1,final$adj.P.Val<0.05
##对于GEO的哮喘数据
#pval_0.05_final <- final[which(final$P.Value < 0.05),]
#select_matrix_adjp=data[rownames(adj_pval_0.05_final),]
#select_matrix_p=data[rownames(pval_0.05_final),]
adj_pval_0.05_final <- final[which(final$adj.P.Val < 0.05),]
logFC_1_final<-final[which(final$logFC>1),]
logFC_f1_final<-final[which(final$logFC< -1 ),]


#取出上调基因和下调基因的名字,匹配然后输出文件
up_gene <- intersect(row.names(logFC_1_final),row.names(adj_pval_0.05_final))
down_gene <- intersect(row.names(logFC_f1_final),row.names(adj_pval_0.05_final))
write.table(up_gene,file ="D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/ENid_GSE19188_limma_UPgene.csv", sep=",", row.names = TRUE, col.names = TRUE, quote =FALSE)
write.table(down_gene,file ="D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/ENid_GSE19188_limma_DOWNgene.csv", sep=",", row.names = TRUE, col.names = TRUE, quote =FALSE)
up<-length(up_gene)
down<-length(down_gene)
up_gene_file<-data[up_gene,]
down_gene_file<-data[down_gene,]

#取出校正后p值小于0.05的差异表达基因名,这些基因的表达数据输出
selected_genes<-data[rownames(adj_pval_0.05_final),]
write.table(selected_genes,file ="D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/ENid_GSE19188_limma_genes_expression.csv", sep=",", row.names = TRUE, col.names = TRUE, quote =FALSE)
#导出差异表达的p值等指标
write.table(final[rownames(adj_pval_0.05_final),],file ="D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/ENid_GSE19188_limma_genes_index.csv", sep=",", row.names = TRUE, col.names = TRUE, quote =FALSE)

# 使用GOSemSim计算整个数据集的语义相似性
library(GOSemSim)
# 安装人类数据包
# BiocManager::install("org.Hs.eg.db")
# 读取gene symbol转化为ENTREZ id后的基因列表
# genes<-read.table("E:/documents/_Papers/Paper_3/DATA/GEO/ENTREZid_GSE31773_DEG_matrix.txt",header=TRUE,sep="\t",quote = "")
# 这里还有问题要确认
genes<-c(rownames(selected_genes))
# genes <- c(genes$ID)
Data <- godata(OrgDb = 'org.Hs.eg.db', keytype = "ENTREZID", ont="BP", computeIC = TRUE)
sim<-mgeneSim(
  genes,
  semData=Data,
  measure = "Rel"
)

remain_genes<-rownames(sim)
write.table (sim,file ="D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/ENid_GSE19188_sim_Matrix.csv", sep=",", row.names = FALSE, col.names =FALSE, quote =FALSE)
write.table (remain_genes,file ="D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/ENid_GSE19188_sim_genes.csv", sep=",", row.names = FALSE, col.names =FALSE, quote =FALSE)

#更改selected_genes为remain_genes在data中对应的数据
selected_genes<-data[remain_genes,]
write.table (selected_genes,file ="D:/_Documents/_Papers/Paper_3/DATA/GEO/GEO_GSE19188/ENid_GSE19188_sim_genes_expression.csv", sep=",", row.names = TRUE, col.names = TRUE, quote =FALSE)


#上下调基因limma结果
#up_gene_limma<-final[up_gene,]
#down_gene_limma<-final[down_gene,]

#write.table(up_gene_limma,"E:/documents/_Papers/Paper_3/DATA/GEO/GSE31773_up_gene.csv",row.names=TRUE,col.names=TRUE,quote=FALSE,sep=",")
#write.table(down_gene_limma,"E:/documents/_Papers/Paper_3/DATA/GEO/GSE31773_down_gene.csva",row.names=TRUE,col.names=TRUE,quote=FALSE,sep=",")
