# enrich GO based on clusterProfiler for a single module

library(clusterProfiler)

C_ENid <- read.table(
  "D:/_Documents/_Papers/Paper_3/DATA/GEO/clusterProfiler_GO/_ENid_BP_GSE19188_Module_1.txt",
  header = FALSE,
  sep = "\t"
)
C_ENid<-as.vector(unlist(C_ENid))


# GO分析
enrichGO_res <- enrichGO(
  gene=C_ENid,
  OrgDb="org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)

summary(enrichGO_res)
dotplot(enrichGO_res,showCategory=10)
# 在S3对象中，一般我使用$来访问一个对象的属性
# 但在S4对象中，我们只能使用@来访问一个对象的属性
write.table(
  enrichGO_res,
  file = "D:/_Documents/_Papers/Paper_3/DATA/GEO/clusterProfiler_GO/GSE19188_GO.txt",
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