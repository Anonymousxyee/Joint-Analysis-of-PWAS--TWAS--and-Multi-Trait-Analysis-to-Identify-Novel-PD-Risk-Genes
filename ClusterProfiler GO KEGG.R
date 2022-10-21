rm(list=ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr) 
library(ggplot2)
ID=read.csv("CNS.csv")

CNS1 = bitr(ID$id, fromType = "SYMBOL", 
            toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")

head(CNS1)

goCNS <- enrichGO(CNS1$ENTREZID,OrgDb = org.Hs.eg.db,
                  ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,keyType = 'ENTREZID',
                  readable =T)

dim(goCNS)

write.csv(goCNS,file="goCNS.csv")

##########################################

barplot(goCNS)

dotplot(goCNS)

cnetplot(goCNS)

kegg <- enrichKEGG(CNS1$ENTREZID, organism = 'hsa', 
                   keyType = 'kegg', pvalueCutoff = 0.05, 
                   pAdjustMethod = 'BH', minGSSize = 3, 
                   maxGSSize = 500, qvalueCutoff = 0.2, 
                   use_internal_data = FALSE)

dim(kegg)

keggCNS<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


write.csv(keggCNS,file = "keggCNS.csv")


barplot(keggCNS)

dotplot(keggCNS)

cnetplot(keggCNS)
