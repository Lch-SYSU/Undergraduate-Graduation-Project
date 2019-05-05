#########################数据整理#########################

# 读入CodonW处理CDS后生成的原始文件
data1 <- read.table(file = 'D:/学习/毕设/数据/hk_gene seq/密码子偏好相关参数1.txt', header = T)
data2 <- read.table(file = 'D:/学习/毕设/数据/hk_gene seq/密码子偏好相关参数2.txt', header = T)

# 重命名data1的基因列名
names(data1)[names(data1) == 'title'] <- c('Gene_description')

# 合并数据集
data3 <- merge(x = data1, y = data2, by = 'Gene_description')

# 将acc num写入新文件
extreme.bias <- subset(x = data3, subset = data3$Nc <= 35, select = Gene_description)
no.bias <- subset(x = data3, subset = data3$Nc >= 60, select = Gene_description)
write.table(x = extreme.bias$Gene_description, file = 'D:\\学习\\毕设\\数据\\hk_gene seq\\ENC小于等于35 acc.txt', quote = F, row.names = F)
write.table(x = no.bias$Gene_description, file = 'D:\\学习\\毕设\\数据\\hk_gene seq\\ENC大于等于60 acc.txt', quote = F, row.names = F)

ENC.less35 <- read.table(file = 'D:\\学习\\毕设\\数据\\hk_gene seq\\ENC小于等于35 ID.txt', header = T)
target_gene_id <- unique(ENC.less35$GeneID)


#########################GO富集#########################
library(clusterProfiler)
library(org.Hs.eg.db)

go_model <- enrichGO(target_gene_id, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,keyType = 'ENTREZID')

###go富集结果barplot图
barplot(go_model,showCategory=20,drop=T)
####go富集结果点图
dotplot(go_model,showCategory=50)
###绘制GO的网络关系图
go.BP <- enrichGO(go_model, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
plotGOgraph(go.BP)
###ont='CC'也可以改为ont='BP'或ont='MF'