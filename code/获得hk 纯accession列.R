hkgene <- read.table('D://学习/毕设/文献/数据/HK_genes.txt')
hkacc <- hkgene[,2]
write.table(x = hkacc, file = 'D://学习/毕设/文献/数据/HK_acc.txt', row.names = F, quote = F)
