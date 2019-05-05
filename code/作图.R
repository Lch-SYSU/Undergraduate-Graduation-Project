#########################数据整理#########################
library('ggplot2')
library('reshape2')
library('corrplot')
library('readr')
library('ggpubr')

# 读入CodonW处理CDS后生成的原始文件
data1 <- read.table(file = 'D:/学习/毕设/数据/data/全基因组/密码子偏好相关参数1.txt', header = T)
data2 <- read.table(file = 'D:/学习/毕设/数据/data/全基因组/密码子偏好相关参数2.txt', header = T)

# 合并数据集
names(data1)[names(data1) == 'title'] <- c('Gene_description')
data3 <- merge(x = data1, y = data2, by = 'Gene_description')

# 添加列GC12、expENC、CAI
data3 <- transform(data3, GC12 = (GC1 + GC2)/2)

data3 <- transform(data3, expENC = 2 + GC3s.x + (29/(GC3s.x**2 + (1-GC3s.x)**2)))
data3 <- transform(data3, expProp = (expENC-Nc)/expENC)

data3 <- subset(data3, select = -CAI)
new_cai <- read.table(file = 'D:/学习/毕设/数据/data/全基因组/CAI.txt', header = T) 
names(new_cai)[names(new_cai) == 'FastaHeader'] <- c('Gene_description')
data3 <- merge(x = data3, y = new_cai, by = 'Gene_description')


#########################GC分布#########################

data.gc <- subset(data3, select = c(Gene_description, GC1, GC2, GC3, GC3s.x, GC.x))
names(data.gc)[names(data.gc) == 'GC.x'] <- c('GC.all')
names(data.gc)[names(data.gc) == 'GC3s.x'] <- c('GC3s')
data.gc2 <- melt(data.gc, id.vars = 'Gene_description', variable.name = 'GC.type', value.name = 'GC.value')

ggplot(data.gc2, aes(x = factor(GC.type), y = GC.value)) + 
  geom_boxplot(outlier.size = 0.1, width = 0.5) +
  theme(axis.line = element_line(colour = 'black')) +
  stat_summary(fun.y = 'mean', geom = 'point', fill = 'white', shape = 23, size = 1.5) +
  xlab(NULL)


#########################ENC分布#########################

ggplot(data3, aes(x = 1, y = Nc)) + 
  geom_violin() +
  geom_boxplot(width = 0.1, fill = 'black', outlier.colour = NA) +
  scale_x_continuous(breaks = NULL) +
  theme(axis.title.x = element_blank()) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65))+
  stat_summary(fun.y = 'mean', geom = 'point', fill = 'white', shape = 23, size = 2)


#########################ENC-plot#########################

# 自定义标准曲线函数
st.fun <- function(gc3){
  2 + gc3 + 29/(gc3**2 + (1-gc3)**2)
}

# 散点图 + 标准曲线
ggplot(data3, aes(x = GC3s.x, y = Nc)) + 
  geom_point(shape = 16, size = 0.1) + 
  stat_function(fun = st.fun, size = 1) +
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65)) +
  theme(axis.line = element_line(colour = 'black')) +
  xlab('GC3s')

# 预期值与实际值的秩和检验
shapiro.test(x = data3$Nc)
shapiro.test(x = data3$expENC)
f <- factor(rep(c('enc','exp'),c(length(data3$Nc),length(data3$expENC))))
bartlett.test(x = c(data3$Nc, data3$expENC), g = f)
wilcox.test(data3$Nc, data3$expENC, paired = F, conf.int = T)

#########################中性图#########################

ggplot(data3, aes(x = GC3s.x, y = GC12)) + 
  geom_point(shape = 16, size = 0.1) + 
  stat_smooth(method = lm) +
  theme(axis.line = element_line(colour = 'black')) +
  annotate('text', label = 'y=0.19x+0.39, R^2=0.32', x=0.4, y=0.7) +
  xlab('GC3s')

cor.test(x = data3$GC12, y = data3$GC3s.x, method = 'spearman', exact = F)
lm1 <- lm(formula = GC12 ~ 1 + GC3s.x, data = data3)
summary(lm1)

#########################对应分析#########################

# 每轴贡献率
inertia <- read.table(file = 'D:\\学习\\毕设\\数据\\data\\对应分析\\Correspondence analysis results from [CodonW] on data 2\\前40轴相对和累积贡献率.txt', header = T)
ggplot(inertia, aes(x = factor(Num.), y = R.Iner.)) +
  theme(axis.line = element_line(colour = 'black')) +
  xlab('轴') + ylab('相对贡献率') +
  theme(text = element_text(family = 'SimSun')) +
  geom_bar(stat = 'identity', colour = 'black', fill = 'red')

# 密码子前两轴分布
codon.coa <- read.csv(file = 'D:\\学习\\毕设\\数据\\data\\对应分析\\Correspondence analysis results from [CodonW] on data 2\\每密码子前四轴分布2.csv', header = T)
ggplot(codon.coa, aes(x = Axis1, y = Axis2)) +
  geom_point()+
  geom_text(aes(label = label), size = 1.5, vjust = -1)

# 密码子第三位碱基的一二轴分布
codonend.coa <- read.csv(file = 'D:\\学习\\毕设\\数据\\data\\对应分析\\Correspondence analysis results from [CodonW] on data 2\\每密码子前四轴分布2.csv', header = T)
ggplot(codonend.coa, aes(x = Axis1, y = Axis2, shape = factor(type))) +
  theme(axis.line = element_line(colour = 'black')) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(size = 2, colour = 'blue') + 
  scale_shape_discrete(name = '密码子第三位碱基', labels = c('A','G','C','U')) +
  theme(text = element_text(family = 'SimSun'))

# 汇总COA数据
cds.coa <- read.table(file = 'D:\\学习\\毕设\\数据\\data\\对应分析\\Correspondence analysis results from [CodonW] on data 2\\每个CDS 4轴.txt', header = T)
names(cds.coa)[names(cds.coa) == 'label'] <- c('Gene_description')
cds.coa2 <- merge(x = cds.coa, y = data3, by = 'Gene_description')

coa.corplot <- subset(cds.coa2, select = c(GC3s.x, GC.x, T3s, C3s, A3s, G3s, Axis1, Nc, Len_aa, CAI))
names(coa.corplot)[names(coa.corplot) == 'GC3s.x'] <- c('GC3s')
names(coa.corplot)[names(coa.corplot) == 'GC.x'] <- c('GC.all')

# 相关矩阵
cor <- cor(coa.corplot)
corrplot(cor, method = 'shade',shade.col = NA, tl.col = 'black', tl.srt = 45, order = 'AOE')

# 每个密码子沿第一轴的GC分布
cds.axis <- read.table('D:\\学习\\毕设\\数据\\data\\对应分析\\Correspondence analysis results from [CodonW] on data 2\\每个CDS 4轴.txt', header = T)
names(cds.axis) <- c('Gene_description', 'Axis1', 'Axis2', 'Axis3', 'Axis4')
cds.axis2 <- merge(x = data3, y = cds.axis, by = 'Gene_description') 

ggplot(cds.axis2, aes(x = Axis1, y = Axis2, colour = -GC3s.x)) +
  geom_point() +
  theme(axis.line = element_line(colour = 'black')) +
  labs(colour = 'GC3s')

#########################最优密码子#########################

# 根据ENC值排序data3，取出首尾10%基因名
data3 <- data3[order(data3$Nc),]
tail.enc <- as.character(head(x = data3$Gene_description, n = nrow(data3) * 0.1))
head.enc <- as.character(tail(x = data3$Gene_description, n = nrow(data3) * 0.1))

write.table(x = head.enc, file = 'D:/学习/毕设/数据/data/极端ENC/ENC极端高 基因.txt', row.names = F, quote = F)
write.table(x = tail.enc, file = 'D:/学习/毕设/数据/data/极端ENC/ENC极端低 基因.txt', row.names = F, quote = F)


#########################获取低ENC高GC基因#########################

Lowenc.Lowgc <- subset(x = data3, subset = GC3s.x<0.3 & Nc<49, select = c(Gene_description, GC3s.x, Nc))
write.table(x = Lowenc.Lowgc$Gene_description, file = 'D:/学习/毕设/数据/data/表达量/GC3s.x0.3 Nc49 acc.txt', quote = F, row.names = F)

Lowenc.Highgc <- subset(x = data3, subset = GC3s.x>0.8 & Nc<36, select = c(Gene_description, GC3s.x, Nc))
write.table(x = Lowenc.Highgc$Gene_description, file = 'D:/学习/毕设/数据/data/表达量/GC3s.x0.8 Nc36 acc.txt', quote = F, row.names = F)

medium.gc <- subset(x = data3, subset = GC3s.x>0.45 & GC3s.x<0.55 & Nc>58,  select = c(Gene_description, GC3s.x, Nc))
write.table(x = medium.gc$Gene_description, file = 'D:/学习/毕设/数据/data/表达量/GC3s.x0.45-0.55 Nc58 acc.txt', quote = F, row.names = F)


#########################不同GC下的表达量比较#########################

# 读入表达量数据1
data1.gc0.3 <- read.csv('D:\\学习\\毕设\\数据\\data\\表达量\\GC0.3 TPM SRR3589956.csv', header = T)
data1.gc0.3[is.na(data1.gc0.3)] = 0
data1.gc0.8 <- read.csv('D:\\学习\\毕设\\数据\\data\\表达量\\GC0.8 TPM SRR3589956.csv', header = T)
data1.gc0.8[is.na(data1.gc0.8)] = 0
data1.gc0.5 <- read.csv('D:\\学习\\毕设\\数据\\data\\表达量\\GC0.5 TPM SRR3589956.csv', header = T)
data1.gc0.5[is.na(data1.gc0.5)] = 0

# 读入表达量数据2
data2.expr <- read_csv('D:/学习/毕设/数据/data/表达量/other rna_celline.csv', col_names = T)

gc0.3.ENSG <- read_tsv('D:/学习/毕设/数据/data/表达量/GC0.3 ENSG.txt', col_names = T)
names(gc0.3.ENSG) <- c('Accession', 'Gene')
whole1.0.3.TPM <- merge(x = data2.expr, y = gc0.3.ENSG, by = 'Gene')
whole2.0.3.TPM <- dcast(data = whole1.0.3.TPM, formula = Gene ~ Sample, value.var = 'Value')

gc0.5.ENSG <- read_tsv('D:/学习/毕设/数据/data/表达量/GC0.5 ENSG.txt', col_names = T)
names(gc0.5.ENSG) <- c('Accession', 'Gene')
whole1.0.5.TPM <- merge(x = data2.expr, y = gc0.5.ENSG, by = 'Gene')
whole2.0.5.TPM <- dcast(data = whole1.0.5.TPM, formula = Gene ~ Sample, value.var = 'Value')

gc0.8.ENSG <- read_tsv('D:/学习/毕设/数据/data/表达量/GC0.8 ENSG.txt', col_names = T)
names(gc0.8.ENSG) <- c('Accession', 'Gene')
whole1.0.8.TPM <- merge(x = data2.expr, y = gc0.8.ENSG, by = 'Gene')
whole2.0.8.TPM <- dcast(data = whole1.0.8.TPM, formula = Gene ~ Sample, value.var = 'Value')

# 箱线图
boxplot(log2(data1.gc0.3$TPM1+1), log2(data1.gc0.5$TPM1+1), log2(data1.gc0.8$TPM1+1), notch = T)
boxplot(log2(whole2.0.3.TPM$`HMC-1`), log2(whole2.0.5.TPM$'HMC-1'), log2(whole2.0.8.TPM$'HMC-1'), notch = T, main = 'HMC-1', names = c('GC0.3', 'GC0.5', 'GC0.8'))

df1 <- data.frame(expr = c(log2(whole2.0.3.TPM$HeLa), log2(whole2.0.5.TPM$HeLa), log2(whole2.0.8.TPM$HeLa)), gctype = rep(c('GC0.3','GC0.5','GC0.8'),c(110,100,110)))
compar <- list(c('GC0.3','GC0.5'),c('GC0.3','GC0.8'),c('GC0.5','GC0.8'))
ggboxplot(df1, x = 'gctype', y = 'expr', title = 'Hela', xlab = 'GC含量', ylab = '表达量/log2(TPM)')



# 56细胞系表达量差异
  # 正态分布检验
Pvalue1 <- c(rep(0,ncol(whole2.0.3.TPM)-1))
for (j in 2:ncol(whole2.0.3.TPM)){
  Pvalue1[j-1] <- shapiro.test(x = log2(whole2.0.3.TPM[,j]))$W
}

  # 方差分析
Pvalue3 <- c(rep(0,ncol(whole2.0.3.TPM)-1))

  # Kruskal-Wallis秩和检验
Pvalue2 <- c(rep(0,3*ncol(whole2.0.3.TPM)-1))
f1 <- rep(1:2, c(110,100))
f2 <- rep(2:3, c(100,110))
for (i in 2:57){
  pwil.test1 <- pairwise.wilcox.test(x = c(whole2.0.3.TPM[,i], whole2.0.5.TPM[,i]), g = f1, p.adjust.method = 'fdr')
  pwil.test2 <- pairwise.wilcox.test(x = c(whole2.0.5.TPM[,i], whole2.0.8.TPM[,i]), g = f2, p.adjust.method = 'fdr')
  Pvalue2[i-1] <- as.vector(pwil.test1$p.value)
  Pvalue2[i] <- as.vector(pwil.test2$p.value)
  
}

# t检验
  # 正态性检验
shapiro.test(log2(data1.gc0.3$TPM1+1))
shapiro.test(log2(data1.gc0.8$TPM1+1))

  # 方差齐性检验
len1 <- sum(complete.cases(data1.gc0.3$TPM1))
len2 <- sum(complete.cases(data1.gc0.8$TPM1))
f <- factor(rep(c("低GC","高GC"), c(len1, len2)))
bartlett.test(x = c(log2(data1.gc0.3$TPM1+1), log2(data1.gc0.8$TPM1+1)), g = f)

  # t检验与Wilcox非参数检验
t.test(x = log2(data1.gc0.3$TPM1+1), y = log2(data1.gc0.8$TPM1+1), var.equal = T)
wilcox.test(x = data1.gc0.3$TPM1, y = data1.gc0.8$TPM1, paired = F, conf.int = T)


#########################不同GC下的氨基酸使用情况比较#########################

aa.usage <- read.table('D:\\学习\\毕设\\数据\\data\\表达量\\result\\relative amino acid usage.txt', header = T)
aa.usage <- aa.usage[order(aa.usage[1,])] 
aa.usage2 <- melt(data = aa.usage, id.vars = 'Genename', variable.name = 'aa', value.name = 'ratio')

ggplot(aa.usage2, aes(x = aa, y = ratio, fill = Genename)) +
  geom_bar(position = 'dodge', stat = 'identity')


#########################不同GC下的密码子使用情况比较#########################

cd.usage1 <- read.csv('D:\\学习\\毕设\\数据\\data\\表达量\\result\\gc0.3 0.5 0.8 RSCU.csv', header = T, fileEncoding = 'gbk')
cd.usage2 <- melt(cd.usage1, id.vars = 'aa.cd', measure.vars = c('X0.3RSCU', 'X0.5RSCU', 'X0.8RSCU', 'hkCDS.RSCU'), variable.name = 'GC.type', value.name = 'RSCU')

# 分面作图。将数据排序，分割为两半，分别作图
cd.usage2 <- cd.usage2[order(cd.usage2$aa.cd),]
cd.usage2.1 <- cd.usage2[1:56,]
cd.usage2.2 <- cd.usage2[57:116,]
cd.usage2.3 <- cd.usage2[117:172,]
cd.usage2.4 <- cd.usage2[173:236,]

ggplot(cd.usage2, aes(x = GC.type, y = RSCU, fill = GC.type)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~aa.cd) +
  theme(axis.text.x = element_blank()) +
  scale_fill_discrete(labels = c('GC0.3', 'GC0.5', 'GC0.8', 'whole HK gene'))

ggplot(cd.usage2.4, aes(x = GC.type, y = RSCU, fill = GC.type)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~aa.cd) +
  theme(axis.text.x = element_blank()) +
  scale_fill_discrete(labels = c('GC0.3', 'GC0.5', 'GC0.8', 'whole HK gene'))

#########################重新审视管家基因列表#########################

# 从组织表达数据入手
tissue.expr1 <- read_csv('D:\\学习\\毕设\\数据\\data\\表达量\\other rna_tissue.csv', col_names = T)
tissue.expr2 <- dcast(data = tissue.expr1, formula = Gene ~ Sample, value.var = 'Value')
  # 标准1：删除0表达的Gene
tissue.expr2[tissue.expr2 == 0] <- NA
tissue.expr2 <- tissue.expr2[complete.cases(tissue.expr2),]

  # 标准2：使ENSG为行名。取log2，保留标准差小于1的Gene。
row.names(tissue.expr2) <- tissue.expr2[,1]
tissue.expr2 <- tissue.expr2[,-1]

tissue.expr2 <- log2(tissue.expr2)
tissue.expr3 <- tissue.expr2
tissue.expr3['sd'] <- apply(tissue.expr2, MARGIN = 1, FUN = sd)
tissue.expr3['mean'] <- apply(tissue.expr2, MARGIN = 1, FUN = mean)
tissue.expr3['max'] <- apply(tissue.expr2, MARGIN = 1, FUN = max)
tissue.expr3['min'] <- apply(tissue.expr2, MARGIN = 1, FUN = min)
tissue.expr3 <- transform(tissue.expr3, minrange2 = abs(min-mean)/2, maxrange2 = abs(max-mean)/2)

tissue.expr3 <- subset(tissue.expr3, sd<1)

  # 标准3：删除实际值减均值的绝对值不小于2的Gene
tissue.expr4 <- subset(tissue.expr3, minrange2<1 & maxrange2<1)
write.table(x = rownames(tissue.expr4), file = 'D:\\学习\\毕设\\数据\\data\\tissue重新筛选得到的HK ENSG.txt', quote = F, row.names = F)

# 从细胞系表达数据入手
cellline.expr1 <- read_csv('D:\\学习\\毕设\\数据\\data\\表达量\\other rna_celline.csv', col_names = T)
cellline.expr2 <- dcast(data = cellline.expr1, formula = Gene ~ Sample, value.var = 'Value')
  # 标准1：删除0表达的Gene
cellline.expr2[cellline.expr2 == 0] <- NA
cellline.expr2 <- cellline.expr2[complete.cases(cellline.expr2),]

  # 标准2：使ENSG为行名。取log2，保留标准差小于1的Gene。
row.names(cellline.expr2) <- cellline.expr2[,1]
cellline.expr2 <- cellline.expr2[,-1]

cellline.expr2 <- log2(cellline.expr2)
cellline.expr3 <- cellline.expr2
cellline.expr3['sd'] <- apply(cellline.expr2, MARGIN = 1, FUN = sd)
cellline.expr3['mean'] <- apply(cellline.expr2, MARGIN = 1, FUN = mean)
cellline.expr3['max'] <- apply(cellline.expr2, MARGIN = 1, FUN = max)
cellline.expr3['min'] <- apply(cellline.expr2, MARGIN = 1, FUN = min)
cellline.expr3 <- transform(cellline.expr3, minrange2 = abs(min-mean)/2, maxrange2 = abs(max-mean)/2)

cellline.expr3 <- subset(cellline.expr3, sd<1)

  # 标准3：删除实际值减均值的绝对值不小于2的Gene
cellline.expr4 <- subset(cellline.expr3, minrange2<1 & maxrange2<1)
write.table(x = rownames(cellline.expr4), file = 'D:\\学习\\毕设\\数据\\data\\cellline重新筛选得到的HK ENSG.txt', quote = F, row.names = F)
