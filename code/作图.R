#########################数据整理#########################
library('ggplot2')
library('reshape2')
library('corrplot')

# 读入CodonW处理CDS后生成的原始文件
data1 <- read.table(file = 'D:/学习/毕设/数据/data/全基因组/密码子偏好相关参数1.txt', header = T)
data2 <- read.table(file = 'D:/学习/毕设/数据/data/全基因组/密码子偏好相关参数2.txt', header = T)

# 合并数据集
names(data1)[names(data1) == 'title'] <- c('Gene_description')
data3 <- merge(x = data1, y = data2, by = 'Gene_description')

# 添加列GC12、CAI
data3 <- transform(data3, GC12 = (GC1 + GC2)/2)

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
  geom_hline(yintercept = 35) +
  scale_y_continuous(breaks = c(25,30,35,40,45,50,55,60,65)) +
  theme(axis.line = element_line(colour = 'black')) +
  xlab('GC3s')


#########################中性图#########################

ggplot(data3, aes(x = GC3s.x, y = GC12)) + 
  geom_point(shape = 16, size = 0.1) + 
  stat_smooth(method = lm) +
  theme(axis.line = element_line(colour = 'black')) +
  annotate('text', label = 'r=0.58, p<2.2e-16', x=0.4, y=0.7) +
  xlab('GC3s')

cor.test(x = data3$GC12, y = data3$GC3s.x, method = 'spearman', exact = F)

#########################对应分析#########################

# 每轴贡献率
inertia <- read.table(file = 'D:\\学习\\毕设\\数据\\data\\对应分析\\Correspondence analysis results from [CodonW] on data 2\\前40轴相对和累积贡献率.txt', header = T)
ggplot(inertia, aes(x = factor(Num.), y = R.Iner.)) +
  theme(axis.line = element_line(colour = 'black')) +
  xlab('轴') + ylab('相对贡献率') +
  geom_bar(stat = 'identity')

# 密码子第三位碱基的一二轴分布
codon.coa <- read.csv(file = 'D:\\学习\\毕设\\数据\\data\\对应分析\\每密码子前四轴分布2.csv', header = T)
ggplot(codon.coa, aes(x = Axis1, y = Axis2, shape = factor(type))) +
  theme(axis.line = element_line(colour = 'black')) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(size = 1) + 
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


#########################最优密码子#########################

# 根据ENC值排序data3，取出首尾10%基因名
data3 <- data3[order(data3$Nc),]
tail.enc <- as.character(head(x = data3$Gene_description, n = nrow(data3) * 0.1))
head.enc <- as.character(tail(x = data3$Gene_description, n = nrow(data3) * 0.1))

write.table(x = head.enc, file = 'D:/学习/毕设/数据/data/极端ENC/ENC极端高 基因.txt', row.names = F, quote = F)
write.table(x = tail.enc, file = 'D:/学习/毕设/数据/data/极端ENC/ENC极端低 基因.txt', row.names = F, quote = F)


#########################获取低ENC高GC基因#########################

Lowenc.Lowgc <- subset(x = data3, subset = GC3s.x<0.3 & Nc<49, select = c(Gene_description, GC3s.x, Nc))
write.table(x = Lowenc.Lowgc$Gene_description, file = 'D:/学习/毕设/数据/data/HELA/GC3s.x0.3 Nc49 acc.txt', quote = F, row.names = F)
Lowenc.Highgc <- subset(x = data3, subset = GC3s.x>0.8 & Nc<36, select = c(Gene_description, GC3s.x, Nc))
write.table(x = Lowenc.Highgc$Gene_description, file = 'D:/学习/毕设/数据/data/HELA/GC3s.x0.8 Nc36 acc.txt', quote = F, row.names = F)

