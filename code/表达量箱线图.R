data.gc0.3 <- read.csv('D:\\学习\\毕设\\数据\\data\\HELA\\GC0.3 TPM SRR3589956.csv', header = T)
data.gc0.8 <- read.csv('D:\\学习\\毕设\\数据\\data\\HELA\\GC0.8 TPM SRR3589956.csv', header = T)
data1 <- na.omit(c(data.gc0.3$TPM1, data.gc0.8$TPM1))

len1 <- sum(complete.cases(data.gc0.3$TPM1))
len2 <- sum(complete.cases(data.gc0.8$TPM1))
f <- factor(rep(c("低GC","高GC"), c(len1, len2)))

boxplot(data1~f)
boxplot(data1~f, outline = F)
