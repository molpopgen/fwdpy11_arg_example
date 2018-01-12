x=read.table("nmuts.txt")

png("qq.png")
qqplot(x$V1,x$V2,xlab="Us",ylab="Jerome")
abline(0,1)
dev.off()

png("xy.png")
plot(x$V1,x$V2,xlab="Us",ylab="Jerome")
abline(0,1)
dev.off()
