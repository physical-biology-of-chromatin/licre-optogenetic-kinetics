time = read.table("licre-fit-example-time.txt")
time = as.numeric(time[,1])
data = read.table("licre-fit-example-data.txt")
solu = read.table("licre-fit-example-solu.txt")

ymax = max(c(as.matrix(data),
	     as.matrix(solu)))
ylim = c(0,ymax)
pdf(file = "licre-fit-example-half.pdf")
# init plot
plot(time,
     data[,7],
     type = "n",
     xlab = "time",
     ylab = "RU",
     ylim = ylim)

# data as points and solu (fits) as lines
for (i in 1:8){
	points(time, data[,i], col = i)
	lines( time, solu[,i], col = i)
}
dev.off()


pdf(file = "licre-fit-example-full.pdf")
# init plot
plot(time,
     data[,9],
     type = "n",
     xlab = "time",
     ylab = "RU",
     ylim = ylim,
     main = "LiCre on Full LoxP")

# data as points and solu (fits) as lines
for (i in 10:15){
	points(time, data[,i], col = i, cex = 0.7, pch = "x")
	lines( time, solu[,i], col = i, lwd = 3)
}
dev.off()


