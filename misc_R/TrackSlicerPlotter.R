#setwd("/home/apa/TrackSlicerPlots")

data <- read.delim("TrackSlicerRData1.txt", sep="\t", fill=T, skip=1, header=F)
meta <- read.delim("TrackSlicerRData2.txt", sep="\t", fill=T, skip=1, header=F)	# gene chr start end strand, ablines, {exon adj.start adj.end} foreach exon

Ldata <- log2(data)
for(i in 2:ncol(Ldata)) { Ldata[is.infinite(Ldata[,i]),i] <- 0 }
Ldata[,1] <- data[,1]

gene <- meta[1,1]
boundaries <- c(meta[3,2],meta[nrow(meta),3])
title <- paste(meta[1,1],", chr",meta[1,2],":",meta[1,3],"-",meta[1,4],"(",meta[1,5],") ± 1Kb",sep="")
Ltitle <- paste(title,", Log2",sep="")
x <- c(1:nrow(data))
maxd <- max(data[,2:ncol(data)])
recmin <- -0.1*maxd
recmax <- 1.1*maxd
cols <- c(2:7)

png(paste(gene, "_Linear.png", sep=""), height=600, width=1200)
plot(x, data[,2], type="l", col=0, ylim=c(recmin,recmax), ylab="depth", xlab="bp", main=title)
for(j in 3:nrow(meta)) { rect(meta[j,2],recmin,meta[j,3],recmax,col="grey85",density=NULL,border=NA) }
abline(h=meta[2,2:3], col=1, lty=2)
abline(v=boundaries, col=1, lwd=2)
for(i in 2:ncol(data)) { lines(x, data[,i], col=cols[i-1]) }
dev.off()

png(paste(gene, "_Log2.png", sep=""), height=600, width=1200)
plot(x, Ldata[,2], type="l", col=0, ylim=c(0,log2(recmax)), ylab="Log2(depth)", xlab="bp", main=Ltitle)
for(j in 3:nrow(meta)) { rect(meta[j,2],-0.1,meta[j,3],log2(recmax), col="grey85", density=NULL, border=NA) }
abline(h=log2(meta[2,2:3]), col=1, lty=2)
abline(v=boundaries, col=1, lwd=2)
for(i in 2:ncol(Ldata)) { lines(x, Ldata[,i], col=cols[i-1]) }
dev.off()


