#!/usr/bin/env Rscript

mingap <- 100
source("/n/projects/apa/R/apa_tools.R")
input <- commandArgs(trailing=TRUE)
lastz <- read.delim(input[1], header=FALSE)

nr <- nrow(lastz)
stats1 <- lastz[nr-1,]
stats2 <- lastz[nr,]
lastz <- as.matrix(lastz[1:(nr-3),])
mode(lastz) <- "numeric"

breaks <- c(1, which(lastz[,7]<0 | lastz[,8] < 0), nrow(lastz)+1)
N <- length(breaks)-1
cols <- rainbow(N+1)[1:N]

ylines <- 2 * N + 3
seqwd <- max(unlist(stats1[1:4]))
imgwd <- 1500
imght <- 400 + N * 20
axtck <- signif(diff(seq(0,seqwd,length=10)),3)[1]
axseq <- seq(0,seqwd,axtck)

png(paste(input[2],"png",sep="."), imgwd, imght); par(cex=1.2)
null.plot(xlim=c(0,seqwd), ylim=c(0,ylines), main=paste("Lastz Alignment Blocks, gaps >=",mingap))
axis(1, at=axseq, labels=axseq, las=1)
rect(1,ylines-0.9, stats1[[2]],ylines, col=1, border=1)  # sequence 1
text(stats1[[2]]/2, ylines-0.45, col="white", font=2, labels=stats1[[7]])
rect(1,0, stats1[[4]],0.9, col=1, border=1)  # sequence 2
text(stats1[[4]]/2, 0.45, col="white", font=2, labels=stats1[[8]])
for (i in 1:N) {
  rect(lastz[breaks[i],1],ylines-i-0.9, lastz[breaks[i+1]-1,2],ylines-i, col=cols[i], border=cols[i])
  rect(lastz[breaks[i],3],N+2-i-0.9, lastz[breaks[i+1]-1,4],N+2-i, col=cols[i], border=cols[i])
  segments(lastz[breaks[i],1],ylines-i-0.9, lastz[breaks[i],3],N+2-i, lty=1, lwd=2, col=cols[i])
  segments(lastz[breaks[i+1]-1,2],ylines-i-0.9, lastz[breaks[i+1]-1,4],N+2-i, lty=1, lwd=2, col=cols[i])
  for (j in (breaks[i]+1):(breaks[i+1]-1)) {
    if (lastz[j,7]>=mingap) {
      rect(lastz[j-1,2],ylines-i-0.9, lastz[j,1],ylines-i, col="white", border="white")
      segments(lastz[j-1,2],ylines-i-0.45, lastz[j,1],ylines-i-0.45, lty=1, lwd=2, col=cols[i])
    }
    if (lastz[j,8]>=mingap) {
      rect(lastz[j-1,4],N+2-i, lastz[j,3],N+2-i-0.9, col="white", border="white")
      segments(lastz[j-1,4],N+2-i-0.45, lastz[j,3],N+2-i-0.45, lty=1, lwd=2, col=cols[i])
    }
  }
}
dev.off()
quit()
