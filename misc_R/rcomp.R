#!/usr/bin/env Rscript

source("/home/apa/apa_tools.R")

live <- read.delim("live.txt", as.is=TRUE, header=FALSE)
dev <- read.delim("dev.txt", as.is=TRUE, header=FALSE)

dat <- list(
    num=data.frame(package=suniq(c(live[,1],dev[,1])), live=0, dev=0, diff=0, lfc=0),
    size=data.frame(package=suniq(c(live[,1],dev[,1])), live=0, dev=0, diff=0, lfc=0)
)
for (i in 1:2) {
    dat[[i]]$live[match(live[,1],dat[[i]][,1])] <- live[,i+1]
    dat[[i]]$dev[match(dev[,1],dat[[i]][,1])] <- dev[,i+1]
    dat[[i]]$diff <- dat[[i]]$dev-dat[[i]]$live
    dat[[i]]$lfc <- log2((dat[[i]]$dev+0.5)/(dat[[i]]$live+0.5))
    dat[[i]] <- dat[[i]][!grepl("^00",dat[[i]][,1]),]
}
nrow(dat[[1]])
sum(dat[[1]]$live==0)
sum(dat[[1]]$dev==0)

pa <- dat[[1]]$live==0 | dat[[1]]$dev==0  # pres/abs
sum(pa)

par(mfrow=c(2,2))
for (i in 1:2) {
    dhist(dat[[i]]$diff, points=TRUE)
    dhist(dat[[i]]$lfc, points=TRUE)
}

## N files
sum(!pa&dat[[1]]$diff!=0); sum(!pa&dat[[1]]$diff!=0)/nrow(dat[[1]])
dat[[1]][!pa&dat[[1]]$diff!=0,]

## total size
sum(!pa&dat[[2]]$diff<0); sum(!pa&dat[[2]]$diff<0)/nrow(dat[[2]])
dat[[2]][!pa&abs(dat[[2]]$diff)>1E5,]



liveA <- read.delim("live_A.txt", sep=" ", as.is=TRUE, header=FALSE)[,c(9,5)]
devA <- read.delim("dev_A.txt", sep=" ", as.is=TRUE, header=FALSE)[,c(9,5)]
liveA[,1] <- sub(".*library/","",liveA[,1])
devA[,1] <- sub(".*library/","",devA[,1])
uA <- sort(union(liveA[,1],devA[,1]))
compA <- data.frame(package=sub("/.*","",uA),file=sub("^.*?/","",uA),live=liveA[match(uA,liveA[,1]),2],dev=devA[match(uA,devA[,1]),2])
for (i in 3:4) compA[[i]][is.na(compA[[i]])] <- 0
compA <- cbind(compA, diff=compA$dev-compA$live, lfc=log2((compA$dev+0.5)/(compA$live+0.5))); head(compA)
as.matrix(sapply(split(compA$diff,compA$package),sum))


par(mfrow=c(2,2))
dhist(compA$diff, points=TRUE)
dhist(compA$lfc, points=TRUE)
dhist(compA$diff, points=TRUE, xlim=c(-1E5,1E5))
dhist(compA$lfc, points=TRUE, xlim=c(-2,2))

x2 <- cbind(
    dat$size[grep("^[Aa]",dat$size[,1]),],
    liveA=round(as.matrix(sapply(split(compA$live,compA$package),sum)/1024),0),
    devA=round(as.matrix(sapply(split(compA$dev,compA$package),sum)/1024),0),
    diffA=round(as.matrix(sapply(split(compA$diff,compA$package),sum)/1024,2))
)
n2 <- nrow(x2)
sum(x2$diff<0);sum(x2$diff<0)/n2
sum(x2$diffA<0);sum(x2$diffA<0)/n2

cAd <- c(LT=sum(compA$diff<0), EQ=sum(compA$diff==0), GT=sum(compA$diff>0))
cAd; cAd/nrow(compA)






