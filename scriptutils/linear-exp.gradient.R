#!/usr/bin/env Rscript

library(methods)
source("/n/projects/apa/R/apa_tools.R")
N <- as.numeric(commandArgs(trailing=TRUE)[[1]])
asHex <- as.logical(commandArgs(trailing=TRUE)[[2]])

x=-2^-c(1:N)*2:(N+1)
#plot(1:N,x)
x=rescale(x,from=range(x),to=c(0,255))
x=paste0("#0000",as.hexmode(x))
#x
#viewPalette(x)


y=sapply(round(seq(0,200,length=N-3),0),function(i){paste(rep(as.character.hexmode(i,width=2),3),collapse="")})
y=c(rep("#000000",4), paste0("#",y[2:length(y)]))
#y
#viewPalette(y)


z=col2rgb(y)
z[3,]=col2rgb(x)[3,]
z=rgb2(z)
#z
#viewPalette(c(z,8))


w=col2rgb(z)
w[2,]=round(seq(0,255,length=N),0)
w=rgb2(w)
#w
#viewPalette(c(w,8))

if (asHex) {
    IM(w,handle=stdout())
} else {
    w2=t(col2rgb(w))
    for (i in 1:N) IM(i,w2[i,],handle=stdout())
}

