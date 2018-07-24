
## http://r.789695.n4.nabble.com/Drawing-lines-in-margins-td897397.html: see grconvertX, grconvertY, par(xpd=TRUE)

chars <- sapply(32:126, function(i) rawToChar(as.raw(i)))
C <- length(chars)

char.strL <- lapply(chars, function(r) sapply(seq(1,30,1), function(i) paste(c(rep(r,i),".",i),collapse="") ))
char.strR <- lapply(chars, function(r) sapply(seq(1,30,1), function(i) paste(c(i,".",rep(r,i)),collapse="") ))
char.strL20 <- sapply(1:C, function(i) paste(c(rep(chars[i],20),".",c(32:126)[i]),collapse="") )
char.strR20 <- sapply(1:C, function(i) paste(c(c(32:126)[i],".",rep(chars[i],20)),collapse="") )

cex <- 1.2
imgw <- 700

dir <- paste0("mar_test.images/")
system(paste("mkdir -p",dir))

png(paste0(dir,"master.png"), 700, 1500)
par(mar=c(0,20,0,0), cex=cex, las=2, xaxs="i", yaxs="i")
#par(mar=c(0,18.5,0,0), cex=1.2, las=2, xaxs="i", yaxs="i")
null.plot(xlim=c(0,20),ylim=c(0,96))
text(-0.4, 95:1, char.strL20, pos=4)
axis(2, 95:1, char.strR20)
abline(v=seq(0,150,1), lty=c(5,rep(1:5,10)), col="#00000055")
abline(v=seq(0,150,1)+0.5, lty=3, col="#0000FF44")
dev.off()



for (r in 1:95) {   # length chars
    IM(r)
    mat <- matrix(0,30,30)
    r2 <- c(32:126)[r]
    k <- 0
    for (j in 1:4) {
        k1 <- ifelse(k+1>9, k+1, paste0("0",k+1))
        k4 <- ifelse(k+4>9, k+4, paste0("0",k+4))
        imgw2 <- ifelse(j==4,imgw+100,imgw)
        png(paste0(dir,"char_",r2,".k_",k1,"-",k4,".cex_",sub("\\.","-",cex),".png"), imgw2, imgw2)
        par(mar=c(k+1,k+2,k+3,k+4), cex=cex, las=2)
        null.plot(xlim=c(1,30), ylim=c(1,30))
        text(c(15,0,15,30), c(0,15,30,15), c(k+1,k+2,k+3,k+4))
        for (a in 1:2) axis(a, 1:30, char.strR[[r]])
        for (a in 3:4) axis(a, 1:30, char.strL[[r]])
        dev.off()
        k <- k+4
    }
}

## ANNOTATE RESULTS IN mar.test.xlsx
