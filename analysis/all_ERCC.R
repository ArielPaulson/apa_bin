
source("U:/apa/R/apa_tools.R")
options(stringsAsFactors=FALSE)

tar <- read.delim("../targets.txt", as.is=T); dim(tar)
#tar <- tar[c(1:360,393:410,361:392),]  # swap UHR_polyA with UHR_Ribo
for (i in 3:4) tar[,i] <- sub("UHR_","UHR",sub("UC_","UC",sub("U_T","UT",sub("Fred_","Fred",sub("Northwestern","Northwest",sub("Ribo","ribo",tar[,i]))))))

tar2 <- read.delim("../supertargets.txt")
tar2[,3:4] <- tar[,3:4]
libN2 <- sub("UHR","",sub("poly","Poly",sub("western","W",gsub("_","",tar2[,4]))))
newN2 <- sub("UHR_","",sub("poly","Poly",tar2[,23]))
ordkits <- qw(Clontech,Miltenyi,NuGen,Sigma,PolyA,Ribo)
ordkits2 <- qw(P1C,P2C,P3C,P1M,P2M,P2N,P3N,P1S,P2S,P3S,PolyA,Ribo)
x <- mat.split(tar2,tar2[,23])
tar.ko <- do.call(rbind, x[mgrep(ordkits2,names(x),T)]); unique(tar.ko[,c(3,23)]); rm(x)


aln <- read.delim("../datasets/ALL_RDS/all_bamreads.txt")[,c(1,6)]
aln <- aln[match(tar[,2],aln[,1]),2]

dat <- as.matrix(read.delim("../datasets/all_ERCC.txt", row.names=1))
len <- dat[,1]
dat <- dat[,!grepl("ERCC",colnames(dat))]
dat <- dat[,match(tar[,2],colnames(dat))] # also gets rid of col 1
colnames(dat) <- tar[,2]
dim(dat)

sum(is.na(dat[68,]))
#lost <- which(is.na(dat[68,]))
sum(is.na(dat))
dat[is.na(dat)] <- 0
sum(rowSums(dat)==0)

rpkm <- 1E9 * t( t(dat/len) / aln )
srpkm <- rpkm[order(rowMedians(rpkm, nonzero=T)),]
lrpkm <- NAfy(log2(srpkm))
lrpkm2 <- NAfy(log2(srpkm))

png("ERCC_by_number.png", 1200, 600); par(mar=c(8,4,4,2), cex=1.2, las=2)
boxplot(t(lrpkm))
dev.off()

ccol <- c(2:6)[match(tar[,5],unique(tar[,5]))]
ccol[tar[,6]>9] <- hsv.tweak(ccol[tar[,6]>9], s=-0.7)
axseq <- sapply(find.runs(as.factor(tar[,3])),mean)
abstart <- do.call(cbind,find.runs(as.factor(tar[,3])))[1,2:length(axseq)]-0.5
png("ERCC_by_sample.png", 3200, 600)
par(mar=c(4,4,4,2), cex=1.2, cex.axis=0.85, las=1, xaxs="i")
boxplot(lrpkm, border=ccol, xaxt="n")
abline(v=abstart, col=8, lty=2)
axis(1, axseq, tick=F, labels=gsub("_","",apply(tar[axseq,4:5], 1, paste, collapse="\n")))
dev.off()

png("ERCC_by_sample.dot.png", 3500, 600); par(mar=c(4,4,4,2), cex=1.2, las=1, xaxs="i")
dotplot(lrpkm, col=palettizer(c(6,4,5,3,7,2),92), xaxt="n", legend=NA, xlab="", ylab="Log2 RPKM", main="92 ERCC controls per replicate, colored by rank of global median")
abline(v=abstart, col=8, lty=2)
axis(1, axseq, tick=F, labels=gsub("_","",apply(tar[axseq,4:5], 1, paste, collapse="\n")))
dev.off()

png("ERCC_by_number.dot.png", 1200, 1000); par(mar=c(8,4,4,2), cex=1.2, las=2)
dotplot(t(lrpkm), col=ccol, pch=rep(1:12,times=c(rep(27,10),18,27)), legend=NA, ylim=c(-5,20), xlab="")
legend("topleft", bty="n", col=ccol[axseq], pch=rep(1:12,times=c(rep(3,10),2,3)), legend=apply(tar[axseq,4:5],1,paste,collapse="_"))
dev.off()

scol <- rep(rainbow(13)[1:12], times=listLengths(find.runs(as.factor(tar[,4]),term=F)))
png("ERCC_by_number.line.png", 1200, 600); par(mar=c(8,4,4,2), cex=1.2, las=2, xaxs="i")
lineplot(t(lrpkm), col=scol, legend=NA, xlab="")
legend("topleft", bty="n", lty=1, col=rainbow(13)[1:12], legend=unique(tar[,4]))
dev.off()

idim <- c(3000,3000); pmar <- c(15,15)
count.cor <- corr.mat(log2(dat), reorder=T, pmar=pmar, imgdim=idim, filename="count.corr.png")
rpm.cor <- corr.mat(log2(t(t(dat)/aln)), reorder=T, pmar=pmar, imgdim=idim, filename="rpm.corr.png")
rpkm.cor <- corr.mat(lrpkm, reorder=F, pmar=pmar, imgdim=idim, filename="rpkm.corr.png")
rpkm.cor.ord <- corr.mat(lrpkm, reorder=T, pmar=pmar, imgdim=idim, filename="rpkm.corr.ord.png")

idim2 <- c(1000,1000); pmar2 <- c(8,8)
rpkm.avg <- aggregate.cols(rpkm, list(tar[,3]), mean)
lrpkm.avg <- aggregate.cols(lrpkm, list(tar[,3]), mean)
rpkm.avg.cor.ord <- corr.mat(lrpkm.avg, reorder=T, pmar=pmar2, imgdim=idim2, cex=1.2, main="ERCC92, Pearson Correlations between Averaged Replicate Log2 RPKMs, Clustered", filename="rpkm.avg.corr.ord.png")
rpkm.avg.cor <- corr.mat(lrpkm.avg, reorder=F, pmar=pmar2, imgdim=idim2, cex=1.2, main="ERCC92, Pearson Correlations between Averaged Replicate Log2 RPKMs", filename="rpkm.avg.corr.png")

seq1 <- seq(1,28,3); 
libs <- unique(tar.ko[,3])
libs <- c(matrix(libs[c(seq1,seq1+1,seq1+2)], ncol=10, byrow=T), libs[c(33:31,35:34)])
ko <- match(libs,colnames(lrpkm.avg))
corr.mat(lrpkm.avg[,ko], names1=unique(newN2)[ko], namse2=unique(newN2)[ko], pmar=pmar2, imgdim=idim2, cex=1.2, main="ERCC92, Pearson Correlations between Averaged Replicate Log2 RPKMs", filename="rpkm.avg.corr_ANON.png")
ko2 <- ko[c(1:30,34:35,31:33)]
corr.mat(lrpkm.avg[,ko2], names1=unique(newN2)[ko2], namse2=unique(newN2)[ko2], pmar=pmar2, imgdim=idim2, cex=1.2, main="ERCC92, Pearson Correlations between Averaged Replicate Log2 RPKMs", filename="rpkm.avg.corr_ANON_reord.png")


dhist(list(as.dist(rpkm.avg.cor)))
abline(v=median(as.dist(rpkm.avg.cor)), col=2)
abline(v=c(rpkm.avg.cor[31:33,34:35]), col=4)
q <- find.quantile(c(rpkm.avg.cor[31:33,34:35]), as.dist(rpkm.avg.cor))

write.table(rpkm.avg.cor, "rpkm.avg.cor.txt", sep="\t", quote=F)

linear <- apply(lrpkm, 2, cor, y=1:92, use="complete.obs")
plot(sort(linear))

ccol2 <- c(2:6)[match(tar[!duplicated(tar[,3]),5],unique(tar[,5]))]
linear2 <- split(linear,tar[,3]); linear2 <- linear2[match(unique(tar[,3]),names(linear2))]
png("ERCC_log2_linearity.png", 600, 800); par(mar=c(4,8,4,2), cex=1.2, las=1, yaxs="i")
boxplot(rev(linear2), border=rev(ccol2), horizontal=T, main="Log2-Linearity of ERCC92 RPKM Values", xlab="Pearson R Value")
abline(h=seq(2,32,3)+0.5, col=8, lty=2)
dev.off()

png("ERCC_log2_linearity_ANON.png", 600, 800); par(mar=c(4,8,4,2), cex=1.2, las=1, yaxs="i")
boxplot(rev(linear2[ko]), names=rev(unique(newN2)[ko]), border=rev(ccol2[ko]), horizontal=T, main="Log2-Linearity of ERCC92 RPKM Values", xlab="Pearson R Value")
abline(h=seq(2,32,3)+0.5, col=8, lty=2)
dev.off()

new.ord <- ko[c(22:30,1:9,10:15,16:21,31:35)]
png("ERCC_log2_linearity_ANON_wide.png", 800, 600); par(mar=c(8,4,4,2), cex=1.2, las=2, xaxs="i")
boxplot(linear2[new.ord], names=unique(newN2)[new.ord], border=ccol2[new.ord], horizontal=F, main="Log2-Linearity of ERCC92 RPKM Values", ylab="Pearson R Value")
abline(v=seq(2,32,3)+1.5, col=8, lty=2)
dev.off()

final.new.ord <- new.ord[c(10:30,1:9,34:35,31:33)]
ablines2 <- c(seq(2,30,3)+1.5, 32.5)
png("ERCC_log2_linearity_ANON_wide_reord.png", 800, 600); par(mar=c(8,4,4,2), cex=1.2, las=2, xaxs="i")
boxplot(linear2[final.new.ord], names=unique(newN2)[final.new.ord], border=ccol2[final.new.ord], horizontal=F, main="Log2-Linearity of ERCC92 RPKM Values", ylab="Pearson R Value")
abline(v=ablines2, col=8, lty=2)
dev.off()


table(tar[dat[68,]>0,4:5])

aov.rpkm <- rpkm[,!grepl("OU",colnames(rpkm))]
aov.src <- as.factor(tar[match(colnames(aov.rpkm),tar[,2]),4])
aov.conc <- as.factor(tar[match(colnames(aov.rpkm),tar[,2]),5])
aov.rep <- as.factor(tar[match(colnames(aov.rpkm),tar[,2]),6])
ercc <- as.factor(rownames(dat))

pre.aov <- data.frame(ercc=rep(ercc,ncol(aov.rpkm)), rep=rep(aov.rep,92), src=rep(aov.src,92), conc=rep(aov.conc,92), rpkm=c(aov.rpkm))

anova.ercc <- function(x) {
	y <- pre.aov[pre.aov[,1]==x,3:5]
	z <- anova(lm( rpkm ~ src + conc + src*conc, data=y, na.action=na.omit ))
	z[[5]][1:3]
}
anova.src <- function(x) {
	y <- pre.aov[pre.aov[,3]==x,c(1,4,5)]
	z <- anova(lm( rpkm ~ ercc + conc + ercc*conc, data=y, na.action=na.omit ))
	z[[5]][1:3]
}
anova.rep <- function(x) {
	y <- pre.aov[pre.aov[,3]==x,c(2,4,5)]
	z <- anova(lm( rpkm ~ rep + conc + rep*conc, data=y, na.action=na.omit ))
	z[[5]][1:3]
}

aov.raw <- t(sapply(ercc,anova.ercc))
aov.adj <- matrix(p.adjust(aov.raw,method="BH"), ncol=3)
min(aov.raw); min(aov.adj)
sig <- which(apply(aov.adj<=0.05,1,any)); length(sig)
aov.adj[sig,]
aov.raw[sig,]
as.character(ercc[sig])

aov.raw <- t(sapply(unique(aov.src),anova.src))
aov.adj <- matrix(p.adjust(aov.raw,method="BH"), ncol=3)
min(aov.raw); min(aov.adj)
sig <- which(apply(aov.adj<=0.05,1,any)); length(sig)
aov.adj[sig,]
aov.raw[sig,]
as.character(unique(aov.src)[sig])

aov.raw <- t(sapply(unique(aov.src),anova.rep))
aov.adj <- matrix(p.adjust(aov.raw,method="BH"), ncol=3)
min(aov.raw); min(aov.adj)
sig <- which(apply(aov.adj<=0.05,1,any)); length(sig)
aov.adj[sig,]
aov.raw[sig,]
as.character(unique(aov.src)[sig]); as.character(unique(aov.src)[-sig])



save.image("all_ERCC.RData")


