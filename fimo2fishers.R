source("~/uapa/R/apa_tools.R")

load("/n/data1/biobase/transfac/current/R/matrix.dat.RData")
vert <- indexes$acc[grep("^V",indexes$ID)]
urls <- rbind(read.delim("transfac.url", sep="\t", header=F, as.is=T),read.delim("jaspar.url", sep="\t", header=F, as.is=T))
ids <- rbind(read.delim("transfac.ids", sep="\t", header=F, as.is=T)[,2:3],read.delim("jaspar.ids", sep="\t", header=F, as.is=T)[,2:3])

### REGULAR MOTIF MEME, TSS-1KB

ff <- read.delim("fimo2fishers.f2f.txt", sep="\t", header=T, as.is=T)

ff.js <- ff[ff[,3] == "jaspar",]
ff.tf <- ff[ff[,3] == "transfac",]
ff.me <- ff[grep("meme",ff[,3]),]
ff.tf <- ff.tf[ff.tf[,4] %in% vert,]  # subset on vertebrate motifs

ff2.tf <- cbind(ff.tf[,c(1:5,10)],fisherize(ff.tf[,6:9]))
ff2.js <- cbind(ff.js[,c(1:5,10)],fisherize(ff.js[,6:9]))
ff2.me <- cbind(ff.me[,c(1:5,10)],fisherize(ff.me[,6:9]))

ff3 <- rbind(ff2.tf, ff2.js, ff2.me)[,c(1:5,7:16,6)]
udb <- unique(ff3[,3]); udb
for (x in udb) { ff3[ff3[,3]==x,15] <- p.adjust(ff3[ff3[,3]==x,14], method="BH") }

min(ff3[,15])
min(ff3[,14])
sum(ff3[,15] <= 0.05)
sum(ff3[,14] <= 0.05)
sum(ff3[,14] <= 0.01)

sig <- which(ff3[,14] <= 0.05)
length(sig)
luniq(ff3[sig,4])
sort(unique(ff3[sig,4]))

ff3[sig,]
WriteXLS2(list(ff3[sig,]), "sig.f2f.motifs.xls")



### REGULAR MOTIF MEME, TSS+-1KB

ff <- read.delim("fimo2fishers.TSS.txt", sep="\t", header=T, as.is=T)
ff[grep("1k$",ff[,3]),3] <- paste(ff[grep("1k$",ff[,3]),3],"meme",sep=".")

ff.js <- ff[ff[,3] == "jaspar",]
ff.tf <- ff[ff[,3] == "transfac",]
ff.me <- ff[grep("meme",ff[,3]),]
ff.tf <- ff.tf[ff.tf[,4] %in% vert,]  # subset on vertebrate motifs

ff2.tf <- cbind(ff.tf[,c(1:5,10)],fisherize(ff.tf[,6:9]))
ff2.js <- cbind(ff.js[,c(1:5,10)],fisherize(ff.js[,6:9]))
ff2.me <- cbind(ff.me[,c(1:5,10)],fisherize(ff.me[,6:9]))

ff3 <- rbind(ff2.tf, ff2.js, ff2.me)[,c(1:5,7:16,6)]
udb <- unique(ff3[,3]); udb
for (x in udb) { ff3[ff3[,3]==x,15] <- p.adjust(ff3[ff3[,3]==x,14], method="BH") }

min(ff3[,15])
min(ff3[,14])
sum(ff3[,15] <= 0.05)
sum(ff3[,14] <= 0.05)
sum(ff3[,14] <= 0.01)

sig <- which(ff3[,14] <= 0.05)
length(sig)
luniq(ff3[sig,4])
sort(unique(ff3[sig,4]))

ff3[sig,]
WriteXLS2(list(ff3[sig,]), "sig.TSS.motifs.xls")



### UTR MIRNA-SITE MEME

ff <- read.delim("fimo2fishers.UTR.txt", sep="\t", header=T, as.is=T)
ff3 <- cbind(ff[,c(1:5,10)],fisherize(ff[,6:9]))[,c(1:5,7:16,6)]

udb <- unique(ff3[,3]); udb
for (x in udb) { ff3[ff3[,3]==x,15] <- p.adjust(ff3[ff3[,3]==x,14], method="BH") }

min(ff3[,15])
min(ff3[,14])
sum(ff3[,15] <= 0.05)
sum(ff3[,14] <= 0.05)
sum(ff3[,14] <= 0.01)

sig <- which(ff3[,14] <= 0.05)
length(sig)
luniq(ff3[sig,4])
sort(unique(ff3[sig,4]))

ff3[sig,]
WriteXLS2(list(ff3[sig,]), "sig.UTR.motifs.xls")
