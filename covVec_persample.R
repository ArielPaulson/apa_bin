
### RUN THIS ON LINUX NOT WINDOWS -- coverage vectors take a lot of memory.

source("/n/projects/apa/R/apa_tools.R")  # for lineplot, row.norm, binify, zerofy, WriteXLS2, slice.list
options(stringsAsFactors=FALSE)

nbins <- 100        # number of bins to break coverage vectors into
drop.short <- TRUE  # drop coverage vectors with < nbins elements to begin with?  Otherwise, get stretched until length = nbins
summ.fun <- mean    # summarize coverages per bin using this function

alias <- commandArgs(trailing=TRUE)
path <- paste(alias,"/",alias,sep="")

## read coverage vector dataset(s) -- change names, add more output files as required

cov.vecs <- scan(paste(path,".covVecs.txt",sep=""), sep="\t", what=list("",0,"","", ""))
cov.vecs[[5]] <- lapply(cov.vecs[[5]], FUN=function(x){as.numeric(unlist(strsplit(x,',')))})   ## split up covg vec string into actual vector
for (i in 1:length(cov.vecs[[1]])) {
    if (length(cov.vecs[[5]][[i]])<cov.vecs[[2]][[i]]) {
       cov.vecs[[5]][[i]] <- c(cov.vecs[[5]][[i]], NA)  ## expecting only a loss of 1 empty bp at 3' end!!!
    }
}
bin.vecs <- lapply(cov.vecs[[5]], FUN=function(x){ binify(zerofy(x),nbins,summ.fun,FALSE,drop.short) })   ## collapse covg vecs into 100 bins each
nonzero <- falsify(sapply(cov.vecs[[5]],function(x){any(x>0)}))   # genes which actually have any coverage
nonzero[cov.vecs[[2]]<nbins] <- FALSE   # also drop any genes less than nbins bp
save(cov.vecs, bin.vecs, nonzero, file=paste(path,".covVecs.Rdata",sep=""))

## also collect other coverage stats

bin.mats <- do.call(rbind, bin.vecs)
rownames(bin.mats) <- cov.vecs[[1]]
bin.mats.z <- row.norm(bin.mats, divSD=T, na.rm=T)
bin.mats.pm <- t(apply(bin.mats, 1, function(x){x/max(x)}))
bin.mats.z[is.na(bin.mats.z)] <- 0
bin.mats.pm[is.na(bin.mats.pm)] <- 0
colnames(bin.mats) <- colnames(bin.mats.z) <- colnames(bin.mats.pm) <- c("Gene\t1",2:100)
write.table(bin.mats, file=paste(path,".covVecs.binMats.txt",sep=""), sep="\t", quote=F)
write.table(bin.mats.z, file=paste(path,".covVecs.binMatsZ.txt",sep=""), sep="\t", quote=F)
write.table(bin.mats.pm, file=paste(path,".covVecs.binMatsPM.txt",sep=""), sep="\t", quote=F)

cov.stats <- matrix( c(
   sapply(cov.vecs[[5]], sum, na.rm=T),     # total sequenced bp assigned to gene
   sapply(cov.vecs[[5]], median, na.rm=T),  # median coverage depth per gene
   sapply(cov.vecs[[5]], max, na.rm=T),  # median coverage depth per gene
   sapply(cov.vecs[[5]], sd, na.rm=T),   # coverage SD per gene
   sapply(cov.vecs[[5]], CV, na.rm=T),   # coverage CV per gene
   cov.vecs[[2]],
   sapply(cov.vecs[[5]], function(x){sum(!is.na(x))}),    # gene/exon space with nonzero coverage
   0,    # nonzero bp percent; fill in later
   sapply(bin.mats, max, na.rm=T),
   rowMedians(bin.mats, na.rm=T),
   rowSDs(bin.mats, na.rm=T),
   rowCVs(bin.mats, na.rm=T)
), length(cov.vecs[[1]]), 12, FALSE, list(cov.vecs[[1]],c("Gene\tSum","Median","Max","SD","CV","TotBp","NZBp","NZPct","BinMed","BinMax","BinSD","BinCV")) )
cov.stats[is.na(cov.stats)] <- 0
cov.stats[is.infinite(cov.stats[,3]),3] <- 0
cov.stats[,8] <- cov.stats[,7]/cov.stats[,6]
write.table(cov.stats, file=paste(path,".covVecs.covStats.txt",sep=""), sep="\t", quote=F)


## summarize by total coverage, avg coverage
totals <- matrix(c(
   colSums(bin.mats, na.rm=T, nonzero=T),
   colMeans(bin.mats, na.rm=T, nonzero=T),
   colMedians(bin.mats, na.rm=T, nonzero=T),
   colMeans(bin.mats.z, na.rm=T, nonzero=T),
   colMedians(bin.mats.z, na.rm=T, nonzero=T)
), 5, 100, TRUE, list(qw(Sum,Mean,Median,ZMean,ZMedian),c("Gene\t1",2:100)))
write.table(totals, file=paste(path,".covVecs.allmeans.txt",sep=""), sep="\t", quote=F)
