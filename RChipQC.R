#!/usr/bin/env Rscript

## Calculates various QC metrics for a ChIP bam file, calls all possible peaks, bins genome signal, cross-correlates, etc, etc...
## command line: /home/apa/local/bin/RChipQC.R <arguments, in order>
## IF VARIABLE READ LENGTH: use mean or median length

## Table of key objects in the session RData that you may want to access:
## covg            : coverage(bam)
## covg.hist       : list with 3 elements: coverage percentiles, full coverage histogram, binned coverage histogram
## CCF             : coverage cross-correlation data, using shift sizes 0-1000
## top.CCF         : same as CCF, but using only the regions of the top 1000 peaks
## islands.by.chr  : GRangesList of all islands (covg>0), with some stats
## totalCoverage   : table with rows=depths and cols=chroms, cells=N bp on chrom with depth D (total=genome size)
## bpAtDepth       : totalCoverage multiplied by the depth: cells=sequenced bases aligned to chrom assigned to positions with depth D (total=sum of all read lengths))
## bin.sums.1000   : sum signal across genome in 1000bp bins (GRanges object, to allow subsetting on regions of interest)
## 
## Objects which may or may not exist, depending on arguments:
## bam             : bam data as a GRanges object
## bin.sums.100    : sum signal across genome in   100bp bins (GRanges object, to allow subsetting on regions of interest)
## bin.sums.10000  : sum signal across genome in 10000bp bins (GRanges object, to allow subsetting on regions of interest)

### Setup
ca <- commandArgs(TRUE)
bamfile <- ca[1]                                             # input BAM file
geno <- ca[2]                                                # provide a genome label for extra chromosome data, or "NA"
minht <- ifelse(length(ca)>2, as.numeric(ca[3]), 1)          # minimum height of peak to report (default 1)
fullpath <- ifelse(length(ca)>3, as.logical(ca[4]), FALSE)   # use full path for prefix (only affects symlinks) (default FALSE)
keepbam  <- ifelse(length(ca)>4, as.logical(ca[5]), FALSE)   # keep 'bam' in RData object
matepair <- ifelse(length(ca)>5, as.logical(ca[6]), FALSE)   # analyze mate-pair behavior (takes a lot more time / RAM)

if (geno=="NA") geno <- NA
if (is.na(minht)) minht <- 10
if (is.na(fullpath)) fullpath <- FALSE
if (is.na(keepbam)) keepbam <- FALSE
if (!file.exists(bamfile)) stop(paste0("BAM file '",bamfile,"' not found!  Halting.\n"))

message("ARGUMENTS:\nBamfile: ",bamfile,"\nGeno: ",geno,"\nMinht: ",minht,"\nFullpath: ",fullpath,"\nKeepbam: ",keepbam,"\nMatepair: ",matepair,"\n")

## http://www.bioconductor.org/help/course-materials/2010/SeattleJan10/day3/ChIPSeq.pdf -- Original Code Concept


### Prefix, temp, etc.
bamfile2 <- ifelse(fullpath, system(paste("readlink -f",bamfile),intern=TRUE), canonicalize(bamfile))
prefix <- sub(".bam$",".RChipQC",bamfile2)
message("Using prefix: ",prefix)
rdata <- paste0(prefix,".RData")
rdata.mp <- paste0(prefix,".matePair.RData")  # if doing mate-pair analysis
rdata.tmp <- paste0(prefix,".tmp.RData")
tmp <- tempdir()
have.cs <- FALSE    # initially
have.mito <- FALSE  # initially


message("\nLoading libraries...")

library(methods)
suppressMessages(require(GenomicRanges))
suppressMessages(require(rtracklayer))
suppressMessages(require(chipseq))
suppressMessages(require(csaw))
source("~/lbin/apa_tools.R")


message("\nAnalyzing BAM...")

### Read BAM / get some basic numbers
bam <- as(import(bamfile2, "bam"),"GRanges")
bam.wd <- as.numeric(width(bam))   # alignment length on genome, which is mostly identical to read length, but will vary due to splicing, indels, etc.  // numeric to avoid integer overflows
readlen <- median(bam.wd)  # taking median read length in case of trimming
aligns <- length(bam.wd)   # alignments in bam
bpAln <- sum(bam.wd)       # total bp sequenced (sum of all read lengths)
upos <- sum(sapply(split(start(bam),as.numeric(seqnames(bam))),luniq))  # unique positions where alignments occur
message("alignments: ",aligns,"\ntot bp align: ",bpAln,"\nunq pos pct: ", round(100*upos/aligns,2),"\nmedian read length: ",readlen)

### Get chrom lengths, and if possible exclude mitochondria
chrsizes <- seqlengths(bam)
chrsizes <- nameless(data.frame(names(chrsizes),c(chrsizes)))
if (!is.na(geno)) {
    chrdata <- paste0("/n/data1/genomes/indexes/",geno,"/extras/",geno,".chrom_data.txt")
    if (!file.exists(chrdata)) stop(paste0("Genome argument '",geno,"' is not a local genome build!  Halting.\n"))
    chrdata <- read.delim(chrdata, as.is=TRUE)
    for (i in 1:ncol(chrdata)) chrdata[falsify(chrdata[,i]==""),i] <- NA
    chrdata <- chrdata[!is.na(chrdata[,1]),]  # only chroms in the actual reference fasta (fasta rank should always be in column 1)
    allchrs <- chrdata[,2]                    # and reference names should always be in column 2
    if (!all(chrsizes[,1] %in% allchrs)) stop(paste0("Some of the chromosomes in bam '",bamfile,"' are not found in genome '",geno,"'!\n"))
    have.cs <- TRUE
    somchrs <- chrdata[chrdata$Assembly=="chromosome" & !chrdata$Circular,1]  # should exclude mito
    mitochr <- chrdata[chrdata$Content=="mitochondrion",1]
    if (length(somchrs)==0) somchrs <- chrsizes[!chrdata$Circular,1]  # if nothing assembled to chr level, then take everything (which isn't plasmid/mito)
} else {
    allchrs <- chrsizes[,1]   # don't know for sure; take all
    is.mito <- grepl("(^MT$|chrM$|chrMT$|mito)",allchrs,TRUE) & chrsizes < 20000  # the usual suspects
    somchrs <- allchrs[!grepl("_",allchrs) & !is.mito]
    mitochr <- allchrs[is.mito]
    if (length(somchrs)==0) somchrs <- setdiff(allchrs, mitochr)  # if none passing filter, take everything non-mito (apparently working in an under-assembled genome)
}
w.somchrs <- match(somchrs,allchrs)
have.mito <- length(mitochr)>0
C <- length(allchrs)


### Stuff that was in the original package-author's code -- retained for informational purposes only.
##readlen <- length(bam@ranges[[1]])  # this was package authors original code -- would it ever have worked??
##emf.bp <- median(emf[match(somchrs,names(emf))]) 
##shift <- (emf.bp/2)-readlen         # no longer works as advertised, if it ever did.  Unnecessary at any rate.
##covg <- coverage(bam, shift=shift)  # don't want shifted version




message("\nCalculating coverage...")

### Coverage track and peak-calling
covg <- coverage(bam)    # bedgraph data in RLEs
covg <- covg[match(allchrs,names(covg))]  # restrict for chroms of interest

### Calculating bin sums at various scales
#system.time({ bin.sums.10000 <- binStats(covg, 10000, "sum") })  # 25 seconds
system.time({ bin.sums.1000  <- binStats(covg, 1000,  "sum") })   # 3.5 minutes
#system.time({ bin.sums.100   <- binStats(covg, 100,   "sum") })  # > 35 minutes!


message("\nGenerating Statistics...")

### Global coverage stats (entire genome)
totalCoverage <- t(zerofy(merge.table.list(lapply(covg,table))))  # matrix: rows=depths; cols=chroms, cells=N bp on chrom at depth
coverages <- sort(as.numeric(rownames(totalCoverage)))  # coverage depths in ascending numeric order
totalCoverage <- totalCoverage[match(coverages,rownames(totalCoverage)),]  # rows in ascending numeric order
sumCoverage <- rowSums(totalCoverage)  # vector: names=depths, values=N bp in genome at depth
L <- length(sumCoverage)  # N unique depths recorded in genome
genome2plus <- sort(rep(coverages[3:L],times=sumCoverage[3:L]), decreasing=TRUE)  # vector representing bp depths where depth > 2
bpCov <- length(genome2plus)+sumCoverage[2]  # total bp in genome with nonzero depth (= 2+ + 1)

### Sub-global coverage stats (somatic golden-path chroms only)
#totalCoverage2 <- totalCoverage[,colnames(totalCoverage) %in% somchrs]  # restrict to SGP chroms
#coverages2 <- as.numeric(rownames(totalCoverage2))  # coverage depths in ascending numeric order
#sumCoverage2 <- rowSums(totalCoverage2)  # vector: names=depths, values=N bp in genome at depth
#L2 <- length(sumCoverage2)  # N unique depths recorded in genome
#genome2plus2 <- sort(rep(coverages2[3:L],times=sumCoverage2[3:L]), decreasing=TRUE)  # vector representing bp depths where depth > 2
#bpCov2 <- length(genome2plus2)+sumCoverage2[2]  # total bp in genome with nonzero depth (= 2+ + 1)

### Jeff's QC stats
bpAtDepth <- totalCoverage * coverages               # as totalCoverage, but showing sum of sequenced bp contributing to positions at given depth
genome <- sum(totalCoverage)                         # size of genome
top1pct.vals <- genome2plus[1:round(bpCov/100,0)]    # depths vector for top 1% of highest-covered bp in genome, sorted decreasing
top1pct.pos <- length(top1pct.vals)                  # N positions comprising 1%
top1pct.seq <- sum(top1pct.vals)                     # sum of depths of top 1% highest-covered bp
rm(top1pct.vals)

message("pct genome covered: ", round(100*bpCov/genome,2))              # percent nonzero bp in genome
#message("pct nzbp in SGP chroms: ", round(100*bpCov2/bpCov,2))          # nonzero bp in SGP chroms as % of nonzero bp in entire genome
message("pct nzbp in top 1% depths: ", round(100*top1pct.seq/bpAln,2))  # percent of bp in top 1% of coverage, as % all nonzero bp

rm(genome2plus, sumCoverage, coverages)
#rm(genome2plus2,sumCoverage2,coverages2,totalCoverage2)


message("\nHistogramming coverage...")   #, MPD distributions...")

##chr.cov <- apply(totalCoverage2, 2, function(x) log10(as.numeric(detable(x[3:L]))) )  # list with depths vector per chrom, depths > 1 only, in log2
##dhist(chr.cov, points=TRUE)  # density histo
covg.bins <- list( c(0,10), c(10,50), c(50,100), c(100,500), c(500,1E3), c(1E3,5E3), c(5E3,1E4), c(1E4,Inf) )
covg.hist <- list(pctls=percentiles(as.numeric(rownames(totalCoverage))))
covg.hist$full <- rowSums(totalCoverage)
covg.hist$full <- data.frame(VERSION="Complete", DEPTH=as.numeric(names(covg.hist$full)), BP=nameless(covg.hist$full))
covg.hist$bins <- c( covg.hist$full[1,3], sapply(covg.bins, function(x) sum(covg.hist$full[covg.hist$full[,2]>x[1] & covg.hist$full[,2]<=x[2],3]) ))
names(covg.hist$bins) <- sub("-Inf","+",c(0, sapply(covg.bins, function(x) paste0(x[1]+1,"-",format(x[2],scientific=FALSE)) )))
covg.hist$bins <- data.frame(VERSION="Binned", DEPTH=names(covg.hist$bins), BP=c(covg.hist$bins))
write.table(rbind(covg.hist$bins, covg.hist$full), paste0(prefix,".coverage_histo.txt"), sep="\t", quote=FALSE, row.names=FALSE)
##par(mfrow=c(1,2)); barplot(covg.hist$bins[,3],names=covg.hist$bins[,2],las=2); dhist(log10(covg.hist$full[,3]))


if (matepair) {
    ## Analyze mate-pair behavior -- requires building a mate-pair table with /home/apa/local/bin/matePairTable
    
    message("\nPreparing Mate-Pair Data...")
    mpt.file <- sub(".bam$",".matePairTable.txt.gz",bamfile)
    mp.cmd <- paste("/home/apa/local/bin/matePairTable",bamfile,mpt.file)
    message(mp.cmd)
    system(paste("/home/apa/local/bin/matePairTable",bamfile))
    
    message("\nReading Mate-Pair Data...")
    system(paste0("zcat ",mpt.file," | cut -f2,3,9,10,16,17 > ",tmp,"/mpt.txt"))
    mpt <- read.delim(paste0(tmp,"/mpt.txt"), as.is=TRUE)
    
    message("\nAnalyzing Mate-Pair Distances...")
    is.orphan1 <- is.na(mpt[,4])
    is.orphan2 <- is.na(mpt[,2])
    is.orphan <- is.orphan1|is.orphan2
    is.cis <- mpt[,1]==mpt[,3] & !is.orphan
    mp.stats <- c(Orphan=sum(is.orphan), Orphan1=sum(is.orphan1), Orphan2=sum(is.orphan2), Cis=sum(is.cis), Trans=sum(!is.cis), table(mpt[,5]))
    print(mp.stats)
    mpp.stats <- round(100*mp.stats/nrow(mpt),2)
    print(mpp.stats)
    mpd.bins <- list( c(0,50), c(50,100), c(100,500), c(500,1E3), c(1E3,5E3), c(5E3,1E4), c(1E4,1E5), c(1E5,1E6), c(1E6,1E7), c(1E7,Inf) )
    mpd.hist.C <- table(abs(mpt[,6]))
    mpd.hist.C <- data.frame(VERSION="Complete", DIST=as.numeric(names(mpd.hist.C)), MATEPAIRS=c(mpd.hist.C))
    mpd.hist.B <- sapply(mpd.bins, function(x) sum(mpd.hist.C[mpd.hist.C[,2]>x[1] & mpd.hist.C[,2]<=x[2],3]) )
    names(mpd.hist.B) <- sub("-Inf","+",sapply(mpd.bins, function(x) paste0(x[1]+1,"-",format(x[2],scientific=FALSE)) ))
    mpd.hist.B <- data.frame(VERSION="Binned", DIST=names(mpd.hist.B), MATEPAIRS=nameless(mpd.hist.B))
    write.table(rbind(mpd.hist.B, mpd.hist.C), paste0(prefix,".mate_pair_dist_histo.txt"), sep="\t", quote=FALSE, row.names=FALSE)
    lampd <- log2(abs(mpt[,6]))
    png(paste0(prefix,".mate_pair_dist_histo.png"), 600, 600); par(cex=1.2, las=1); dhist(lampd, main="Log2 Mate Pair Dist"); abline(v=log2(1000), col=2); dev.off()
    lampd.saddle <- saddlepoint(lampd)
    message("SADDLE POINT: ",2^lampd.saddle," (linear)")
    ##par(mfrow=c(1,2)); barplot(mpd.hist.B[,3],names=mpd.hist.B[,2],las=2); dhist(log10(mpd.hist.C[,3]))
    
    orient.dist <- table(data.frame(Orient=mpt[,5],Align=sign(mpt[,6])))
    orient.dist.pct <- round(100*orient.dist/nrow(mpt),2)
    
    message("\nAnalyzing Mate-Pair Connectivity...")
    conn <- as.matrix(table(mpt[,c(1,3)]))  # chr-chr connecting mate pairs
    conn <- zerofy(conn[match(allchrs,rownames(conn)),match(allchrs,colnames(conn))])
    dimnames(conn) <- list(allchrs,allchrs)
    minMb <- conn
    mb <- chrsizes[,2]/1E6
    for (i in 1:nrow(minMb)) for (j in 1:ncol(minMb)) minMb[i,j] <- ifelse (i==j, mb[i], min(mb[i],mb[j]))
    conn2 <- round(conn/minMb,0)  # chr-chr connections, per MB of smaller chr ("ppMb")
    conn.diag <- diag(conn)       # cis-connectivity
    conn2.diag <- diag(conn2)     # cis-connectivity, ppMb

    conn.chr <- match(somchrs, rownames(conn))
    chr.conn <- conn[conn.chr,conn.chr]
    chr.conn2 <- conn2[conn.chr,conn.chr]
    lcc.ord <- hclust(dist(zapdiag(log2(chr.conn+1))),"average")$order
    lcc2.ord <- hclust(dist(zapdiag(log2(chr.conn2+1))),"average")$order
    if (have.mito) {
        conn.chrM <- match(mitochr, rownames(conn))
        chrM.conn <- conn[conn.chrM,]    # just the mito row
        chrM.conn2 <- conn2[conn.chrM,]  # "
    }
    
    png(paste0(prefix,".chrom_conn.raw.ord.png"), 700, 600)
    mipu(zapdiag(log2(chr.conn+1)), pal="Reds0", main="Log2 Trans-Connecting Pairs")
    dev.off()
    png(paste0(prefix,".chrom_conn.raw.hc.png"), 700, 600)
    mipu(zapdiag(log2(chr.conn[lcc.ord,lcc.ord]+1)), pal="Reds0", main="Log2 Trans-Connecting Pairs (clustered)")
    dev.off()
    png(paste0(prefix,".chrom_conn.norm.ord.png"), 700, 600)
    mipu(zapdiag(log2(chr.conn2+1)), pal="Reds0", main="Log2 Trans-Connecting Pairs, Norm")
    dev.off()
    png(paste0(prefix,".chrom_conn.norm.hc.png"), 700, 600)
    mipu(zapdiag(log2(chr.conn2[lcc2.ord,lcc2.ord]+1)), pal="Reds0", main="Log2 Connecting Pairs, Norm (clustered)")
    dev.off()
    
    cbars1 <- cbind(Cis=diag(chr.conn), Trans=rowSums(chr.conn)-diag(chr.conn))  #, Mito=chrM.conn[conn.chr])
    cbars1p <- 100*cbars1/rowSums(cbars1)
    ctcol <- c("dodgerblue","tomato")
    png(paste0(prefix,".chrom_cistrans_bars.png"), 900, 900)
    par(mfrow=c(2,1), cex=1.2, las=1)
    barplot(t(cbars1),beside=TRUE,col=ctcol,main="Connectivity Per Chrom, Read Pairs")
    legend("topright", bty="n", pch=15, pt.cex=1.5, col=ctcol, legend=c("Cis","Trans"))
    b <- barplot(t(cbars1p),beside=TRUE,col=ctcol,main="Connectivity Per Chrom, as Pct",ylim=c(0,100))
    text(b[1,], cbars1p[,1], paste0(round(cbars1p[,1],0),"%"), col="white", pos=1, cex=0.9)
    dev.off()
    
    ## Complete connectivity tables (no binning)
     ccon.hist.C <- table(conn.diag)                # cis-connectivities for all chrs
    ccon2.hist.C <- table(conn2.diag)               # cis-connectivities for all chrs, ppMb
     tcon.hist.C <- table(conn[lower.tri(conn)])    # trans-connectivities for all chr pairs
    tcon2.hist.C <- table(conn2[lower.tri(conn2)])  # trans-connectivities for all chr pairs, ppMb
     ccon.hist.C <- data.frame(VERSION="Complete", TYPE="Raw",  CONN="Cis",   VALUE=as.numeric(names(ccon.hist.C)),  MATEPAIRS=nameless(c(ccon.hist.C)))
    ccon2.hist.C <- data.frame(VERSION="Complete", TYPE="Norm", CONN="Cis",   VALUE=as.numeric(names(ccon2.hist.C)), MATEPAIRS=nameless(c(ccon2.hist.C)))
     tcon.hist.C <- data.frame(VERSION="Complete", TYPE="Raw",  CONN="Trans", VALUE=as.numeric(names(tcon.hist.C)),  MATEPAIRS=nameless(c(tcon.hist.C)))
    tcon2.hist.C <- data.frame(VERSION="Complete", TYPE="Norm", CONN="Trans", VALUE=as.numeric(names(tcon2.hist.C)), MATEPAIRS=nameless(c(tcon2.hist.C)))
    
    conn.bins <- list( c(0,10), c(10,50), c(50,100), c(100,500), c(500,1E3), c(1E3,5E3), c(5E3,1E4), c(1E4,5E4), c(5E4,1E5), c(1E5,1E6), c(1E6,Inf) )
    conn.hist.B.names <- sub("-Inf","+",c(0, sapply(conn.bins, function(x) paste0(x[1]+1,"-",format(x[2],scientific=FALSE)) )))

    ## Binned connectivity tables
     ccon.hist.B <- c(ccon.hist.C[1,5],  sapply(conn.bins, function(x) sum( ccon.hist.C[ ccon.hist.C[,4]>x[1] &  ccon.hist.C[,4]<=x[2],5]) ))
    ccon2.hist.B <- c(ccon2.hist.C[1,5], sapply(conn.bins, function(x) sum(ccon2.hist.C[ccon2.hist.C[,4]>x[1] & ccon2.hist.C[,4]<=x[2],5]) ))
     tcon.hist.B <- c(tcon.hist.C[1,5],  sapply(conn.bins, function(x) sum( tcon.hist.C[ tcon.hist.C[,4]>x[1] &  tcon.hist.C[,4]<=x[2],5]) ))
    tcon2.hist.B <- c(tcon2.hist.C[1,5], sapply(conn.bins, function(x) sum(tcon2.hist.C[tcon2.hist.C[,4]>x[1] & tcon2.hist.C[,4]<=x[2],5]) ))
     ccon.hist.B <- data.frame(VERSION="Binned", TYPE="Raw",  CONN="Cis",   VALUE=conn.hist.B.names, MATEPAIRS=ccon.hist.B)
    ccon2.hist.B <- data.frame(VERSION="Binned", TYPE="Norm", CONN="Cis",   VALUE=conn.hist.B.names, MATEPAIRS=ccon2.hist.B)
     tcon.hist.B <- data.frame(VERSION="Binned", TYPE="Raw",  CONN="Trans", VALUE=conn.hist.B.names, MATEPAIRS=tcon.hist.B)
    tcon2.hist.B <- data.frame(VERSION="Binned", TYPE="Norm", CONN="Trans", VALUE=conn.hist.B.names, MATEPAIRS=tcon2.hist.B)
    
    write.table(
        rbind(ccon.hist.B, ccon2.hist.B, tcon.hist.B, tcon2.hist.B, ccon.hist.C, ccon2.hist.C, tcon.hist.C, tcon2.hist.C),
        paste0(prefix,".chrom_conn_histo.txt"), sep="\t", quote=FALSE, row.names=FALSE
    )
    
    matepair.obj <- qw(mpt, is.cis, is.orphan, is.orphan1, is.orphan2, conn, conn2, minMb)  # just the 'big' objects from the mate-pair analysis
    save(list=matepair.obj, file=rdata.mp)  # save to separate RData file
    rm(list=matepair.obj)                   # do not retain in main RData file
}




save.image(rdata.tmp)  # islands and post-islands code is a common source of crashes

message("\nFinding Coverage Islands...")

### Call islands
islands.by.chr <- coverage2islands(covg, aligns, readlen, mindepth=1L, prefix="RChipQC_peak_")
Nislands <- sapply(islands.by.chr,length)
islands.by.chr <- islands.by.chr[Nislands>0]
Nislands <- sum(Nislands)
CI <- length(islands.by.chr)  # some chrs may have been removed now

message("\nPreparing Peaks Table...")

### "Peaks" = islands above minimum height (if sample is really bad, this may result in some max-peaks not being reported in all-peaks table!)
all.peaks <- as.data.frame(do.call(rbind, lapply(islands.by.chr, function(x) cbind(gr2bed(x), Width=width(x), mcols(x)) )))
all.peaks <- all.peaks[all.peaks$MaxHeight>=minht,]  # peaks w/ reportable height

### Finding outlier peaks per chromosome -- tallest, widest, most-reads
get.max.peaks <- function(i, key) {
    ## requires gr2bed() from apa_tools.R
    gr <- islands.by.chr[[i]]
    at <- wd <- width(gr)
    if (key!="Width") at <- attributes(gr)$elementMetadata[[key]]
    m <- max(at)
    w <- if (key=="Width") { which(wd==m) } else { which(attributes(gr)$elementMetadata[[key]]==m) }
    as.data.frame(cbind(gr2bed(gr[w]),Width=wd[w],attributes(gr)$elementMetadata[w,]))
}
pre.max.peaks <- list(
    tallest=do.call(rbind, lapply(1:CI, get.max.peaks, "MaxHeight")),
    widest=do.call(rbind, lapply(1:CI, get.max.peaks, "Width")),
    largest=do.call(rbind, lapply(1:CI, get.max.peaks, "TotalBp"))
)
max.peaks <- unique(do.call(rbind, pre.max.peaks))
max.peaks <- cbind(max.peaks, Crowns=sapply(max.peaks[,4], function(n) paste(names(pre.max.peaks)[sapply(pre.max.peaks,function(x) n %in% x[,4])],collapse=",") ))
max.peaks <- rownameless(max.peaks[order(max.peaks[,4]),])
colnames(max.peaks)[1:4] <- c("Chrom","Start","End","Name")
max.values <- c("MaxHeight","Width","TotalBp"); names(max.values) <- max.values
max.values <- sapply(max.values, function(x) max(max.peaks[[x]]) )
rm(pre.max.peaks)



message("\nCross-Correlating on Genome...")

### Get best shift size based on cross-correlation, then second-best shift size, in case the first equals the read length
CCF <- cross.correlate(bamfile2, readlen, 1000)
if (CCF$best$shift.size[1]==readlen) {
    message("best x-cor shift size equals the read length: x-cor failed!")
    if (CCF$best$shift.size[2]>readlen*2+1) {
        message("  but, second-best x-cor shift size found at ",CCF$best$shift.size[2],"bp")
    } else {
        message("  second-best x-cor search also failed!")
    }
} else {
    message("best x-cor shift size ",round(CCF$best$xcor.value[1],4)," at ",CCF$best$shift.size[1],"bp")
}



message("\nCross-Correlating on Top Peaks...")

### Re-do cross-correlation using only the 1000 highest peaks that are wider than 2*readlen -- similar to MACS2 modeling
top.peaks <- all.peaks[all.peaks$Width>2*readlen,]
top.peaks <- head(top.peaks[rev(order(top.peaks$MaxHeight)),1:4], 1000)
top.CCF <- cross.correlate(bamfile2, readlen, 1000, top.peaks)
if (top.CCF$best$shift.size[1]==readlen) {
    message("best top-peak x-cor shift size equals the read length: top-peak x-cor failed!")
    if (top.CCF$best$shift.size[2]>readlen*2+1) {
        message("  but, second-best top-peak x-cor shift size found at ",top.CCF$best$shift.size[2],"bp")
    } else {
        message("  second-best top-peak x-cor search also failed!")
    }
} else {
    message("best top-peak x-cor shift size ",round(top.CCF$best$xcor.value[1],4)," at ",top.CCF$best$shift.size[1],"bp")
}



message("\nWriting and Saving...")

### Build and write stats table
sci.cov <- function(x) ifelse(x>=0.01, sprintf("%0.3f",x), format(x,scientific=TRUE,digits=3))
sci.pct <- function(x) ifelse(x>=0.0001, sprintf("%0.3f",100*x), format(100*x,scientific=TRUE,digits=3))  # percent 'x' must be on [0,1] NOT [0,100]
comma <- function(x) format(x,big.mark=",",scientific=FALSE)

stats.table <- do.call(rbind, list(
    c("Bam file", bamfile2 ),
    c("Read Length", readlen ),
    c("Best X-Cor Shift", CCF$best$shift.size[1] ),
    c("Best X-Cor Value", round(CCF$best$xcor.value[1],4) ),
    c("2nd-Best X-Cor Shift", CCF$best$shift.size[2] ),
    c("2nd-Best X-Cor Value", round(CCF$best$xcor.value[2],4) ),
    c("Top-Peak Best X-Cor Shift", top.CCF$best$shift.size[1] ),
    c("Top-Peak Best X-Cor Value", round(top.CCF$best$xcor.value[1],4) ),
    c("Top-Peak 2nd-Best X-Cor Shift", top.CCF$best$shift.size[2] ),
    c("Top-Peak 2nd-Best X-Cor Value", round(top.CCF$best$xcor.value[2],4) ),
    c("Alignments", aligns ),
    c("Genomic Coverage", sci.cov(as.numeric(aligns)*as.numeric(readlen)/genome) ),
    c("Genome Bp", genome ),
    c("Nonzero Pos", bpCov ),
    c("Nonzero Pos %", sci.pct(bpCov/genome) ),
    c("Total Bp Aligned", bpAln ),
    c("Pos in Top 1% of Depths", top1pct.pos ),
    c("Bp in Top 1% Pos", top1pct.seq ),
    c("% Bp Aligned in Top 1% Pos",sci.pct(top1pct.seq/bpAln) ),
    c("Unique Alignment Starts", upos ),
    c("Unique Alignment %", sci.pct(upos/aligns) ),
    c("Min Height of Reportable Island", minht ),
    c("Total Islands", Nislands ),
    c("Islands >= Min Height", nrow(all.peaks) ),
    c("% Islands >= Min Height", sci.pct(nrow(all.peaks)/Nislands) ),
    c("Max Island Height", max.values[["MaxHeight"]] ),
    c("Max Island Width", max.values[["Width"]] ),
    c("Max Island Bp", max.values[["TotalBp"]] )
))
write.table(stats.table, paste(prefix,"stats.txt",sep="."), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

### Output peak-calls tables
write.table(all.peaks, gzfile(paste0(prefix,".peaks.txt.gz")), sep="\t", quote=FALSE, row.names=FALSE)
write.table(max.peaks, paste(prefix,"chrom_max_peaks.txt",sep="."), sep="\t", quote=FALSE, row.names=FALSE)
rm(all.peaks,max.peaks)

### Cleanup, save and exit
if (!keepbam) rm(bam)
rm.apa.tools()
save.image(rdata)
system(paste("rm -f",rdata.tmp))
message(bamfile," complete!\n")
quit()




#stats.table <- do.call(rbind, list(
#    c("Bam file", bamfile2),
#    c("Read Length", comma(readlen)),
#    c("Best X-Cor Shift", CCF$best$shift.size[1]),
#    c("Best X-Cor Value", round(CCF$best$xcor.value[1],4)),
#    c("2nd-Best X-Cor Shift", CCF$best$shift.size[2]),
#    c("2nd-Best X-Cor Value", round(CCF$best$xcor.value[2],4)),
#    c("Top-Peak Best X-Cor Shift", top.CCF$best$shift.size[1]),
#    c("Top-Peak Best X-Cor Value", round(top.CCF$best$xcor.value[1],4)),
#    c("Top-Peak 2nd-Best X-Cor Shift", top.CCF$best$shift.size[2]),
#    c("Top-Peak 2nd-Best X-Cor Value", round(top.CCF$best$xcor.value[2],4)),
#    c("Alignments", comma(aligns)),
#    c("Genomic Coverage", sci.cov(as.numeric(aligns)*as.numeric(readlen)/genome) ),
#    c("Genome Bp", comma(genome)),
#    c("Nonzero Pos", comma(bpCov)),
#    c("Nonzero Pos %", sci.pct(bpCov/genome) ),
#    c("Total Bp Aligned", comma(bpAln)),
#    c("Pos in Top 1% of Depths", comma(top1pct.pos)),
#    c("Bp in Top 1% Pos", comma(top1pct.seq)),
#    c("% Bp Aligned in Top 1% Pos",sci.pct(top1pct.seq/bpAln) ),
#    c("Unique Alignment Starts", comma(upos)),
#    c("Unique Alignment %", sci.pct(upos/aligns) ),
#    c("Min Height of Reportable Island", comma(minht)),
#    c("Total Islands", comma(Nislands)),
#    c("Islands >= Min Height", comma(nrow(all.peaks))),
#    c("% Islands >= Min Height", sci.pct(nrow(all.peaks)/Nislands)),
#    c("Max Island Height", comma(max.values[["MaxHeight"]])),
#    c("Max Island Width", comma(max.values[["Width"]])),
#    c("Max Island Bp", comma(max.values[["TotalBp"]]))
#))




