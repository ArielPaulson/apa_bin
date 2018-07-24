#!/usr/bin/env Rscript

source("/home/apa/apa_tools.R")

ca <- commandArgs(trailing=TRUE)
opdir <- fixpath(ca[1])
sample <- ca[2]
geno <- ca[3]
cores <- as.numeric(ca[4])
outpref <- ca[5]     # OPTIONAL
idr1file <- ca[6]    # OPTIONAL
idr2file <- ca[7]    # OPTIONAL

## opdir="/n/core/Bioinformatics/analysis/Krumlauf/bdk/cbio.bdk.115/data/pipeline/output/"; sample="Sox2_RA"; geno="mm10"; cores=60; outpref=NA; idr1file=NA; idr2file=NA
## opdir="/n/core/Bioinformatics/analysis/Krumlauf/nip/cbio.nip.101/data/chipseq/output/"; sample="Hoxa1_D1_Spondin"; geno="mm10"; cores=60; outpref=NA; idr1file=NA; idr2file=NA 

sampdir <- paste0(opdir,sample,"/")
bamdir <- paste0(sampdir,"bams/")
idrdir <- paste0(sampdir,"idr/")
if (!dir.exists(idrdir)) stop("Expected IDR directory '",idrdir,"' does not exist!\n")    
idr1dir <- paste0(idrdir,"IDR1/")
idr2dir <- paste0(idrdir,"IDR2/")
compdir <- paste0(idrdir,"comp/")
if (!dir.exists(compdir)) system(paste("mkdir",compdir))

if (is.na(outpref)) outpref <- "idr1_idr2_compare."
if (!grepl("\\.$",outpref)) outpref <- paste0(outpref,".")
outpref <- paste0(compdir,outpref)
RData <- paste0(outpref,"RData")

checkpoint <- 0
#IM(RData)
#IM(file.exists(RData))
#outpref1=outpref
if (file.exists(RData)) {
#    message("Reloading from last checkpoint...")
#    load(RData)
}
#outpref=outpref1
#stop(checkpoint," ",outpref)


if (checkpoint < 1) {
    
    spacer <- 5
    window <- 5000
    IDRs <- c("IDR1","IDR2")
    sigs  <- c(0.01, 0.05)  # default significant limits for each IDR
    sigs2 <- c(0.05, 0.10)  # new relaxed thresholds to use in combination with FC>=5
    
    if (is.na(idr1file)) idr1file <- paste0(idr1dir,"true.IDR-overlapped-peaks.expanded.txt")
    if (is.na(idr2file)) idr2file <- paste0(idr2dir,"true.idrValues.expanded.txt")
    
    ## p-value breakpoints; classification uses (N,N+1]
    ## if any p=0, they will get included in tranch 1
    ## tranch 11 is only for those with p=1 exactly (p>0.9999999999999, basically)
    ## tranch 12, if it appears, is only for NA values (i.e. peak found only in the other IDR)
    p.tranches <- c(0,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.25,0.5,0.9999999999999,1)
    P <- length(p.tranches)  # not -1, because we will add the NA tranch later
    p.tranch.cols <- c("red","darkorange","gold2","chartreuse","green3","cyan","dodgerblue","blue","purple3","sienna4","black","grey")
    p.tranch.names <- sapply(1:(P-1), function(i) paste(p.tranches[i],"< p <=",p.tranches[i+1],collapse="") )
    p.tranch.names[1] <- sub("<","<=",p.tranch.names[1])
    for (i in 1:2) p.tranch.names[i+4] <- paste0(p.tranch.names[i+4]," (IDR",i," default sig limit)")
    p.tranch.names[10] <- sub("= 0.9+"," 1",p.tranch.names[10])
    p.tranch.names[11] <- "p = 1"
    p.tranch.names[12] <- "Not Found"
    short.p <- 7    # short-list heatmap: tranches 1-short.p only
    short.p.cols <- p.tranch.cols[1:short.p]
    
    idr1 <- read.delim(idr1file,as.is=TRUE)
    idr2 <- read.delim(idr2file,as.is=TRUE)
    
    ## Create BED datasets
    idr1.bed <- sort.bed(data.frame(idr1[,1:4], score=10^-idr1[,12], strand="+"))
    idr2.bed <- sort.bed(data.frame(idr2[,1:4], score=10^-idr2[,12], strand="+"))
    colnames(idr1.bed) <- colnames(idr2.bed) <- c("chrom","start","end","name","score","strand")
    ## write BED files
    idr1.bedfile <- paste0(compdir,"IDR1.bed")
    idr2.bedfile <- paste0(compdir,"IDR2.bed")
    write.bed(idr1.bed, idr1.bedfile)
    write.bed(idr2.bed, idr2.bedfile)
    ## run mergePeaks
    mp.cmd <- paste0("/home/apa/local/bin/mergePeaks -oc ",compdir,"merge -n ",sample,"_IDR_merge -f 500 -g ",geno," ",idr1.bedfile," ",idr2.bedfile)
    message(mp.cmd)
    system(mp.cmd)
    merge <- read.mergePeaks(paste0(compdir,"merge"))
    print(merge[1:3])
    
    ## IDR1-IDR2 merged BED file
    final.bed <- merge$peaks[,1:6]
    final.xbedfile <- paste0(compdir,"merge/final.5k.bed")
    final.xbed <- write.heatmap.bed(final.bed, final.xbedfile, geno, window=window)
    
    ## midpoint+-5k BED files
    idr1.xbedfile <- paste0(compdir,"IDR1.5k.bed")
    idr2.xbedfile <- paste0(compdir,"IDR2.5k.bed")
    ## create/write midpoint+-5k BED datasets
    idr1.xbed <- write.heatmap.bed(idr1.bed[,1:6], idr1.xbedfile, geno, window=window)
    idr2.xbed <- write.heatmap.bed(idr2.bed[,1:6], idr2.xbedfile, geno, window=window)
    
    ## initial BED data, sorted increasing by p-value
    idr1.ordbed <- idr1.bed[order(idr1.bed[,5]),]
    idr2.ordbed <- idr2.bed[order(idr2.bed[,5]),]
    ## any peaks from IDR1 which WERE found in IDR2, and vice versa
    idr1.common <- unique(merge$source2merge$Merge[merge$source2merge$Source=="IDR1"])
    idr2.common <- unique(merge$source2merge$Merge[merge$source2merge$Source=="IDR2"])
    ## any peaks from IDR1 which were NOT found in IDR2, and vice versa
    idr1.unique <- idr1.ordbed[idr1.ordbed[,4] %in% merge$source2merge$Peak[!(merge$source2merge$Merge %in% idr2.common)],]
    idr2.unique <- idr2.ordbed[idr2.ordbed[,4] %in% merge$source2merge$Peak[!(merge$source2merge$Merge %in% idr1.common)],]
    
    ## total-IDR datasets, from IDR1 or IDR2 point of view
    ## i.e. IDR1 = IDR1 peaks + unique-to-IDR2 peaks; IDR2 = IDR2 peaks + unique-to-IDR1 peaks
    ## 7th col is p-value tranch
    idr1.ordbed2 <- cbind(idr1.ordbed,sig=0)
    idr2.ordbed2 <- cbind(idr2.ordbed,sig=0)
    if (nrow(idr2.unique)>0) idr1.ordbed2 <- rbind(idr1.ordbed2, cbind(idr2.unique,sig=NA))
    if (nrow(idr1.unique)>0) idr2.ordbed2 <- rbind(idr2.ordbed2, cbind(idr1.unique,sig=NA))
    
    ## assign p-value tranches
    for (p in 1:P) {
        ## tranches are pmin < pvalue <= pmax
        pmin <- ifelse(p==1,-1,p.tranches[p])  # since breakpoint 1 == 0, make sure to include any (super-rare!) p=0 in it...
        pmax <- p.tranches[p+1]
        w1 <- which(!is.na(idr1.ordbed2$sig) & idr1.ordbed2$score>pmin & idr1.ordbed2$score<=pmax)
        w2 <- which(!is.na(idr2.ordbed2$sig) & idr2.ordbed2$score>pmin & idr2.ordbed2$score<=pmax)
        idr1.ordbed2$sig[w1] <- p
        idr2.ordbed2$sig[w2] <- p
    }
    idr1.ordbed2$sig[is.na(idr1.ordbed2$sig)] <- 12  # tranch 12, if any exist
    idr2.ordbed2$sig[is.na(idr2.ordbed2$sig)] <- 12
    
    ## midpoint+-5k BED files for total-IDR
    idr1.xordbed2file <- paste0(compdir,"IDR1_ord2.5k.bed")
    idr2.xordbed2file <- paste0(compdir,"IDR2_ord2.5k.bed")
    ## create/write midpoint+-5k BED datasets for total-IDR
    idr1.xordbed2 <- write.heatmap.bed(idr1.ordbed2[,1:6], idr1.xordbed2file, geno, window=window, no.mito=TRUE, debug=TRUE)
    idr2.xordbed2 <- write.heatmap.bed(idr2.ordbed2[,1:6], idr2.xordbed2file, geno, window=window, no.mito=TRUE, debug=TRUE)
    
    ## bed data MUST be restricted to only those rows in the heatmaps!!!
    hm.bed <- list(
        IDR1=list(
            bed=idr1.ordbed2[idr1.ordbed2$name %in% idr1.xordbed2$name[idr1.xordbed2$OK],],
            xbed=cbind(idr1.xordbed2,FC=idr1$FoldChange[match(idr1.xordbed2$name,idr1$Name)])[idr1.xordbed2$OK,],  # add FC as col 7
            xbedf=idr1.xordbed2file
        ),
        IDR2=list(
            bed=idr2.ordbed2[idr2.ordbed2$name %in% idr2.xordbed2$name[idr2.xordbed2$OK],],
            xbed=cbind(idr2.xordbed2,FC=idr2$FoldChange[match(idr2.xordbed2$name,idr2$Name)])[idr2.xordbed2$OK,],  # add FC as col 7
            xbedf=idr2.xordbed2file
        )
    )
    
    bams <- c(
        pool=paste0(bamdir,"pooled.bam"),
        true1=paste0(bamdir,"true-1.bam"),
        true2=paste0(bamdir,"true-2.bam"),
        input=paste0(bamdir,"input.bam")
    )
    B <- length(bams)
    for (i in 1:B) if (!file.exists(sub("bam$","idxstats.txt",bams[i]))) system(paste("samtools idxstats",bams[i],">",sub("bam$","idxstats.txt",bams[i])))
    reads <- sapply(bams, function(x){ y=read.delim(sub("bam$","idxstats.txt",x),as.is=TRUE,header=FALSE); sum(y[y[,1]!="*",3]) })
    bams.exist <- file.exists(bams)
    
    checkpoint <- 1
    message("Saving Checkpoint 1...")
    save.image(RData)  # initial savepoint
}



if (checkpoint < 2) {
    
    print(data.frame(BAM=bams,EXISTS=bams.exist,READS=reads))
    if (any(!bams.exist)) stop("Some bam files are missing!\n")
    
    require(CoverageView)
    
    mats <- logrpm <- LFC <- new.list(IDRs,elem=new.list(names(bams)))
    for (i in 1:2) {
        for (j in 1:B) {
            message(Sys.time()," Processing ",IDRs[i]," ",bams[j])
            mats[[i]][[j]] <- suppressMessages(t(cov.matrix(CoverageBamFile(bams[j],reads_mapped=reads[j]), coordfile=hm.bed[[i]]$xbedf, no_windows=100, num_cores=cores)))
            rownames(mats[[i]][[j]]) <- hm.bed[[i]]$xbed$name
            logrpm[[i]][[j]] <- round(log2((1E6*mats[[i]][[j]]/reads[j])+1),3)
        }
    }
    
    for (i in 1:2) {
        LFC[[i]] <- LFC[[i]][1:(B-1)]
        for (j in 1:(B-1)) LFC[[i]][[j]] <- log2( (2^logrpm[[i]][[j]]) / (2^logrpm[[i]][[4]]) )
    }
    
    checkpoint <- 2
    message("Saving Checkpoint 2...")
    save.image(RData)  # post-CoverageView savepoint
}



if (checkpoint < 3) {
    
    attw <- 20
    fccw <- 10
    main <- paste0( rep(c("LogRpm","Zscore","LFC"),times=c(B,B,B-1)), ".", rep(names(mats[[1]]),3)[2:(3*B)-1] )
    spacers <- lapply(mats, function(x) matrix(NA,nrow(x[[1]]),spacer) )
    att.cols <- spacers; for (i in 1:2) att.cols[[i]] <- rowname(colrep(hm.bed[[i]]$bed$sig,attw),hm.bed[[i]]$bed$name)  # add att.cols = p.tranches // USE $BED NOT $XBED
    
    block.norm <- function(x) {  # 'x' is a cbound matrix of heatmap columns, OF THE SAME TYPE i.e. all logrpm, zscore, or LFC, NOT mixed
        x <- x-min(x,na.rm=TRUE)             # min val -> 0
        x <- x/quantile(x,0.999,na.rm=TRUE)  # max val -> some percent of the 99.9th pctl (may exceed 1)
        threshold(x,1,"gt")                  # max val -> 1
    }
    p.order <- function(x, i, t0=TRUE) {   # 'x' is a list of heatmap matrices; 'i' is 1 or 2
        x <- if (t0) { lapply(x, function(y) zerofy(threshold(y,0,"lt")) ) } else { lapply(x, function(y) zerofy(y) ) }
        y <- if (length(x)>1) { coverage.matrix.cbind(x, spacer) } else { cbind(x[[1]],spacer) }
        y <- block.norm(y)
        att <- att.cols[[i]][real(match(rownames(y),rownames(att.cols[[i]]))),]   # restrict to y rows // put in y order
        ##IM(i,dim(y),dim(spacers[[i]]),dim(att))
        z <- cbind(y,spacers[[i]],att)
        z[real(match(hm.bed[[i]]$bed$name,rownames(y))),]    # put in p-sorted order (use row order of $bed, not $xbed)
    }
    
    hm <- hm.att <- hm.fcatt <- shortlist <- new.list(names(LFC))
    for (i in 1:2) {
        i.idr <- if (i==1) { idr1 } else { idr2 }
        fcm <- rowname(i.idr[,c("FoldChange","Rep1.FoldChange","Rep2.FoldChange")], i.idr$Name)
        ##fcm <- rowname(i.idr[,c("Pvalue","Rep1.Pvalue","Rep2.Pvalue")], i.idr$Name)
        ##fcm <- rowname(colrep(i.idr$IDR.Pglobal,3), i.idr$Name)
        fcm <- fcm[,rep(1:3,each=fccw)]
        fcm <- as.matrix(fcm[match(hm.bed[[i]]$xbed$name, rownames(fcm)),])
        print(listLengths(LFC[[i]]))
        hm[[i]] <- list( p.order(logrpm[[i]],i,FALSE), p.order(lapply(logrpm[[i]],zscore),i), p.order(LFC[[i]],i), p.order(list(fcm),i) )
        hm[[i]][[3]] <- hm[[i]][[3]][,1:(ncol(hm[[i]][[3]])-attw)]  # drop att.cols from LFC part; adjacent fcm will bring them back (after fcm)
        offset <- 0
        for (j in 1:4) {
            nc <- ncol(hm[[i]][[j]])
            if (j==3) {
                ## no att cols
            } else {
                hm.att[[i]] <- c(hm.att[[i]], offset+c((nc-attw):nc))
                if (j==4) hm.fcatt[[i]] <- (offset+1):(offset+ncol(fcm))
            }
            offset <- offset+nc+spacer
        }
        hm[[i]] <- coverage.matrix.cbind(hm[[i]], spacer)
        hm[[i]] <- cbind( att.cols[[i]][real(match(rownames(hm[[i]]),rownames(att.cols[[i]]))),], spacers[[i]], hm[[i]])
        hm.att[[i]] <- c(1:attw, hm.att[[i]]+attw+spacer)
        hm.fcatt[[i]] <- hm.fcatt[[i]]+attw+spacer
        shortlist[[i]] <- which(falsify(hm[[i]][,ncol(hm[[i]])]<=short.p))
    }
    
    ## based on an 11-matrix heatmap (4 logrpm, 4 zscore, 3 LFC)
    ## which was written to PNG @ 100 pixels wide, then measured by hand
    ##  1/22 of width = margin
    ## 21/22 of width = heatmap
    ## ~45.5 pixels per 22nd
    ## heatmap cols: 1100 (88%) matrix, 80 (6.4%) att, 30 (2.4%) req spacer, 40 (3.2%) opt spacer
    ## N matrices = N*100+(N-1)*5
    ## N heatmap cols wide = N*100+(N-1)*5 + 80+30
    ## 11 matrices @ 1 pixel per column = (11*100+(11-1)*5 + 80+30) / (21/22) = 1320 PNG width, which is acceptable
    
    N <- 3*B-1
    message("N = ",N)
    imgwd  <- round( (N*100+(N-1)*5+110)/(21/22) * 0.76 , 0 )  # scaling down by 24% (e.g. 1320 -> 1003)
    imght  <- round(nrow(hm[[1]])*0.08,0)
    imght  <- 200 + ifelse(imght<=500, 500, ifelse(imght<=1800, imght, 1800)) 
    imght2 <- 200 + sapply(shortlist, function(x) max(c( 500, round(length(x)*0.1,0) )) )  # LENGTH 2
    
    checkpoint <- 3
    message("Saving Checkpoint 3...")
    save.image(RData)  # final savepoint
}


## Image writing is past all checkpoints

message("Imaging...")

png(paste0(compdir,"color_key.png"), 550, 500)
viewPalette(p.tranch.cols,labels=p.tranch.names,values=FALSE,cex=1.2)
dev.off()

idr1vs <- idr1[,c("Pvalue","IDR.Pglobal","FoldChange")]
idr2vs <- idr2[,c("Pvalue","IDR.Pglobal","FoldChange")]
idr1sig <- idr1vs[,2]>=-log10(sigs2[1]) & idr1vs[,3]>=5
idr2sig <- idr2vs[,2]>=-log10(sigs2[2]) & idr2vs[,3]>=5

png(paste0(compdir,"idr_vs_p_fc.png"), 1000, 1000)
par(mfrow=c(2,2), las=1, cex=1.2)

plot(idr1vs[,3:2], col=0, xlab="FoldChange", ylab="IDR PGlobal", main="IDR1 vs Mean FC")
points(idr1vs[!idr1sig,3:2], col=1); points(idr1vs[idr1sig,3:2], col=3)
abline(v=c(2,5), lty=2:1, lwd=2, col=5); abline(h=-log10(c(sigs[1],sigs2[1])), lty=2:1, lwd=2, col=2)
legend("bottomright", bty="n", col=c(3,1), pch=1, legend=c(sum(idr1sig),sum(!idr1sig)))

plot(idr2vs[,3:2], xlab="FoldChange", ylab="IDR PGlobal", main="IDR2 vs Mean FC")
points(idr2vs[!idr2sig,3:2], col=1); points(idr2vs[idr2sig,3:2], col=3)
abline(v=c(2,5), lty=2:1, lwd=2, col=5); abline(h=-log10(c(sigs[2],sigs2[2])), lty=2:1, lwd=2, col=2)
legend("bottomright", bty="n", col=c(3,1), pch=1, legend=c(sum(idr2sig),sum(!idr2sig)))

plot(idr1vs[,1:2], xlab="Pvalue", ylab="IDR PGlobal", main="IDR1 vs Mean P")
points(idr1vs[!idr1sig,1:2], col=1); points(idr1vs[idr1sig,1:2], col=3)
abline(v=-log10(c(0.05,1E-4)), lty=2:1, lwd=2, col=5); abline(h=-log10(c(sigs[1],sigs2[1])), lty=2:1, lwd=2, col=2)
legend("bottomright", bty="n", col=c(3,1), pch=1, legend=c(sum(idr1sig),sum(!idr1sig)))

plot(idr2vs[,1:2], xlab="Pvalue", ylab="IDR PGlobal", main="IDR2 vs Mean P")
points(idr2vs[!idr2sig,1:2], col=1); points(idr2vs[idr2sig,1:2], col=3)
abline(v=-log10(c(0.05,1E-4)), lty=2:1, lwd=2, col=5); abline(h=-log10(c(sigs[2],sigs2[2])), lty=2:1, lwd=2, col=2)
legend("bottomright", bty="n", col=c(3,1), pch=1, legend=c(sum(idr2sig),sum(!idr2sig)))
dev.off()

for (i in 1:2) {
    nsig <- sum(hm.bed[[i]]$xbed$score <= sigs[i])
    nsig2 <- sum(hm.bed[[i]]$xbed$score <= sigs2[i] & hm.bed[[i]]$xbed$FC >= 5)
    sig.blurb <- paste0("(",nsig," <= ",sigs[i]," | ",nsig2," <= ",sigs2[i]," & FC >= 5)")
    attrib <- list(list(cols=hm.att[[i]], palette=c()), list(cols=hm.fcatt[[i]], palette=palettizer("KW",256)))
    
    IM(i,"full")
    x <- nameless(hm[[i]])
    attrib[[1]]$palette <- p.tranch.cols[sort(unique(x[,ncol(x)]))]
    subtitle <- paste0(sample,": ",nrow(x)," IDR",i," peaks ",sig.blurb)
    png(paste0(compdir,"IDR",i,"_full.png"), imgwd, imght)
    mipu(x, pal="Reds0", pmar=c(6,1), cex=1.2, main=c(main,"FC"), stagger.main=2, sub=subtitle, show.scale=FALSE, attrib=attrib)
    dev.off()
    
    IM(i,"short")
    x <- nameless(hm[[i]][shortlist[[i]],])
    attrib[[1]]$palette <- p.tranch.cols[sort(unique(x[,ncol(x)]))]
    subtitle <- paste0(sample,": ",nrow(x)," short-list IDR",i," peaks ",sig.blurb)
    png(paste0(compdir,"IDR",i,"_short.png"), imgwd, imght2[[i]])
    mipu(x, pal="Reds0", pmar=c(6,1), cex=1.2, main=c(main,"FC"), stagger.main=2, sub=subtitle, show.scale=FALSE, attrib=attrib)
    dev.off()
}


   # mipu(x[1:1000,], pal="Reds0", pmar=c(1,1), cex=1.2, main=main, stagger.main=2, sub=subtitle, show.scale=FALSE, att.cols=hm.att[[i]], att.pal=att.cols)


