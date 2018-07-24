#!/usr/bin/env Rscript

library(edgeR)
library(WriteXLS)

cmds <- commandArgs(trailing=TRUE)
input <- cmds[[1]]
genes <- cmds[[2]]
design <- cmds[[3]]

roots <- c(unix="/n/projects", windows="U:")
path <- roots[match(.Platform$OS.type,names(roots))]
source(paste(path,"apa/R/apa_tools.R",sep="/"))

inp <- read.delim(input, as.is=TRUE)
gdat <- read.delim(genes, as.is=T)
header <- slice.list(read.xlsx.simple("BySample.xlsx", header=FALSE), 1, "rows")
dat <- read.xlsx.simple("BySample.xlsx")
for (i in 1:length(dat)) colnames(dat[[i]]) <- unlist(header[[i]][1,])
NC <- ncol(dat[[1]])

loc.s <- 2:NC      # local columns, summary page
loc.n <- loc.s+2   # local columns, ncRNA pages

ncrna.dat <- cbind(CLASS=qw(rRNA,snoRNA,snRNA,tRNA,miRNA,misc_RNA,Mt_rRNA,Mt_tRNA),COLOR=c(3,4,5,6,2,7,1,8))
ncrna.dat <- ncrna.dat[ncrna.dat[,1] %in% dat[[1]][2:nrow(dat[[1]]),1],]  # remove any classes not present in analysis

x <- which(apply(is.na(dat[[1]]),1,all))
pct.offset <- x[which(diff(x)==1)]+1


### TOTAL READS: COUNTS AND PERCENTS

read.plot <- function(rows, cols, main) {
    colors <- c(1,4,3,7,2)
    x <- zerofy.numeric(as.matrix(dat[[1]][rows,cols]))
    rownames(x) <- dat[[1]][rows,1]
    rownames(x)[1] <- "TOTAL"
    N <- length(cols)
    barplot(x, beside=T, col=colors, main=main, xlim=c(0,length(x)+length(cols)+6))
    legend("topright", pt.cex=1.5, pch=15, bty="n", col=colors, legend=rownames(x))
}

png("total_alignments.png", max(c(500, 250+150*length(loc.s))), 1000)
par(mar=c(5,5,4,2), mfrow=c(2,1), cex=1.2, las=1, xaxs="i")
ord <- c(1,3:5,2)
read.plot(ord, loc.s, "READS: COUNTS")
read.plot(ord+pct.offset, loc.s, "READS: PERCENTS")
dev.off()


### TOTAL ALIGNMENTS: COUNTS AND PERCENTS

ncRNA.plot <- function(rows, cols, main, log) {
    x <- zerofy.numeric(as.matrix(dat[[1]][rows,cols]))
    rownames(x) <- dat[[1]][rows,1]
    if (log) x <- log10(x)  # cheap way of doing this
    b <- barplot(x, beside=T, col=ncrna.dat[,2], main=main, xlim=c(0,length(x)+length(cols)+6))
    legend("topright", pt.cex=1.5, pch=15, bty="n", col=ncrna.dat[,2], legend=rownames(x))
    invisible(b)
}

png("ncRNA_alignments.png", max(c(500, 200+150*length(loc.s))), 1000)
par(mar=c(5,3,4,2), mfrow=c(2,1), cex=1.2, las=1, xaxs="i")
ord <- match(ncrna.dat[,1],dat[[1]][1:pct.offset,1])
ncRNA.plot(ord, loc.s, "READS: LOG10 COUNTS", TRUE)
ncRNA.plot(ord+pct.offset, loc.s, "READS: PERCENTS", FALSE)
dev.off()


### NCRNA COUNTS BY CLASS

class.ct <- suppressWarnings(lapply(mat.split(dat[[2]][,c(2,loc.n)], dat[[2]][,1]), function(x){ zerofy.numeric(rownamed.matrix(x)) }))
class.ct <- class.ct[match(ncrna.dat[,1],names(class.ct))]
class.rp <- suppressWarnings(lapply(mat.split(dat[[3]][,c(2,loc.n)], dat[[3]][,1]), function(x){ zerofy.numeric(rownamed.matrix(x)) }))
class.rp <- class.rp[match(ncrna.dat[,1],names(class.rp))]
class.n <- sapply(class.ct, nrow)
class.n

count.plot <- function(classes, rpkm=FALSE) {
    class.range <- matrix(0, nrow(ncrna.dat), 2, F, list(ncrna.dat[,1],qw(Min.Exp,Max.Exp)))
    par(mfrow=c(2,4), cex=1.3, las=1, xaxs="i", family="mono")
    for (i in 1:nrow(ncrna.dat)) {
        if (rpkm) {
            x <- log2(classes[[i]]+1)
            main <- paste(names(classes)[i], "Log2 RPKMs")
        } else {
            x <- log10(classes[[i]]+1)
            main <- paste(names(classes)[i], "Log10 Reads")
        }
        axseq <- seq(0, max(x), length=10)
        xmin <- min(1, 1-nrow(x)*0.01)
        null.plot(xlim=c(xmin,nrow(x)), ylim=c(0,max(x)), main=main, xlab="Individual Genes", ylab="")
        for (k in 1:ncol(x)) lines(1:nrow(x), rev(sort(x[,k])), col=k+1)
        range.exp <- range(apply(x,2,function(y){ sum(y>0) }))
        class.range[i,] <- range.exp
        abline(v=range.exp[2])
        text.pos <- ifelse(range.exp[2]<nrow(x)/2, 4, 2)
        max.exp <- nchar(max(range.exp))
        fmt <- paste("Per-Sample Expression:\n%",max.exp,"i Max (%0.2f%s)\n%",max.exp,"i Min (%0.2f%s)",sep="")
        exp.text <- sprintf(fmt, range.exp[2], 100*range.exp[2]/nrow(x), '%', range.exp[1], 100*range.exp[1]/nrow(x), '%')
        text.y <- ifelse(sum(colMin(x)>max(x)/2)/ncol(x) > 0.5, 0.1, 0.8)
        text(range.exp[2], max(x)*text.y, pos=text.pos, labels=exp.text)
        axis(1, c(1,nrow(x)), labels=c(1,nrow(x)))
        axis(2, at=axseq, labels=round(axseq,2), las=1)
    }
    invisible(class.range)
}

png("class_count_distributions.png", 2000, 1000)
class.range <- count.plot(class.ct)   # or count.plot(class.rp, TRUE) for RPKM distribs
dev.off()

# AS BARPLOT
class.range <- class.range[!is.na(class.range[,1]),]
x <- 100*rowMeans(class.range)/class.n
png("expr_pct_bars.png", 750, 500); par(cex=1.2, las=1)
b <- barplot(x, col=ncrna.dat[,2], ylim=c(0,100), main="Percent of Class Expressed, All Samples")
pm <- error.bars(b, x, 100*class.range/class.n)
dev.off()


### NCRNA HEATMAPS

hm.paths <- c("heatmaps_expr","heatmaps_expr_clust","heatmaps_trend","heatmaps_trend_clust")
for (i in 1:4) system(paste("mkdir",hm.paths[i]))

core.heatmap <- function(x, imgbase, title, rpkm, row.clust=TRUE, path=NULL) {
    x <- if (rpkm) { log2(x+1) } else { log10(x+1) }
    x <- x[rowSums(x)>0,,drop=FALSE]
    # EUCLIDEAN
    imgpath <- paste(ifelse(length(path)==0,hm.paths[1],path),"/",imgbase,".expr.png",sep="")
    heat.map(x, r.dend=row.clust, aspect=c(2,1,12), pal="Reds", main=title, imgname=imgpath)
    imgpath <- paste(ifelse(length(path)==0,hm.paths[2],path),"/",imgbase,".expr.clust.png",sep="")
    heat.map(x, r.dend=row.clust, c.dend=TRUE, aspect=c(2,1,12), pal="Reds", main=title, imgname=imgpath)
    # PEARSON
    z <- row.norm(x, divSD=TRUE)
    imgpath <- paste(ifelse(length(path)==0,hm.paths[3],path),"/",imgbase,".trend.png",sep="")
    heat.map(z, r.dend=row.clust, aspect=c(2,1,12), pal="CKY", metric="pearson", main=paste(title,", Z-Scores",sep=""), imgname=imgpath)
    imgpath <- paste(ifelse(length(path)==0,hm.paths[4],path),"/",imgbase,".trend.clust.png",sep="")
    heat.map(z, r.dend=row.clust, c.dend=TRUE, aspect=c(2,1,12), pal="CKY", metric="pearson", main=paste(title,", Z-Scores",sep=""), imgname=imgpath)
}

count.heatmap <- function(classes, rpkm=FALSE) {
    for (i in 1:length(classes)) {
        cname <- names(classes)[i]
        if (rpkm) {
            main <- paste(cname, "Log2 RPKMs")
            imgbase <- paste("heatmap", cname, "RPKM", sep=".")
        } else {
            main <- paste(cname, "Log10 Reads")
            imgbase <- paste("heatmap", cname, "count", sep=".")
        }
        core.heatmap(classes[[i]], imgbase, main, rpkm, TRUE)
    }
}

count.heatmap(class.rp, rpkm=TRUE)


## QC PLOTS

x <- do.call(rbind,class.ct)
rownames(x) <- paste(detable(class.n),rownames(x),sep="|")
x2 <- x[rowSums(x)>0,]
cm <- corr.mat(log2(x2+1), plot=TRUE, main=paste("Sample Corrs, Log2(Counts+1),",nrow(x2)," Expr. Genes"), filename="sample_corrs.png", pmar=c(10,10), pal="RYG", imgdim=c(700,600))

hc <- hclust(distance(log2(t(x2)+1)),"average")
png("sample_tree.png", 700, 600)
par(mar=c(4,4,4,10), cex=1.2, las=1)
plot(as.dendrogram(hc), horiz=TRUE, main=paste("Log2(Counts+1),",nrow(x2)," Expr. Genes, Pearson Dist."))
dev.off()

png("count_densities.png", 700, 600)
par(cex=1.2, las=1)
nzc <- apply(x, 2, function(y){ log10(y[y>0]) })
dhist(nzc, main="Log10 Nonzero Counts Per Gene", legend="topright")
dev.off()



## INITIAL SAVE (in case DE blows up)
save.image("analysis.RData")


### DIFFERENTIAL EXPRESSION PREP

if (!is.na(design)) {
    
    des <- read.delim(design, header=FALSE, as.is=TRUE)
    des <- mat.split(des[,2:ncol(des)], des[,1])
    C <- nrow(des$CONTRAST)
    et <- top <- sig <- new.list(des$CONTRAST[,1])
    
    ### DIFFERENTIAL EXPRESSION CALCULATION
        
    if ("REPLICATES" %in% names(des)) {
        
        outdir <- "EdgeR"   # use EdgeR, if replicates
        
        pre.group <- new.list(des$REPLICATES[,1])
        for (i in 1:length(pre.group)) pre.group[[i]] <- unlist(strsplit(des$REPLICATES[i,2],";"))
        pre.group <- breakout(pre.group, reverse=TRUE)
        group <- factor(pre.group[match(colnames(x),pre.group[,2]),1])
        if (any(is.na(group))) stop("Cannot mix replicated and non-replicated samples into one analysis, sorry!  Cannot proceed with DE testing.\n")
        
        x <- x[,order(group)]
        group <- group[order(group)]
        
        y <- DGEList(counts=x,group=group)
        y <- calcNormFactors(y)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        for (i in 1:C) {
            et[[i]] <- exactTest(y, pair=match(rev(des$CONTRAST[i,2:3]),levels(group)))
            et[[i]]$table <- cbind(et[[i]]$table, AdjPValue=p.adjust(et[[i]]$table[,3],"BH"), decideTestsDGE=decideTestsDGE(et[[i]], p=0.05, adjust="BH"))
            top[[i]] <- topTags(et[[i]], nrow(x))
            top[[i]]$table <- top[[i]]$table[,c(1:4,6,5)]
            sig[[i]] <- top[[i]]$table[top[[i]]$table$decideTestsDGE!=0,]
        }
        
        norm.ct.avg <- aggregate.cols(y$pseudo.counts, list(group), mean)
        norm.ct.sd <- aggregate.cols(y$pseudo.counts, list(group), sd)
        
        DE.blurb <- "Samples were replicated and analysis was conducted in EdgeR."
        DE.readme <- "This spreadsheet contains EdgeR significant-DE data for all contrasts (BH adjusted p-value <= 0.05)."
        DE.title <- "DE Genes at BH-adjusted p <= 0.05"
        
    } else {
        
        outdir <- "2FC"   # no replicates, do 2-fold-change
        y <- RPM.norm(x, M=mean(colSums(x)))  # RPMs, scaled to mean sample read count
        y2 <- y+min(nonzero(y))/2  # pseudocount-adjusted
        
        for (i in 1:C) {
            top[[i]] <- new.list("table")
            w <- y2[,match(unlist(des$CONTRAST[i,2:3]), colnames(x))]
            lfc <- log2(w[,1]/w[,2])
            top[[i]]$table <- data.frame(logFC=lfc, logCPM=log2(rowMeans(w)), Pvalue=1, AdjPValue=1, FDR=1, decideTestsDGE=sign(lfc))
            top[[i]]$table$decideTestsDGE[ abs(top[[i]]$table$logFC) < 1 ] <- 0
            sig[[i]] <- top[[i]]$table[ top[[i]]$table$decideTestsDGE!=0, ]
        }
        
        norm.ct.avg <- norm.ct.sd <- y2 # there are no means
        norm.ct.sd[norm.ct.sd>0] <- 0   # there are no SDs
        
        DE.blurb <- "Samples were not replicated, analysis is based on 2-fold-change cutoff."
        DE.readme <- "This spreadsheet contains DE data for all contrasts (absolute FC >= 2)."
        DE.title <- "DE Genes at abs(FC) >= 2"
    }

    system(paste("mkdir",outdir))

    
    # ALL-CONTRASTS SPREADSHEET
    readme <- matrix(c("This spreadsheet contains DE data for all contrasts.","",
               DE.blurb,"",
               "","",
               "Explanation of Columns:","",
               "","",
               "Fasta.Header","Fasta sequence header",
               "X MEAN","Scaled read count mean for X replicates",
               "Y MEAN","Scaled read count mean for Y replicates",
               "X STDEV","Scaled read count standard deviation for X replicates",
               "Y STDEV","Scaled read count standard deviation for Y replicates",
               "logFC","Log-fold-change in means, adjusted for count dispersion",
               "logCPM","Log-adjusted-counts-per-million (mean for both samples in comparison), adjusted for count dispersion",
               "PValue","Raw exact-test p-value",
               "AdjPValue","BH-adjusted p-value",
               "FDR","False discovery rate",
               "DE","Differential-expression flag: 1=up, -1=down, 0=not DE",
               "Symbol","Gene symbol (or other name)",
               "Gene.ID","Gene ID, (Ensembl, or NCBI Accession)",
               "Transcript.ID","Transcript ID (Ensembl, or NCBI GI)",
               "Chromosome","Chromosome Name",
               "Start","Transcript start",
               "End","Transcript end",
               "Strand","Transcript strand",
               "Biotype","Gene biotype",
               "Exons","Transcript exon count",
               "Length","Transcript length",
#               "GC","Transcript G/C percent",
#               "Host ID","Gene ID for ncRNA host gene",
#               "Host Symbol","Gene symbol for ncRNA host gene",
#               "Host Relation","Host-ncRNA relationship: can be ",
               "Status","Gene status",
               "Description","Gene description"
               ), ncol=2, byrow=TRUE)
    
    top2 <- lapply(1:C, function(i){
        x <- top[[i]]$table
        w <- match(rownames(x), rownames(norm.ct.avg))
        pair <- match(des$CONTRAST[i,2:3],colnames(norm.ct.avg))
        y <- cbind(x, norm.ct.avg[w,pair], norm.ct.sd[w,pair], gdat[match(rownames(x),gdat[,1]),])
        colnames(y)[6:10] <- c("DE",paste(colnames(y)[7:10], qw(MEAN,MEAN,SD,SD), sep="."))
        y[,c(11,7:10,1:6,13:19,12,20:ncol(y))]
    })
    names(top2) <- sub("/","vs",names(top))
    for (s in c("/","\\\\","\\?","\\*",":","\\]","\\[")) names(top2) <- gsub(s,"_",names(top2))
    outfile <- paste0(outdir,"/All_Contrasts.xls")
    system(paste("rm -f",outfile))  # guarantee file is removed prior to writing (basically, hedge against file being open in Excel)
    WriteXLS2(c(README=list(readme),top2), outfile, BoldHeaderRow=TRUE, FreezeRow=1)
    
    if (sum(sapply(sig,nrow))==0) {
        
        IM("No differentially expressed genes!\n")
        
    } else {

        # SIG-GENES SPREADSHEET
        readme2 <- readme
        readme2[1,1] <- DE.readme
        sig2 <- lapply(top2, function(x){ y=x[order(x[,6],decreasing=TRUE),]; y[y$DE!=0,] })
        WriteXLS2(c(README=list(readme2),sig2), paste0(outdir,"/DE_genes.xls"), BoldHeaderRow=TRUE, FreezeRow=1)        
        
        # DE GENES BARPLOT
        png(paste0(outdir,"/DE_genes.png"), max(c(400,200+100*length(sig))), 500)
        par(cex=1.2, las=1)
        sigdir <- rbind(
            POS=NAfy(sapply(top, function(x){ sum(x$table$decideTestsDGE==1) }), include=0),
            NEG=-NAfy(sapply(top, function(x){ sum(x$table$decideTestsDGE==-1) }), include=0)
            )
        barplot(sigdir[1,,drop=FALSE], col="#FFAAAA", ylim=range(c(0,sigdir), na.rm=TRUE), main=DE.title)
        barplot(sigdir[2,], col="#AAAAFF", add=TRUE, names=NA)
        abline(h=0)
        dev.off()
        
        # MA PLOTS
        idims <- MA.scatter.plotdim(length(top))
        ylim <- range(real.only(unlist(slice.list(slice.list(top,"table"),1,"cols"))))
        xlim <- range(real.only(unlist(slice.list(slice.list(top,"table"),2,"cols"))))
        
        png(paste0(outdir,"/MA_plots.png"), idims[4], idims[3])
        par(mfrow=idims[1:2], family="mono", cex=1.3, las=1)
        for (i in 1:C) {
            z <- top[[i]]$table
            plot(z[,2:1], col=8, xlim=xlim, ylim=ylim, xlab="LogCPM", ylab="LogFC", main=paste(names(top)[i],"MA Plot"))
            abline(h=0)
            up <- which(z$decideTestsDGE==1)
            dn <- which(z$decideTestsDGE==-1)
            points(z[up,2:1,drop=FALSE], col=2)
            points(z[dn,2:1,drop=FALSE], col=3)
            maxdig <- max(nchar(c(sum(!is.infinite(z[,2])),length(up),length(dn))))
            ltext <- apply(cbind(c(sum(!is.infinite(z[,2])),length(up),length(dn)), qw(Expressed,SigUp,SigDown)), 1, function(s){ sprintf(paste0("%",maxdig,"s %s"),s[1],s[2]) })
            legend("topright", bty="n", pch=1, col=1:3, legend=ltext)
        }
        dev.off()
        
        # DE HEATMAP
        all.DE <- suniq(unlist(lapply(sig, rownames)))
        all.rp <- do.call(rbind, class.rp)
        DE.rp <- all.rp[match(sub("^.*?\\|","",all.DE), rownames(all.rp)),,drop=FALSE]
        rownames(DE.rp) <- all.DE
        core.heatmap(DE.rp, "heatmap.All-DE.RPKM", "All DE Genes, Log2 RPKMs", rpkm=TRUE, row.clust=FALSE, path=outdir)
    }
}


### FINAL SAVE, EXIT
save.image("analysis.RData")
quit()



