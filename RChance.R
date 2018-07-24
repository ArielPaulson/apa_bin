#!/usr/bin/env Rscript
source("/n/projects/apa/R/apa_tools.R")

## Call "/home/apa/local/bin/RChance.R <outprefix> <sample1> <sample2> ... <sampleN>
##      where each <sample> has format "[IP|input]=alias,/path/to/sample.RChipQC.binMeans.RData"





## Setup

## Get command line
ca <- commandArgs(trailing=TRUE)
outprefix <- sub("\\.?RChance$","",ca[1])  # if 'RChance' included in outprefix, remove it, because it gets added on the line below
RDatas <- sort(ca[2:length(ca)])   ## expecting "IP=alias,/path/to/.RChipQC.binMeans.RData" or "input=alias,/path/to/.RChipQC.binMeans.RData" strings

outprefix <- ifelse(falsify(file.info(outprefix)$isdir), paste0(outprefix,"/RChance"), paste0(outprefix,".RChance"))
message(paste("\nOUTPREFIX:",outprefix))

## Was '--clobber' indicated?
w.clob <- which(RDatas=="--clobber")
if (length(w.clob)>0) {
    clobber <- TRUE
    RDatas <- RDatas[-w.clob]
} else {
    clobber <- FALSE
}
message(paste("CLOBBER:",clobber))

## Checkpoint outputs
emdata <- paste0(outprefix,".emergency.RData")  # first session save object; exists only transiently, unless early failure in object creation
initdata <- paste0(outprefix,".initial.RData")  # second session save object; exists only until script finishes, then replaced by 'finaldata'
finaldata <- paste0(outprefix,".RData")         # session save object

## Quit if run appears complete, and not clobbering
if (!clobber & file.exists(finaldata)) {
    message(paste0("Final save file '",finaldata,"' already exists!  Stopping..."))
    quit()
}    

## Are we restarting from a saved session?

if (!clobber & file.exists(initdata)) {
    
    ## Starting from existing mature RData object
    
    this.RDatas <- RDatas
    load(initdata)  # overwrites 'RDatas'
    
    if (all(RDatas==this.RDatas)) {
        ## ok to proceed
        message(paste("Restarting from saved session '",initdata,"'",sep=""))
    } else {
        ## sample list has changed: cannot proceed from save point
        stop("Sample list has changed: cannot proceed from saved session '",initdata,"'!",sep="")
    }
    rm(this.RDatas)  # don't want it getting saved to later objects!
    
    from.scratch <- FALSE
    from.emerg <- FALSE
    from.initial <- TRUE
    
} else if (!clobber & file.exists(emdata)) {
     
    ## Starting from existing emergency RData object
    ## GUARANTEED TO FAIL, unless you have _just_ edited this script to fix the error
    
    this.RDatas <- RDatas
    load(emdata)  # overwrites 'RDatas'
    
    if (all(RDatas==this.RDatas)) {
        ## ok to proceed
        message(paste("Restarting from emergency session '",emdata,"'",sep=""))
    } else {
        ## sample list has changed: cannot proceed from save point
        stop("Sample list has changed: cannot proceed from emergency session '",emdata,"'!",sep="")
    }
    rm(this.RDatas)  # don't want it getting saved to later objects!
    
    from.scratch <- FALSE
    from.emerg <- TRUE
    from.initial <- FALSE
    
} else {
    
    ## Normal start
    
    ## Parse sample strings
    names(RDatas) <- extract.names.from.paths(RDatas)
    N <- length(RDatas)
    RDatas2 <- data.frame(Type=rep("",N),Alias=rep("",N),File=rep("",N))
    for (i in 1:N) {
        x <- unlist(strsplit(RDatas[i],"[=,]"))  # ideally, length-3 vector: IP|input, alias, RData
        if (length(x)==2) x <- c(x[1], sub(".RChipQC.binMeans.RData$","",sub(".*/","",x[2])), x[2])  # or length-2 vector... IP|input, RData: extract alias from RData name (must be "/path/to/alias.RChipQC.binMeans.RData")
        if (!file.exists(x[3])) stop(paste0("RData file '",x[3],"' does not exist!  Stopping.\n"))
        real.file <- system(paste("readlink -f",x[3]),intern=TRUE)
        if (file.exists(real.file)) x[3] <- real.file
        RDatas2[i,] <- x
    }
    
    ## Fix IP/input cases; test for bad samples
    RDatas2[grepl("^IP$",RDatas2$Type,ignore.case=TRUE),1] <- "IP"
    RDatas2[grepl("^input$",RDatas2$Type,ignore.case=TRUE),1] <- "input"
    NIP <- sum(RDatas2$Type=="IP")  # multi-IP check
    NIN <- sum(RDatas2$Type=="input")  # multi-input check
    bad <- which(!(RDatas2$Type %in% c("IP","input")))  # unassignable to IP or input?
    RDatas2 <- RDatas2[c(which(RDatas2$Type=="IP"),which(RDatas2$Type=="input")),]
    
    ## Print sample table to screen
    message("\nSample table:")
    message(paste(colnames(RDatas2),collapse="\t"))
    for (i in 1:N) message(paste(RDatas2[i,],collapse="\t"))
    
    ## Die if bad samples
    if (length(bad)>0) {
        message("\nThe following sample strings could not be assigned to IP or input:")
        for (i in bad) message(RDatas[i])
        stop("Cannot proceed.")
    }
    
    from.scratch <- TRUE
    from.emerg <- FALSE
    from.initial <- FALSE
    
}





#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#####################################################################   FUNCTIONS   #####################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################





set.palettes <- function (ndat) {
    
    ## 'ndat' is the output of chance.normalize.samples()
    
    LIP <- length(ndat$IP$norm)
    LIN <- length(ndat$input$norm)
    
    IP.pal <- c("red","gold","darkorange","brown2")  # up to 4 non-interpolated colors (IPs)
    IPC.pal <- "darkred"  # consensus-IP color
    inp.pal <- c("green3","cyan","dodgerblue","slateblue2")  # up to 4 non-interpolated colors (inputs)
    inpC.pal <- "blue"    # consensus-input color
    
    if (LIP>NIP) {
        ## multiple IPs, so consensus was generated
        if (LIP-1<=length(IP.pal)) {
            IP.pal <- c(IPC.pal, IP.pal[1:(LIP-1)])
        } else {
            IP.pal <- c(IPC.pal, palettizer(IP.pal, LIP-1))
        }
        IP.lty <- c(2, rep(1, LIP-1))
    } else {
        ## single IP, so consensus was NOT generated
        if (LIP<=length(IP.pal)) {
            IP.pal <- IP.pal[1:LIP]
        } else {
            IP.pal <- palettizer(IP.pal, LIP)
        }
        IP.lty <- rep(1, LIP)
    }

    if (LIN>NIN) {
        ## multiple inputs, so consensus was generated
        if (LIN-1<=length(inp.pal)) {
            inp.pal <- c(inpC.pal, inp.pal[1:(LIN-1)])
        } else {
            inp.pal <- c(inpC.pal, palettizer(inp.pal, LIN-1))
        }
        inp.lty <- c(2, rep(1, LIN-1))
    } else {
        ## single input, so consensus was NOT generated
        if (LIN<=length(inp.pal)) {
            inp.pal <- inp.pal[1:LIN]
        } else {
            inp.pal <- palettizer(inp.pal, LIN)
        }
        inp.lty <- rep(1, LIN)
    }
    
    list(IP.pal=IP.pal, IPC.pal=IPC.pal, inp.pal=inp.pal, inpC.pal=inpC.pal, IP.lty=IP.lty, inp.lty=inp.lty)
    
}



chance.normalize.samples <- function(IP.list, input.list) {
    
    ## Normalization and cumulative-summing of samples.
    ## This generates the curves for the IP-input enrichment plots.
    ## 'IP.list' and 'input.list' are lists of numeric vectors of 1k bin means along genome; all vectors must have same length.  (These can be found in RChipQC.binMeans.RData objects)
    ##   Alternatively, they can be matrices with the vectors as columns
    
    consensus <- TRUE
    if (is(IP.list,"matrix")) IP.list <- as.list(as.data.frame(IP.list))
    if (is(input.list,"matrix")) input.list <- as.list(as.data.frame(input.list))
    
    LIP <- length(IP.list)
    LIN <- length(input.list)
    bins <- length(IP.list[[1]])
    IM(LIP,"IP samples and",LIN,"input samples, using",bins,"bins")
    
    IP.csp <- IP.list
    for (i in 1:LIP) {
        IM("Normalizing",names(IP.list)[i])
        names(IP.list[[i]]) <- 1:bins
        IP.list[[i]] <- 1E6*sort(IP.list[[i]])/sum(as.numeric(IP.list[[i]]))  # sorted RPMs
        IP.csp[[i]] <- cumsum(as.numeric(IP.list[[i]]))
        IP.csp[[i]] <- IP.csp[[i]]/max(IP.csp[[i]])
    }
    if (consensus) {
        if (length(IP.list)==1) {
            # no consensus for single samples
        } else {
            IM("Generating IP consensus")
            IP.cons <- colMeans(do.call(rbind,IP.list))
            IP.cons.csp <- cumsum(as.numeric(IP.cons))
            IP.cons.csp <- IP.cons.csp/max(IP.cons.csp)
            IP.list <- c(list(IP.consensus=IP.cons), IP.list)
            IP.csp <- c(list(IP.consensus=IP.cons.csp), IP.csp)
        }
    }
    
    input.csp <- input.list
    for (i in 1:LIN) {
        IM("Normalizing",names(input.list)[i])
        names(input.list[[i]]) <- 1:bins
        input.list[[i]] <- 1E6*sort(input.list[[i]])/sum(as.numeric(input.list[[i]]))  # sorted RPMs
        input.csp[[i]] <- cumsum(as.numeric(input.list[[i]]))
        input.csp[[i]] <- input.csp[[i]]/max(input.csp[[i]])
    }
    if (consensus) {
        if (length(input.list)==1) {
            # no consensus for single samples
        } else {
            IM("Generating input consensus")
            input.cons <- colMeans(do.call(rbind,input.list))
            input.cons.csp <- cumsum(as.numeric(input.cons))
            input.cons.csp <- input.cons.csp/max(input.cons.csp)
            input.list <- c(list(input.consensus=input.cons), input.list)
            input.csp <- c(list(input.consensus=input.cons.csp), input.csp)
        }
    }
    
    list(IP=list(raw=IP.list, norm=IP.csp), input=list(raw=input.list, norm=input.csp))
}



chance.enrichment.curve.plot <- function(ndat, palettes, consensus=FALSE, cex=1, lwd=1) {
    
    ## plot a Chance-style IP-vs-input cumulative enrichment plot
    ## 'ndat' is the output of chance.normalize.samples()
    ## 'palettes' is the output of set.palettes()
    ## 'consensus=TRUE' plots IP and input consensus lines
    
    IM("Calculating Enrichment")
    cutoffs <- sapply(ndat$IP$norm, function(y) nameless(which.max(abs(y-ndat$input$norm[[1]]))) )
    bins <- length(ndat$IP$norm[[1]])
    LIP <- length(ndat$IP$norm)
    LIN <- length(ndat$input$norm)
    nmax <- max(nchar(c(names(ndat$IP$norm),names(ndat$input$norm))))
    diagonal <- seq(0,1,length.out=length(ndat$input$norm[[1]]))
    sum.diagonal <- sum(diagonal)
    
    IPI <- 1:LIP; names(IPI) <- names(ndat$IP$norm)
    INI <- 1:LIN; names(INI) <- names(ndat$input$norm)
    
    above.below <- function(x,thresh) sum(x[thresh:length(x)])/sum(x)
    total.enrich <- function(x) (sum.diagonal-sum(x))/sum.diagonal

    if (LIN>NIN) {
        cutoffs.inp <- c(NA, sapply(ndat$input$norm[2:LIN], function(y) nameless(which.max(abs(y-ndat$input$norm[[1]]))) ))
    } else {
        cutoffs.inp <- rep(NA, LIN)
    }
    names(cutoffs.inp) <- names(ndat$input$norm)
    
    L2=list(
        bins=bins,
        enrich.cutoff.global=cutoffs[1],
        enrich.cutoff.local=c(cutoffs,cutoffs.inp),
        enrich.bins.pct.global=(bins-cutoffs[1])/bins,
        enrich.bins.pct.local=(bins-c(cutoffs,cutoffs.inp))/bins,
        Input.enrich.signal.pct.global=sapply(INI, function(i) above.below(ndat$input$norm[[i]],cutoffs[1])),
        Input.enrich.signal.pct.local=sapply(INI, function(i) above.below(ndat$input$norm[[i]],cutoffs[i])),
        IP.enrich.signal.pct.global=sapply(IPI, function(i) above.below(ndat$IP$norm[[i]],cutoffs[1])),
        IP.enrich.signal.pct.local=sapply(IPI, function(i) above.below(ndat$IP$norm[[i]],cutoffs[i])),
        Input.pct.total.enrich=sapply(ndat$input$norm, total.enrich),
        IP.pct.total.enrich=sapply(ndat$IP$norm, total.enrich),
        IP.pct.enrich.over.input=c()
    )
    mean.Input.enrich <- ifelse(LIN>NIN, L2$Input.pct.total.enrich[1], mean(L2$Input.pct.total.enrich))
    L2$IP.pct.enrich.over.input <- L2$IP.pct.total.enrich/mean.Input.enrich-1
    
#    legend.text <- c(
#        sapply(1:LIP, function(i){ x=L2$IP.enrich.signal.pct[i];    gsub("_(,|\\))","%\\1",sprintf(paste0("%-",nmax,"s (%4.1f_, %4.1f_, %4.1f_)"),names(x),100*L2$IP.enrich.signal.pct[i],100*L2$IP.pct.total.enrich[i],100*L2$IP.pct.enrich.over.input[i])) }),
#        sapply(1:LIN, function(i){ x=L2$Input.enrich.signal.pct[i]; gsub("_(,|\\))","%\\1",sprintf(paste0("%-",nmax,"s (%4.1f_, %4.1f_)"),names(x),100*L2$Input.enrich.signal.pct[i],100*L2$Input.pct.total.enrich[i])) })
#        ## Input enrichment deviation from mean calculated, but not reported
#    )
#    legend.text <- c(legend.text, paste0("Signal/Bkg cutoff (",sprintf("%0.1f",100*L2$cutoff.pct.bins),"%)"))
    
    if (consensus) {
        use.IP <- 1:LIP
        use.inp <- 1:LIN
    } else {
        if (LIP>NIP) {
            use.IP <- 2:LIP
        } else {
            use.IP <- 1:LIP
        }
        if (LIN>NIN) {
            use.inp <- 2:LIN
        } else {
            use.inp <- 1:LIN
        }
    }

    nmax <- max(nchar(names(c(L2$IP.enrich.signal.pct.global[use.IP],L2$Input.enrich.signal.pct.global[use.inp]))))
    legend.text <- c(
        sapply(use.IP, function(i){
            x <- L2$IP.enrich.signal.pct.global[i]
            paste0(
                sprintf(paste0("%-",nmax,"s"),names(x)),
                gsub("_","%",sprintf(" (%0.1f_, %0.1f_)",100*L2$IP.pct.total.enrich[i],100*L2$IP.pct.enrich.over.input[i]))
            )
        }),
        sapply(use.inp, function(i){
            x <- L2$Input.enrich.signal.pct.global[i]
            paste0(
                sprintf(paste0("%-",nmax,"s"),names(x)),
                gsub("_","%",sprintf(" (%0.1f_)",100*L2$Input.pct.total.enrich[i]))
            )
        })
        ## Input enrichment deviation from mean calculated, but not given in plot
    )
    legend.text <- c(legend.text, paste0("Enrichment Cutoff (",sprintf("%0.1f",100*(1-L2$enrich.bins.pct.global)),")"))
    
    par(family="mono", las=1, cex=cex)
    null.plot(xlim=c(0,bins), ylim=c(0,1), xlab="Percentage of Bins", ylab="Percent of Total Signal", main="Cumulative Enrichment in Each Sample")
    for (i in rev(use.IP)) lines(1:bins, ndat$IP$norm[[i]], col=palettes$IP.pal[i], lty=palettes$IP.lty[i], lwd=2)   # reverse order of everything so consensus prints last
    for (i in rev(use.inp)) lines(1:bins, ndat$input$norm[[i]], col=palettes$inp.pal[i], lty=palettes$inp.lty[i], lwd=2)
    abline(0, 1/length(ndat$input$norm[[1]]), col=8)
    abline(v=cutoffs[1], col="grey25", lwd=2, lty=3)
#    text(cutoffs[1], 0, col="grey25", pos=2, label=cutoffs[1])
    axis(1, at=seq(0,bins,length.out=11), labels=seq(0,100,10))
    axis(2, at=seq(0,1,0.1), labels=seq(0,100,10))
#    par(family="mono", cex=0.9)
    mtext("Sample (Total Enrichment over Null, IP Enrichment over Input)", 3, 0, cex=0.9, col="grey45")
    legend("topleft", bty="n", lwd=2, col=c(palettes$IP.pal[use.IP],palettes$inp.pal[use.inp],1), lty=c(palettes$IP.lty[use.IP],palettes$inp.lty[use.inp],3), legend=legend.text)
    
    invisible(L2)
}



chance.enrichment.overlap.plot <- function(ndat, edat, palettes, consensus=FALSE, cex=1, label.cex=1.5, pmar=c(10,10)) {
    
    ## plot a Chance-style multi-IP differential enrichment plot
    ## 'ndat' is the output of chance.normalize.samples()
    ## 'edat' is the output of chance.enrichment.curve.plot()
    ## 'palettes' is the output of set.palettes()
    ## 'consensus=TRUE' plots IP and input consensus lines
    ##
    ## UNDER CONSTRUCTION - runs, but have not recapitulated original Chance behavior
    
    LIP <- length(ndat$IP$norm)
    LIN <- length(ndat$input$norm)
    use.IP <- 1:LIP
    use.inp <- 1:LIN
    if (!consensus) {
        if (LIP>NIP) use.IP <- use.IP[2:LIP]
        if (LIN>NIN) use.inp <- use.inp[2:LIN]
    }
    LIP <- length(use.IP)
    LIN <- length(use.inp)
    
    LTOT <- LIP+LIN
    xtemp <- c(ndat$IP$norm[use.IP], ndat$input$norm[use.inp])
    allnames <- names(xtemp)
    bins <- edat$bins
    cutoff <- edat$enrich.cutoff.global
    ebins.N <- length(cutoff:bins)
    
    pair.same <- pair.diff <- new.list(qw(pooled,possible,total), elem=matrix(0, LTOT, LTOT, FALSE, list(allnames,allnames)))
    
    label.col <- pair.same[[1]]
    label.col[label.col==0] <- 1

    for (i in 1:LTOT) {
        ebins.i <- as.numeric(names(xtemp[[i]])[cutoff:bins])                          # genomic bins which are enriched in sample 1 (bin signal >= cutoff)
        for (j in 1:LTOT) {
            ebins.j <- as.numeric(names(xtemp[[j]])[cutoff:bins])                      # genomic bins which are enriched in sample 2 (bin signal >= cutoff)
            pooled.DE <- length(unique(c(ebins.i,ebins.j)))                            # unique enriched genomic bins from both samples
            common.DE <- length(intersect(ebins.i,ebins.j))                            # enriched genomic bins shared by both samples
#            IM(i,j,ebins.N,length(ebins.i),length(ebins.j),pooled.DE,common.DE,pooled.DE-common.DE)
            pair.same$pooled[i,j] <- round(100*common.DE/pooled.DE,2)                  # enriched genomic bins shared by both samples, as percent unique enriched bins in both samples
            pair.diff$pooled[i,j] <- round(100*(pooled.DE-common.DE)/pooled.DE,2)      # enriched genomic bins NOT shared by both samples, as percent unique enriched bins in both samples
            pair.same$possible[i,j] <- round(100*common.DE/ebins.N,2)                  # enriched genomic bins shared by both samples, as percent all enriched bins
            pair.diff$possible[i,j] <- round(100*(pooled.DE-common.DE)/(2*ebins.N),2)  # enriched genomic bins NOT shared by both samples, as percent all enriched bins
            pair.same$total[i,j] <- round(100*common.DE/bins,2)                        # enriched genomic bins shared by both samples, as percent all bins
            pair.diff$total[i,j] <- round(100*(pooled.DE-common.DE)/bins,2)            # enriched genomic bins NOT shared by both samples, as percent all bins
        }
    }
    diag(pair.same$pooled) <- 100  # inputs will have 0; set to 100 for consistency's sake
    
    par(mfcol=c(2,3), las=2)
    mipu(pair.same$pooled, pal=c(5,6), show.scale=FALSE, col.limits=c(0,100), main="Shared Enriched Bins, % Pooled", pmar=pmar, cex=cex, label.text=apply(pair.same$pooled,2,function(x) paste0(x,"%")), label.col=label.col, label.cex=cex)
    mipu(pair.diff$pooled, pal=c(5,6), show.scale=FALSE, col.limits=c(0,100), main="Disjoint Enriched Bins, % Pooled", pmar=pmar, cex=cex, label.text=apply(pair.diff$pooled,2,function(x) paste0(x,"%")), label.col=label.col, label.cex=cex)
    mipu(pair.same$possible, pal=c(5,6), show.scale=FALSE, col.limits=c(0,100), main="Shared Enriched Bins, % Possible", pmar=pmar, cex=cex, label.text=apply(pair.same$possible,2,function(x) paste0(x,"%")), label.col=label.col, label.cex=cex)
    mipu(pair.diff$possible, pal=c(5,6), show.scale=FALSE, col.limits=c(0,100), main="Disjoint Enriched Bins, % Possible", pmar=pmar, cex=cex, label.text=apply(pair.diff$possible,2,function(x) paste0(x,"%")), label.col=label.col, label.cex=cex)
    mipu(pair.same$total, pal=c(5,6), show.scale=FALSE, col.limits=c(0,100), main="Shared Enriched Bins, % Total", pmar=pmar, cex=cex, label.text=apply(pair.same$total,2,function(x) paste0(x,"%")), label.col=label.col, label.cex=cex)
    mipu(pair.diff$total, pal=c(5,6), show.scale=FALSE, col.limits=c(0,100), main="Disjoint Enriched Bins, % Total", pmar=pmar, cex=cex, label.text=apply(pair.diff$total,2,function(x) paste0(x,"%")), label.col=label.col, label.cex=cex)
    
    invisible(list(pair.same=pair.same, pair.diff=pair.diff))
}



shared.peak.matrix <- function(x) {
    
    ## alternative to chance.enrichment.overlap.plot
    ## UNDER CONSTRUCTION
    
    
}



multi.CCF.plot <- function(all.CCF, all.stats, palettes, cex=1, lwd=1) {
    
    ## 'all.CCF' is a matrix of c-saw CCF vectors, made from the incoming RData objects, ONE ROW PER SAMPLE / EXPECTS 1001 COLUMNS
    ## 'all.stats' is a matrix of RChipQC.R stats, made from the incoming RData objects, ONE ROW PER SAMPLE
    ## 'palettes' is the output of set.palettes()
    
    samp.ord <- c(which(RDatas2$Type=="IP"),which(RDatas2$Type=="input"))
    RL <- all.stats$int[rownames(all.stats$int)=="Read Length",]  # read lengths per sample
    
    N <- nrow(all.CCF)
    all.CCF2 <- 100*all.CCF/apply(all.CCF,1,max)
    all.CCFL <- t(sapply(1:N, function(i) lowess(all.CCF2[i,], f=RL[i]/2002, delta=0)$y ))  # f = read length / 2, as % ncol all.CCF
    max.xcor <- apply(all.CCF, 1, max)
    max.shift <- apply(all.CCF, 1, which.max) - 1  # shift sizes start at 0
    
    all.pal <- rep("", N)
    LIP <- length(ndat$IP$norm)
    LIN <- length(ndat$input$norm)
    if (LIP>NIP) {
        all.pal[which(RDatas2$Type=="IP")] <- palettes$IP.pal[2:length(palettes$IP.pal)]
    } else {
        all.pal[which(RDatas2$Type=="IP")] <- palettes$IP.pal
    }
    if (LIN>NIN) {
        all.pal[which(RDatas2$Type=="input")] <- palettes$inp.pal[2:length(palettes$inp.pal)]
    } else {
        all.pal[which(RDatas2$Type=="input")] <- palettes$inp.pal
    }
    
    ## set legend text & abline(s)
    ## if all same RL, one abline (black), one legend line
    ## if diff RLs, one abline per sample (colored), legend text per sample = mCCF + RL
    ## stats per sample: read len (if variable), max x, mCCF, 
    if (length(unique(RL))==1) {
        ## hoped-for situation
        sprintf.format <- paste0("%-",max(nchar(RDatas2$Alias)),"s: %",nchar(max(max.shift)),"g: %0.4f")
        ltext <- sapply(samp.ord, function(i) sprintf(sprintf.format, RDatas2[i,"Alias"], max.shift[i], max.xcor[i]))
        subtitle <- "Sample : Best Shift : Best X-Cor"
        ltext <- c(ltext, "Read Length")
        ablines <- RL[1]
        abcols <- 1
    } else {
        ## apparently, bad experiment design or sample pairing
        sprintf.format <- paste0("%-",max(nchar(RDatas2$Alias)),"s: %",nchar(max(RL)),"i: %",nchar(max(max.shift)),"g: %0.4f")
        ltext <- sapply(samp.ord, function(i) sprintf(sprintf.format, RDatas2[i,"Alias"], RL[i], max.shift[i], max.xcor[i]))
        subtitle <- "Sample : Read Length : Best Shift : Best X-Cor"
        ltext <- c(ltext, "Read Lengths")
        ablines <- RL
        abcols <- all.pal
    }   
    Lab <- length(ablines)
    
    par(family="mono", xaxs="i", las=1, cex=cex)
    null.plot(xlim=c(-1,1000), ylim=c(0,100), xlab="X-Corr Shift Size (bp)", ylab="X-Corr as Pct of Max", main="Sample Cross-Correlations")
    axis(1, at=seq(0,1000,100))
    axis(2, at=seq(0,100,10))
    abline(v=ablines+1, col=abcols[1:Lab], lty=3)
    text(ablines, seq(0,10*N/3,length.out=Lab+1)[1:Lab], labels=paste0(ablines,"bp"), pos=4)
    for (i in 1:N) {
        lines(1:ncol(all.CCF2), all.CCF2[i,], col=all.pal[i])
        points(RL[i]+1, all.CCF2[i,RL[i]+1], col=all.pal[i])
        wmi <- which.max(all.CCF2[i,])
        if (wmi!=RL[i]+1) points(wmi, all.CCF2[i,wmi], col=all.pal[i], pch=2)
    }
    mtext(subtitle, 3, 0, cex=0.9, col="grey45")
    legend("topright", bty="n", col=c(all.pal[samp.ord],1), lty=c(rep(1,N),3), legend=ltext)
    
    all.stats.plus <- all.stats
    bxv <- which(rownames(all.stats$int)=="Best X-Cor Shift")
    all.stats.plus$int <- all.stats$int[c(1:bxv,bxv:nrow(all.stats$int)),]
    rownames(all.stats.plus$int)[bxv+1] <- "Second-Best X-Cor Shift"
    truemax.shift <- sapply(1:nrow(all.CCF), function(i) which.max(all.CCF[i,(2*RL[i]+1):ncol(all.CCF)])+(2*RL[i]) ) - 1
    truemax.shift[truemax.shift<=(2*RL[i]+11)] <- max.shift[i]  # if new shift size = threshold(+10bp), then don't bother
    all.stats.plus$int[bxv+1,] <- truemax.shift
    colnames(all.stats.plus$int) <- RDatas2$Alias
    invisible(all.stats.plus)
}



write.stats <- function(edat, all.stats) {
    
    ## basically, integrate 'edat' into 'all.stats' and write single table
    ## 'edat' is the output of chance.enrichment.curve.plot()
    ## 'all.stats' is a matrix of RChipQC.R stats, made from the incoming RData objects, ONE ROW PER SAMPLE

    s <- rownameless(do.call(cbind, lapply(all.stats, function(x) as.data.frame(t(x)) )))
    colnames(s) <- sub("^(char|int|float)\\.","",colnames(s))
    orig.stats.ord2 <- c(colnames(s)[1:8],orig.stats.ord[4:length(orig.stats.ord)])
    s <- s[,match(orig.stats.ord2,colnames(s))]
    is.IP <- which(s$Type=="IP")
    is.inp <- which(s$Type=="input")
    
    s2 <- rownameless(cbind(
        s,
        TotalBins=edat$bins,
        EnrichedSignalPct=NA,
        EnrichedSignalPctGlobal=NA,
        TotalEnrichmentPctOverNull=NA,
        IPEnrichmentPctOverInput=NA,
        EnrichedBinCutoff=NA,
        EnrichedBinPct=NA,
        EnrichedBinCutoffGlobal=NA,
        EnrichedBinPctGlobal=NA
    ))
    
    s2$EnrichedBinCutoff <- edat$enrich.cutoff.local[match(s$Alias,names(edat$enrich.cutoff.local))]
    s2$EnrichedBinPct <- edat$enrich.bins.pct.local[match(s$Alias,names(edat$enrich.bins.pct.local))]
    s2$EnrichedSignalPct[is.IP] <- edat$IP.enrich.signal.pct.local[match(s$Alias[is.IP],names(edat$IP.enrich.signal.pct.local))]
    s2$EnrichedSignalPct[is.inp] <- edat$Input.enrich.signal.pct.local[match(s$Alias[is.inp],names(edat$Input.enrich.signal.pct.local))]
    s2$EnrichedSignalPctGlobal[is.IP] <- edat$IP.enrich.signal.pct.global[match(s$Alias[is.IP],names(edat$IP.enrich.signal.pct.global))]
    s2$EnrichedSignalPctGlobal[is.inp] <- edat$Input.enrich.signal.pct.global[match(s$Alias[is.inp],names(edat$Input.enrich.signal.pct.global))]
    s2$TotalEnrichmentPctOverNull[is.IP] <- edat$IP.pct.total.enrich[match(s$Alias[is.IP],names(edat$IP.pct.total.enrich))]
    s2$TotalEnrichmentPctOverNull[is.inp] <- edat$Input.pct.total.enrich[match(s$Alias[is.inp],names(edat$Input.pct.total.enrich))]
    s2$IPEnrichmentPctOverInput[is.IP] <- edat$IP.pct.enrich.over.input[match(s$Alias[is.IP],names(edat$IP.pct.enrich.over.input))]
    s2$EnrichedBinCutoffGlobal[is.IP] <- edat$enrich.cutoff.global
    s2$EnrichedBinPctGlobal[is.IP] <- edat$enrich.bins.pct.global
    for (i in (ncol(s)+1):ncol(s2)) if (any(grepl("\\.",s2[,i]))) s2[,i] <- round(100*s2[,i],2)  # only round float columns!!!
    
    WriteXLS2(list(gsub(" ","",t(rownamed.matrix(s2,"character")))), paste0(outprefix,".stats.xls"), FreezeRow=1, FreezeCol=1, BoldHeaderRow=TRUE, row.names=TRUE)
    write.table2(s2, paste0(outprefix,".stats.txt"))
    
    ## write multi-IP matrix?  multi-IP stats?  multi-input stats?
}

## Stats from enrichment plot 'L2' object (globally, 'edat')
#    L2=list(
#        bins=bins,
#        enrich.cutoff.global=cutoffs[1],
#        enrich.cutoff.local=cutoffs,
#        enrich.bins.pct.global=cutoffs[1]/bins,
#        enrich.bins.pct.local=cutoffs/bins,
#        Input.enrich.signal.pct.global=sapply(INI, function(i) above.below(ndat$input$norm[[i]],cutoffs[1])),
#        Input.enrich.signal.pct.local=sapply(INI, function(i) above.below(ndat$input$norm[[i]],cutoffs[i])),
#        IP.enrich.signal.pct.global=sapply(IPI, function(i) above.below(ndat$IP$norm[[i]],cutoffs[1])),
#        IP.enrich.signal.pct.local=sapply(IPI, function(i) above.below(ndat$IP$norm[[i]],cutoffs[i])),
#        Input.pct.total.enrich=sapply(ndat$input$norm, total.enrich),
#        IP.pct.total.enrich=sapply(ndat$IP$norm, total.enrich),
#        IP.pct.enrich.over.input=c()
#    )

## Stats from RChipQC.R 'stats.table' object (globally, 'all.stats')
#    c("Read Length", readlen),
#    c("Best X-Cor Shift", best.shift),
#    c("Best X-Cor Value", round(best.xcor,4)),
#    c("Alignments", aligns),
#    c("Genomic Coverage", sci.cov(as.numeric(aligns)*as.numeric(readlen)/genome) ),
#    c("Genome Bp", genome),
#    c("Nonzero Pos", bpCov),
#    c("Nonzero Pos %", sci.pct(bpCov/genome) ),
#    c("Total Bp Aligned", bpSeq),
#    c("Pos in Top 1% of Depths", top1pct.pos),
#    c("Bp in Top 1% Pos", top1pct.seq),
#    c("% Bp Aligned in Top 1% Pos",sci.pct(top1pct.seq/bpSeq) ),
#    c("Unique Alignment Starts", upos),
#    c("Unique Alignment %", sci.pct(upos/aligns) ),
#    c("Min Height of Reportable Island", minht),
#    c("Total Islands", sum(sapply(islands,nrow))),
#    c("Islands >= Min Height", nrow(allpeaks)),
#    c("% Islands >= Min Height", sci.pct(nrow(allpeaks)/sum(sapply(islands,nrow)))),
#    c("Tallest Island(s)", all.tallest$MaxHt[1]),
#    c("Widest Island(s)", all.widest$Width[1]),
#    c("Largest Island(s)", all.largest$TotBp[1])





#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#######################################################################   MAIN   ########################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################





library(methods)


if (from.scratch) {
    
    ## Starting from scratch
    
    library(GenomicRanges)
    
    all.covg <- all.depth <- all.peaks <- all.CCF <- all.stats <- new.list(RDatas2[,2])
    orig.stats.ord <- c()
    keep.list <- ls()
    
    this.ca <- ca
    for (i in 1:N) {
        true.i <- i
        message(paste(RDatas2[i,"Alias"],":",system("date",intern=T)))
        load(RDatas2[i,3])   ## will overwrite 'i' // loads 'bin.means.100', 'bin.means.1000', 'bin.means.10000'
        i <- true.i
        all.covg[[i]] <- attributes(bin.means.1000)$elementMetadata$mean*1000  # 1kb-mean*1000=sum (except for the last in on the chromosome)
        all.depth[[i]] <- rowSums(bpAtDepth)
        all.CCF[[i]] <- CCF
        all.peaks[[i]] <- islands.by.chr
        x <- gsub(",","",rbind(c("IP",ifelse(RDatas2[i,"Type"]=="IP",1,0)), stats.table))
        bam.line <- grepl(".bam$",x[,2],TRUE)
        float.lines <- grepl("\\.",x[,2]) & !bam.line
        int.lines <- !float.lines & !bam.line
        all.stats[[i]] <- list(
            ## 3 elements for 3 datatypes: character, integer, float
            char=rbind2(t(RDatas2[i,]),rownamed.matrix(x[bam.line,,drop=FALSE],"character")),
            int=rownamed.matrix(x[which(int.lines),]),
            float=rownamed.matrix(x[which(float.lines),])
        )
        if (i==1) orig.stats.ord <- stats.table[,1]
        suppressWarnings(rm(bin.sums,allpeaks,bpAtDepth,CCF,L2,stats.table))  #,readlen,genome,aligns,upos,bpSeq,bpCov,L2,top1pct.pos,top1pct.seq))
    }
    ca <- this.ca
    rm(bin.means.100,bin.means.1000,bin.means.10000)
    
    message("Saving emergency session...")
    save.image(emdata)  # in case later operations blow up
    message("Saved.")
    
}


if (!from.initial) {
    
    ## Fresh start, or emergency data exists; convert to final output form if possible
    
    all.covg <- do.call(cbind, all.covg)
    all.depth <- t(merge.table.list(all.depth))
    all.CCF <- do.call(rbind, all.CCF)
    all.stats <- lapply(1:3, function(i) do.call(cbind, slice.list(all.stats,i)) )
    names(all.stats) <- c("char","int","float")
    for (i in 1:3) colnames(all.stats[[i]]) <- RDatas2$Alias
    
    message("Saving initial session...")
    rm(list=setdiff(ls(),keep.list))  # remove anything extra that came in from loading the RChipQC.binMeans.RData objects
    save.image(initdata)  # in case later operations blow up
    system(paste("rm -f",emdata))  # if we got here, then we can remove emergency data
    message("Saved.")
    
}


## Main ops at last!!

cm <- corr.mat(all.covg, plot=TRUE, main="Sample Correlations, Binned Genome Means", imgdim=c(700,600), pmar=c(10,10), filename=paste0(outprefix,".corr.mat.png"))

ndat <- chance.normalize.samples(
    as.list(as.data.frame(all.covg[,which(RDatas2$Type=="IP"),drop=FALSE])),
    as.list(as.data.frame(all.covg[,which(RDatas2$Type=="input"),drop=FALSE]))
)

palettes <- set.palettes(ndat)

png(paste0(outprefix,".cross-correlations.png"), 1000, 600)
all.stats <- multi.CCF.plot(all.CCF, all.stats, palettes, cex=1.2)  # NEW ROW ADDED TO ALL.STATS
dev.off()

png(paste0(outprefix,".enrichment.curves.png"), 600, 600)
edat <- chance.enrichment.curve.plot(ndat, palettes, cex=1.2)
dev.off()

png(paste0(outprefix,".enrichment.overlaps.png"), 1500, 1000)
ddat <- chance.enrichment.overlap.plot(ndat, edat, palettes, cex=1.2)
dev.off()

## Output of shared.peak.matrix:
##  - not yet written -

write.stats(edat, all.stats)

message("Saving final session...")
save.image(finaldata)
system(paste("rm -f",initdata))

message(paste("RChance.R",outprefix,"complete!"))

