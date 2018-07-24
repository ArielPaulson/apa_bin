#!/usr/bin/env Rscript

## SOME OF THIS IS STILL UNDER CONSTRUCTION (or outright experimental)

## MAJOR FIXME: (ADDME): ABILITY TO HANDLE > 1 ALT PER AA (since each AA is 3 NT)

source("/home/apa/apa_tools.R")

aa.props <- list(
    NP    = list( WEAK=FALSE, COLOR="olivedrab3",    NAME="Nonpolar",           AA=qw(G,A,I,L,V,P,M,F,W,Y) ),
    NPAL  = list( WEAK=FALSE, COLOR="green3",        NAME="Nonpolar-Aliphatic", AA=qw(G,A,I,L,V,P)         ),
    POLU  = list( WEAK=FALSE, COLOR="darkorange",    NAME="Polar-Uncharged",    AA=qw(C,S,T,N,Q)           ),
    ACID  = list( WEAK=FALSE, COLOR="red",           NAME="Polar-Acidic",       AA=qw(D,E)                 ),
    BASE  = list( WEAK=FALSE, COLOR="blue",          NAME="Polar-Basic",        AA=qw(H,R,K)               ),
    AROM  = list( WEAK=FALSE, COLOR="magenta",       NAME="Aromatic",           AA=qw(F,W,Y,H)             ),
    SULF  = list( WEAK=FALSE, COLOR="yellow",        NAME="Sulfur-Bearing",     AA=qw(C,M)                 ),
    PHOS  = list( WEAK=FALSE, COLOR="cyan",          NAME="Phosphorylatable",   AA=qw(S,T,Y)               ),
    AMIDE = list( WEAK=FALSE, COLOR="purple3",       NAME="Amide-Bearing",      AA=qw(W,N,Q)               ),
    ASX   = list( WEAK=TRUE,  COLOR="mediumorchid3", NAME="Asp/Asn",            AA=qw(D,N)                 ),
    GLX   = list( WEAK=TRUE,  COLOR="dodgerblue",    NAME="Glu/Gln",            AA=qw(E,Q)                 ),
    ILV   = list( WEAK=TRUE,  COLOR="darkgreen",     NAME="Leu/Ile/Val",        AA=qw(L,I,V)               ),
    OXO   = list( WEAK=TRUE,  COLOR="navy",          NAME="Oxo-Bearing",        AA=qw(D,E,N,Q)             ),
    NITRO = list( WEAK=TRUE,  COLOR="grey40",        NAME="Nitrogen-Bearing",   AA=qw(H,R,K,W,N,Q)         ),
    SMALL = list( WEAK=TRUE,  COLOR="grey75",        NAME="Small",              AA=qw(G,A,V,P,S,T,C)       ),
    POLAR = list( WEAK=TRUE,  COLOR="gold2",         NAME="Polar",              AA=qw(C,S,T,D,E,N,Q,H,R,K) )
)
aa.cons.dat <- data.frame(
    FGC=c("red","darkorange","green3","darkgreen","blue","purple4","black"),
    BLURB1=c("Identity, 100%","Identity, >=75%","Strong Chem Property, 100%","Strong Chem Property, >=75%","Weak Chem Property, 100%","Non-Conserved, but Present","Present / Absent"),
    BLURB2=c("100% Identity","75%+ Identity","100% Property","75%+ Property","100% Weak Property","No Conservation","Present/Absent")
)
aa.colors <- data.frame(
    LETTER=c("-",qw(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,X,Y),"*"),
    FGC=qw(grey,darkorange,black,white,white,black,darkorange,violet,forestgreen,white,forestgreen,black,dodgerblue,purple3,dodgerblue,white,black,black,green3,black,black,black,red),
    BGC=qw(NA,NA,yellow,red,red,violet,NA,blue,NA,blue,NA,yellow,NA,NA,NA,blue,cyan,cyan,NA,violet,grey,violet,black),
    BLURB=c("Gap","Small Nonpolar","Cys/Met","Acid","Acid","Aromatic","Small Nonpolar","Basic-Aromatic","Medium Nonpolar","Basic","Medium Nonpolar","Cys/Met","Amide","Proline","Amide","Basic","Ser/Thr","Ser/Thr","Medium Nonpolar","Aromatic","Unknown","Aromatic","Stop")
)

dna.props <- list(PUR=qw(A,G), PYR=qw(C,T))
dna.cons.dat <- data.frame(
    FGC=c("red","darkorange","green3","darkgreen","purple4","black"),
    BLURB1=c("Identity, 100%","Identity, >=75%","Transition, 100%","Transition, >=75%","Non-conserved, but Present","Present / Absent"),
    BLURB2=c("100% Identity","75%+ Identity","100% Transition","75%+ Transition","No Conservation","Present/Absent")
)
dna.colors <- data.frame(
    LETTER=c("-",qw(A,C,G,T,N)),
    FGC=qw(NA,"red","gold2","green3","blue","grey")
)


args <- commandArgs(trailing=TRUE)
align <- args[[1]]  # only required arg
imgname <- ifelse(identical("NA",args[[2]]), sub(".[mfast]+$",".png",align), args[[2]])
width <- ifelse(identical("NA",args[[3]]), 150, as.numeric(args[[3]]))
ascons <- ifelse(identical("NA",args[[4]]), FALSE, TRUE)
colors <- ifelse(identical("NA",args[[5]]), 'residues', args[[5]])   # if specified, can be 'residues' (default), 'conservation', or 'properties'
markrow <-
    if (identical("NA",args[[6]])) {
        NA
    } else {
        unlist(lapply( strsplit(args[[6]],",")[[1]], function(x){ y=suppressWarnings(as.numeric(unlist(strsplit(x,"-")))); if (length(y)>1) { y[1]:y[2] } else { y } } ) )
    }
markpos <-
    if (identical("NA",args[[7]])) {
        NA
    } else {
        unlist(lapply( strsplit(args[[7]],",")[[1]], function(x){ y=suppressWarnings(as.numeric(unlist(strsplit(x,"-")))); if (length(y)>1) { y[1]:y[2] } else { y } } ) )
    }
markalt <-
    if (identical("NA",args[[8]])) {
        NA
    } else {
        unlist(strsplit(args[[8]],","))
    }
reord <-
    if (identical("NA",args[[9]])) {
        NA
    } else {
        unlist(lapply( strsplit(args[[9]],",")[[1]], function(x){ y=suppressWarnings(as.numeric(unlist(strsplit(x,"-")))); if (length(y)>1) { y[1]:y[2] } else { y } } ) )
    }
maxflank <- ifelse(identical("NA",args[[10]]), Inf, as.numeric(args[[10]]))
alphabet <- ifelse(identical("NA",args[[11]]), "AA", args[[11]])
tabular <- falsify(as.logical(args[[12]]))
keyonly <- falsify(as.logical(args[[13]]))
makepdf <- falsify(as.logical(args[[14]]))
verbose <- falsify(as.logical(args[[15]]))

if (!(alphabet %in% qw(DNA,AA))) stop("'alphabet' must be either 'AA' or 'DNA'!\n")

message("\nAlign: ",align,"\nImg: ",imgname,"\nWidth: ",width,"\nAsCons: ",ascons,"\nColors: ",colors,"\nMarkPos: ",paste(markpos,collapse=";"),"\nMarkrow: ",paste(markrow,collapse=";"),"\nReord: ",paste(reord,collapse=";"),"\nMaxFlank: ",maxflank,"\nTabular: ",tabular,"\nKeyOnly: ",keyonly,"\nMakePDF: ",makepdf,"\n")

if (keyonly) {
    NRC <- nrow(col.mat)
    ypos <- NRC:1
    imght <- ifelse(colors == "residues", 1000, 600)
    prefix <- paste0("colorMultialign.key.",alphabet,"-",colors)
    if (makepdf) {
        pdf(paste0(prefix,".pdf"), 5, imght/10) #; par(cex=1.4)
    } else {
        png(paste0(prefix,".png"), 500, imght); par(cex=1.4)
    }
    if (colors == "residues") {
        if (alphabet == "AA") {
            col.mat <- res.colors
            null.plot(xlim=c(0,10), ylim=c(1,NRC), main="colorMultialign Color Key: 'residue'")
            rect(1,ypos-0.5, 2,ypos+0.5, col=col.mat[,3], border=NA)  # background color rectangles for each letter
            text(x=1.5, y=ypos, labels=col.mat[,1], col=col.mat[,2], font=2)
            rect(3,ypos-0.5, 9,ypos+0.5, col=col.mat[,3], border=NA)  # background color rectangles for each blurb
            text(x=3.5, y=ypos-0.05, labels=col.mat[,4], col=col.mat[,2], font=2, pos=4)
        } else {
            col.mat <- dna.colors
            null.plot(xlim=c(0,10), ylim=c(1,NRC), main="colorMultialign Color Key: 'base'")
            rect(1,ypos-0.5, 2,ypos+0.5, col=col.mat[,2], border=NA)  # background color rectangles for each letter
            text(x=1.5, y=ypos, labels=col.mat[,1], col=col.mat[,1], font=2)
        }
    } else if (colors == "conservation") {
        if (alphabet == "AA") {
            col.mat <- aa.cons.dat
            null.plot(xlim=c(0,10), ylim=c(1,NRC), main="colorMultialign Color Key: 'conservation'")
            rect(1,ypos-0.3, 2,ypos+0.3, col=col.mat[,1], border=NA)  # background color rectangles for each letter
            text(x=3, y=ypos-0.05, labels=col.mat[,2], col=col.mat[,1], font=2, pos=4)
        } else {
            col.mat <- dna.cons.dat
            null.plot(xlim=c(0,10), ylim=c(1,NRC), main="colorMultialign Color Key: 'conservation'")
            rect(1,ypos-0.3, 2,ypos+0.3, col=col.mat[,1], border=NA)  # background color rectangles for each letter
            text(x=3, y=ypos-0.05, labels=col.mat[,2], col=col.mat[,1], font=2, pos=4)
        }
    } else if (colors == "properties") {
        if (alphabet == "AA") {
            stop("color scheme 'properties' not ready yet!\n")
        } else {
            stop("No color scheme 'properties' for alphabet 'DNA'!\n")
        }
    } else {
        stop("color scheme must be one of 'residues', 'conservation', or 'properties'!\n")
    }
    dev.off()
    quit()
}

mfa <- t(sapply(read.fasta(align),function(x) unlist(strsplit(x,'')) ))
mfa <- mfa[,!apply(mfa=="-",2,all)]  # 'muscle' aligner gets stupid sometimes and returns columns pure-gap columns, i.e. with no NTs/AAs (fortunately, this is very rare)

do.markrow <- ifelse(is.na(markrow[1]), FALSE, TRUE)
do.markpos <- ifelse(is.na(markpos[1]), FALSE, TRUE)
do.markalt <- ifelse(is.na(markalt[1]), FALSE, TRUE)
do.classalt <- ifelse(do.markalt & tabular & do.markrow & length(markrow)==1, TRUE, FALSE)

if (do.markpos) {
    markpos2 <- rep("",ncol(mfa))
    if (do.markalt) {
        markpos2[markpos] <- markalt  # the alternates
    } else {
        markpos2[markpos] <- "#"  # the default
    }
}

if (!is.na(reord[1])) {
    ## if reordering, rows may be dropped. Remove any all-gap columns that result.
    mfa <- mfa[reord,]
    if (do.markrow) markrow <- match(markrow, reord)
    keep.pos <- rep(TRUE,ncol(mfa))
    for (i in 1:ncol(mfa)) {
        if (all(mfa[,i]=="-")) keep.pos[i] <- FALSE
    }
    mfa <- mfa[,keep.pos]
    if (do.markrow) markrow <- match(markrow, reord)
    if (do.markpos) {
        markpos3 <- markpos2[keep.pos]
        if (length(markpos3) != length(markpos2)) {
            lostpos <- markpos2[!keep.pos]
            message("WARNING: 'markpos' vector included positions which were 100% gap! ",paste(lostpos,collapse=","),"\n")
        }
        markpos2 <- markpos3
        markpos <- which(markpos2=="*")
    }
}

NC <- ncol(mfa)
NR <- nrow(mfa)
supermaj <- NR*0.75

## original sequence lengths
len <- as.matrix(apply(mfa!="-",1,sum))
if (verbose) print(len)

## matrix of non-gap identity positions between sequences
ident <- matrix(NA, NR, NR, FALSE, list(rownames(mfa),rownames(mfa)))
for (i in 1:NR) {
    for (j in 1:NR) {
        ident[i,j] <- sum(mfa[i,]==mfa[j,] & mfa[i,]!="-")
    }
}
if (verbose) print(ident)


colN=list(rep(0,NC))
colC=list(rep("",NC))
cons <- cons2 <- cons.col <- rep("",NC)
classify <- as.data.frame(c(Pos=list(1:NC),as.list(as.data.frame(t(mfa))),Ref.Class=colN,Ref.Comments=colC,Alt.Obs=colC,Alt.Class=colN,Alt.Comments=colC))
## most of the fields in 'classify' cater to snpReveal output, as this program was developed as the MA visualizer/annotator tools for the snpReveal script



classify.AA <- function(AAs) {
    
    tab <- table(AAs[AAs!="-"])
    consensus <- names(tab)[which.max(tab)]
    sum.ident <- sum(AAs==consensus)
    why <- c()
    
    if (sum.ident==NR) {
        ## all identical to residue
        class.n <- 1
    } else if (any(AAs=="-")) {
        ## residue position not found in all proteins
        class.n <- 7
    } else {
        ## all present, but not 100% identical, so look for other trends
        sum.prop      <- sapply(props,      function(x) sum(AAs %in% x$AA) )
        sum.prop.weak <- sapply(props.weak, function(x) sum(AAs %in% x$AA) )
        
        if (any(sum.prop==NR)) {
            ## not identically conserved, but some set of key chemical properties are fully conserved
            class.n <- 3
            why <- names(props)[sum.prop==NR]
        } else if (sum.ident>=supermaj) {
            ## neither identically nor chemically conserved 100%, but,
            ## supermajority exists for identity, AND residue is part of that supermajority
            class.n <- 2
        } else if (any(sum.prop>=supermaj)) {
            ## supermajority does not exist for identity, but,
            ## supermajority exists for some chemical property, AND residue is part of that supermajority
            class.n <- 4
            why <- names(props)[sum.prop>=supermaj]
        } else if (any(sum.prop.weak==NR)) {
            ## supermajority does not exist for any strong property, but,
            ## all residues belong to a "weak property" conservation class, see defs above
            class.n <- 5
            why <- names(props.weak)[sum.prop.weak==NR]
        } else {
            ## then, not conserved at all, but at least present in 100% proteins
            class.n <- 6
        }
    }
    
    label <- aa.cons.dat$BLURB2[class.n]
    if (length(why)>0) {
        why <- why[order(sapply(aa.props[why],function(x) length(x$AA) ))]  # list candidate conservation classes from smallest (most rigorous) to largest (least rigorous)
        label <- paste0(label,": ",paste(why,collapse=","))
    }
    
    list(CLASS=class.n, LABEL=label, COLOR=aa.cons.dat$FGC[class.n], CONS=consensus, IDENT=sum.ident)
}


classify.DNA <- function(NTs) {
    
    tab <- table(NTs[NTs!="-"])
    consensus <- names(tab)[which.max(tab)]
    sum.ident <- sum(NTs==consensus)
    why <- c()
    
    if (sum.ident==NR) {
        ## all identical to base
        class.n <- 1
    } else if (any(NTs=="-")) {
        ## base position not found in all transcripts
        class.n <- 6
    } else {
        ## all present, but not 100% identical, so look for other trends
        sum.prop <- sapply(props, function(x) sum(NTs %in% x) )
        
        if (any(sum.prop==NR)) {
            ## not identically conserved, but all belong to a single transition pair (A-G or C-T)
            class.n <- 3
            why <- names(props)[sum.prop==NR]
        } else if (sum.ident>=supermaj) {  
            ## neither identically nor transitionally conserved 100%, but,
            ## supermajority exists for identity
            class.n <- 2
        } else if (any(sum.prop>=supermaj)) {
            ## supermajority does not exist for identity, but,
            ## supermajority exists for one transition pair
            class.n <- 4
            why <- names(props)[sum.prop>=supermaj]
        } else {   # not conserved at all, but at least present in 100% sequences
            class.n <- 5
        }
    }
    
    label <- aa.cons.dat$BLURB2[class.n]
    if (length(why)>0) {
        why <- sort(why)  # don't need to sort DNA classes by size; there are only two and they are both the same size
        label <- paste0(label,": ",paste(why,collapse=","))
    }
    
    list(CLASS=class.n, LABEL=label, COLOR=dna.cons.dat$FGC[class.n], CONS=consensus, IDENT=sum.ident)
}



if (alphabet == "AA") {
    
    ## AMINO-ACID ALPHABET HANDLING
    res.colors <- aa.colors
    props <- aa.props[!sapply(aa.props,"[[","WEAK")]
    props.weak <- aa.props[sapply(aa.props,"[[","WEAK")]
    P  <- length(props)
    Pw <- length(props.weak)
    
    for (i in 1:NC) {
        ## first classify the reference residues
        AAs <- mfa[,i]
        x <- classify.AA(AAs)
        if (length(x$CONS)==0) print(x)
        cons[i] <- x$CONS
        cons.col[i] <- x$COLOR
        classify[i,NR+2] <- x$CLASS
        classify[i,NR+3] <- x$LABEL
        if (ascons) mfa[AAs==cons[i],i] <- "."
        
        ## then classify the alt residues, if any
        if (do.classalt & i %in% markpos) {
            alt <- markalt[which(markpos==i)]
            AAs2 <- AAs
            if (length(alt)==1) {
                AAs2[markrow] <- alt    # replace ref with alt and re-classify
                x <- classify.AA(AAs2)
                classify[i,NR+4] <- ifelse(alt %in% AAs, "Yes", "No")  # alt observed in this position in other sequence(s)  # NOT AAs2, THAT ALWAYS CONTAINS ALT
                classify[i,NR+5] <- x$CLASS
                classify[i,NR+6] <- x$LABEL
            } else {
                ## ## THIS IS A PROBLEM -- MULTI-ALT POSITIONS
                AAs2[markrow] <- paste(alt,collapse=",")    # replace ref with BOTH alts, but for now, this means we cannot classify...
                ##x <- classify.AA(AAs2)  # not yet
                classify[i,NR+4] <- paste(c("No","Yes")[0+c(alt %in% AAs)], collapse=",")  # alt observed in this position in other sequence(s)  # NOT AAs2, THAT ALWAYS CONTAINS ALT
                classify[i,NR+5] <- NA
                classify[i,NR+6] <- "Multi-Alt"
            }
            alt.done <- FALSE
        }
    }
    
} else if (alphabet == "DNA") {
    
    ## DNA ALPHABET HANDLING
    res.colors <- dna.colors
    props <- dna.props
    
    for (i in 1:NC) {
        ## first classify the reference bases
        NTs <- mfa[,i]
        x <- classify.DNA(NTs)
        cons[i] <- x$CONS
        cons.col[i] <- x$COLOR
        classify[i,NR+2] <- x$CLASS
        classify[i,NR+3] <- x$LABEL
        if (ascons) mfa[NTs==cons[i],i] <- "."
        
        ## then classify the alt bases, if any
        if (do.classalt & i %in% markpos) {
            alt <- markalt[which(markpos==i)]
            NTs2 <- NTs
            if (length(alt)==1) {
                NTs2[markrow] <- alt    # replace ref with alt and re-classify
                x <- classify.NT(NTs2)
                classify[i,NR+4] <- ifelse(alt %in% NTs, "Yes", "No")  # alt observed in this position in other sequence(s)  # NOT NTs2, THAT ALWAYS CONTAINS ALT
                classify[i,NR+5] <- x$CLASS
                classify[i,NR+6] <- x$LABEL
            } else {
                ## ## THIS IS A PROBLEM -- MULTI-ALT POSITIONS
                NTs2[markrow] <- paste(alt,collapse=",")    # replace ref with BOTH alts, but for now, this means we cannot classify...
                ##x <- classify.NT(NTs2)  # not yet
                classify[i,NR+4] <- paste(c("No","Yes")[0+c(alt %in% NTs)], collapse=",")  # alt observed in this position in other sequence(s)  # NOT AAs2, THAT ALWAYS CONTAINS ALT
                classify[i,NR+5] <- NA
                classify[i,NR+6] <- "Multi-Alt"
            }
            alt.done <- FALSE
        }
    }
    
}

if (tabular) write.table(classify, paste0(imgname,".table.txt"), sep="\t", quote=FALSE, row.names=FALSE)

if (ascons) {
    mfa <- rbind2(mfa,matrix(cons,nrow=1))
    rownames(mfa)[NR+1] <- "CONSENSUS"
    NR <- NR + 1
    mfa[mfa=="-"] <- ""
}

mfa <- rbind2( matrix(unlist(lapply(seq(0,NC,10), function(x){ y=unlist(strsplit(as.character(x),'')); c(y,rep("",10-length(y))) }))[1:NC+1],nrow=1), mfa )
mfa[1,1] <- 1
NR <- NR + 1
rownames(mfa)[1] <- "POSITION"
if (do.markrow) markrow <- markrow+1   # +1 because POSITION row has been added

if (do.markpos) {
    IM(dim(mfa))
    IM(length(markpos2))
    mfa <- rbind2(mfa, markpos2)
    NR <- NR + 1
    rownames(mfa)[NR] <- "KEY RESIDUES"
}

mfa <- cbind("",mfa); colnames(mfa) <- 0:NC  # the zero position, not displayed; a formatting crutch
NC <- NC + 1
cons.col <- c(1,cons.col)
if (ncol(mfa)<width) width <- ncol(mfa)

cmat <- bmat <- matrix(NA, NR, NC)
if (colors == "residues") {
    if (alphabet == "AA") {
        for (i in 1:nrow(res.colors)) {
            cmat[mfa==res.colors[i,1]] <- res.colors[i,2]
            bmat[mfa==res.colors[i,1]] <- res.colors[i,3]
        }
    } else {
        for (i in 1:nrow(res.colors)) {
            cmat[mfa==res.colors[i,1]] <- res.colors[i,2]
        }
    }
} else if (colors == "conservation") {
    if (alphabet == "AA") {
        for (i in 1:NC) cmat[,i] <- cons.col[i]
    } else {
        for (i in 1:NC) cmat[,i] <- cons.col[i]
    }
} else if (colors == "properties") {  # not ready
    if (alphabet == "AA") {
#        for (i in 1:length(aa.prop.colors)) cmat[,i] <- cons.col[i]
    } else {
#        stop()
    }
}

black.rows <- 1
if (do.markpos) black.rows <- c(black.rows, NR)  # add marked-pos row to the black-rows list
for (i in black.rows) {
    cmat[i,truthify(mfa[i,]!="")] <- 1  # first, last rows do not use custom colors
    bmat[i,] <- rep(NA,NC)
}

cmat[mfa=="."] <- 1
bmat[mfa=="."] <- NA





gap <- 2
expand <- c(0.5,0.5)

cex <- 1.2
ppa <- c(16,12)  # ifelse(NR<=20,12,30))  # x, y pixels-per-aspect
ppc <- 5.4 + round((cex-1)/0.2,0)
rowname.pix <- ppc * max(nchar(rownames(mfa)))
mar2 <- rowname.pix*1.1/10

if (is.na(maxflank)) {
    
    windows <- list(1:ncol(mfa))
     mfa <- list(mfa)
    cmat <- list(cmat)
    bmat <- list(bmat)
    B <- ceiling(NC/width)
    blocks <- cbind( width*(0:(B-1))+1, width*1:B )
    if (blocks[B,2] > NC) blocks[B,2] <- NC
    blocks <- list(blocks)
    
} else {
    
    windows <- lapply(markpos, function(i) (i-maxflank):(i+maxflank) )
     mfa <- lapply(windows, function(x)  mfa[,x+1] )  # +1 because col 1 of 'mfa' is 0
    cmat <- lapply(windows, function(x) cmat[,x+1] )
    bmat <- lapply(windows, function(x) bmat[,x+1] )
    blocks <- lapply(windows, function(x) {
        lx <- length(x)
        if (lx > width) {
            B <- ceiling(lx/width)
            b <- cbind( width*(0:(B-1))+1, width*1:B )
        } else {
            B <- 1
            b <- cbind( lx*(0:(B-1))+1, lx*1:B )
        }
        if (b[B,2] > lx) b[B,2] <- lx
        b
    })
    
}

B <- lapply(blocks,nrow)
W <- length(windows)
init <- ylines <- blocks

for (w in 1:W) {
    x <- cbind(BLOCK=rep(1:B[[w]],each=NR+gap), MFA=rep(c(1:NR,rep(NA,gap)),B[[w]]), PNG=((NR+gap)*B[[w]]):1)
    x[falsify(x[,2]==1),3] <- x[falsify(x[,2]==1),3] + 0.25  # bump up position row slightly
    y <- do.call(rbind, lapply(1:B[[w]], function(b) mfa[[w]][x[x[,1]==1,2],1:2,drop=FALSE] ))
    y[,1] <- 0;  y[,2] <- x[,3];  mode(y) <- "numeric"
    ylines[[w]] <- x
    init[[w]] <- y
}

if (do.markrow) {
      marked <- lapply(ylines, function(x) which(x[,2] %in% markrow) )
    unmarked <- lapply(ylines, function(x) setdiff(1:nrow(x), marked) )
} else {
    unmarked <- lapply(ylines, function(x) 1:nrow(x) )
}





save.image(paste0(imgname,".RData"))

for (w in 1:W) {
    
    imgname.w <- ifelse( is.na(maxflank), imgname, sub("png$",paste0(min(windows[[w]]),"-",max(windows[[w]]),".png"),imgname) )
    lw <- length(windows[[w]])
    if (lw<=width) {
        img.ncol <- lw
        ## drop the block-spacer rows, if only one block in the image
        init[[w]] <- init[[w]][1:(nrow(init[[w]])-2),]
        init[[w]][,2] <- init[[w]][,2]-2
        ylines[[w]] <- ylines[[w]][1:(nrow(ylines[[w]])-2),]
        ylines[[w]][,3] <- ylines[[w]][,3]-2
        unmarked[[w]] <- unmarked[[w]][1:(length(unmarked[[w]])-2)]
    } else {
        img.ncol <- width
    }
    pngwd <- rowname.pix+img.ncol*ppa[1]
    #print(ylines[[w]]); print(init[[w]])
    
    if (makepdf) {
        pngwd <- round(10*pngwd/pnght,1)
        pnght <- 10
        message("Width: ",pngwd,", Height: ",pnght)
        pdf(sub("png$","pdf",imgname.w), pngwd, pnght)
        par(mar=c(1,mar2,4,1), family="mono", xaxs="i", yaxs="i", cex=0.75)
    } else {
        pnght.w <- ifelse(W==1,50,70)+nrow(ylines[[w]])*ppa[2]
        message("Width: ",pngwd,", Height: ",pnght.w)
        png(imgname.w, pngwd, pnght.w)
        par(mar=c(1,mar2,4,1), family="mono", xaxs="i", yaxs="i", cex=cex)
    }
    ##null.plot(xlim=c(0,img.ncol+1), ylim=c(min(ylines[[w]][,3])-1,max(ylines[[w]][,3])+1))
    
    plot(init[[w]], pch="", main="", xlab="", ylab="", axes=FALSE, xlim=c(0,img.ncol+1), ylim=c(min(ylines[[w]][,3])-1,max(ylines[[w]][,3])+1))

    #IM(w, nrow(ylines[[w]]), length(marked[[w]]), length(unmarked[[w]]))
    if (do.markrow) axis(2, at=ylines[[w]][marked[[w]],3], tick=FALSE, labels=rownames(init[[w]])[marked[[w]]], las=2, font=2)
    axis(2, at=ylines[[w]][unmarked[[w]],3], tick=FALSE, labels=rownames(init[[w]])[unmarked[[w]]], las=2, font=1)
    axis(2, at=ylines[[w]][,3], tick=TRUE, labels=rep("",nrow(ylines[[w]])))  # don't want ticks bolded, so do them all separately

    for (b in 1:B[[w]]) {
        irange <- (blocks[[w]][b,1]):(blocks[[w]][b,2])
        iwidth <- length(irange)
        irange2 <- 1:iwidth
        mfab <-  mfa[[w]][,irange,drop=FALSE]
        cmatb <- cmat[[w]][,irange,drop=FALSE]
        bmatb <- bmat[[w]][,irange,drop=FALSE]
        yl <- ylines[[w]][ylines[[w]][,1]==b,]
        for (n in 1:nrow(yl)) {
            r <- yl[n,2]
            y <- yl[n,3]
            if (is.na(r)) next
            symb <- mfab[r,,drop=FALSE]
            rowcols <- dotcols <- cmatb[r,,drop=FALSE]
            if (ascons) {
                dots <- rep(NA,iwidth)
                dots[mfab[r,,drop=FALSE]=="."] <- 20  # pch for dot positions
                ##dotcols[mfab[r,]!="."] <- 0  # do not display
                rowcols[mfab[r,,drop=FALSE]=="."] <- 0  # ditto
            }
            if (colors=="residues") {
                use <- which(!is.na(bmatb[r,,drop=FALSE]))
                for (i in use) rect(irange2[i]-expand[1],y-expand[2], irange2[i]+expand[1],y+expand[2], col=bmatb[r,i], border=NA)  # background color rectangles for each letter
            }
            text(x=irange2, y=rep(y,iwidth), labels=symb, col=rowcols, font=2)
            if (ascons) points(x=irange2, y=rep(y,iwidth), pch=dots, col=dotcols, cex=0.2)  # consensus "dot"
        }
    }
    dev.off()
}

readme <- c(
    "\nExplanation of classifications and colors used by colorMultiAlign:",
    "\nAmino Acid Groups and Attributes:\nPROPERTY\tLONG.NAME\tWEAK?\tCOLOR\tRESIDUES",
    sapply(1:length(aa.props), function(i){ x=aa.props[[i]]; paste(c(names(aa.props)[i], x$NAME, x$WEAK, x$COLOR, paste(sort(x$AA),collapse="")), collapse="\t") }),
    "\nAmino Acid Conservation Levels:\nLEVEL\tCONSERVATION\tCOLOR",
    sapply(1:length(aa.cons.dat[[1]]), function(i) paste(c(i,aa.cons.dat[[2]][[i]],aa.cons.dat[[1]][[i]]), collapse="\t") ),
    "\nAmino Acid Plot Colors:\nAA\tFG.COL\tBG.COL\tNOTES",
    sapply(1:length(aa.colors[[1]]), function(i) paste(translate(sapply(aa.colors,"[",i),"NA","white"), collapse="\t") ),
    "\nNucleotide Groups (Transition Groups):\nPROPERTY\tLONG.NAME\tBASES",
    sapply(1:length(dna.props), function(i){ x=dna.props[[i]]; paste(names(dna.props)[i], ifelse(i==1,"Purines","Pyrimidines"), x, collapse="\t") }),
    "\nNucleotide Conservation Levels:\nLEVEL\tCONSERVATION\tCOLOR",
    sapply(1:length(dna.cons.dat[[1]]), function(i) paste(c(i,dna.cons.dat[[2]][[i]],dna.cons.dat[[1]][[i]]), collapse="\t") ),
    "\nNucleotide Plot Colors:\nNT\tFG.COL\tBG.COL",
    sapply(1:length(dna.colors[[1]]), function(i) paste(translate(c(sapply(dna.colors,"[",i),"white"),"NA","white"), collapse="\t") ),
    ""
)
write.vector(readme, paste0(imgname,".readme.txt"))
