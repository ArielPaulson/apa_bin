
## Some libraries

#suppressMessages(require(GenomicRanges))
suppressMessages(require(Biostrings))   # BLOSUM100
data(BLOSUM100)

## Some objects

snpReveal_subunit  <- "/home/apa/local/bin/scriptutils/snpReveal3_subunit.R"
snpReveal_finalize <- "/home/apa/local/bin/scriptutils/snpReveal3_finalize.R"

angstroms <- c(4,3)  # for PDB distance filtering: c( max distance to report, max distance to call interaction )

codons <- rbind(
    c('TTT','F'),c('TTC','F'),c('TTA','L'),c('TTG','L'),
    c('TCT','S'),c('TCC','S'),c('TCA','S'),c('TCG','S'),
    c('TAT','Y'),c('TAC','Y'),c('TAA','*'),c('TAG','*'),
    c('TGT','C'),c('TGC','C'),c('TGA','*'),c('TGG','W'),
    c('CTT','L'),c('CTC','L'),c('CTA','L'),c('CTG','L'),
    c('CCT','P'),c('CCC','P'),c('CCA','P'),c('CCG','P'),
    c('CAT','H'),c('CAC','H'),c('CAA','Q'),c('CAG','Q'),
    c('CGT','R'),c('CGC','R'),c('CGA','R'),c('CGG','R'),
    c('ATT','I'),c('ATC','I'),c('ATA','I'),c('ATG','M'),
    c('ACT','T'),c('ACC','T'),c('ACA','T'),c('ACG','T'),
    c('AAT','N'),c('AAC','N'),c('AAA','K'),c('AAG','K'),
    c('AGT','S'),c('AGC','S'),c('AGA','R'),c('AGG','R'),
    c('GTT','V'),c('GTC','V'),c('GTA','V'),c('GTG','V'),
    c('GCT','A'),c('GCC','A'),c('GCA','A'),c('GCG','A'),
    c('GAT','D'),c('GAC','D'),c('GAA','E'),c('GAG','E'),
    c('GGT','G'),c('GGC','G'),c('GGA','G'),c('GGG','G')
)
#c('TAA','#'),c('TAG','$'),c('TGA','%'),

TLAs <- rbind(
    c('ALA','A'),c('CYS','C'),c('ASP','D'),c('GLU','E'),
    c('PHE','F'),c('GLY','G'),c('HIS','H'),c('ILE','I'),
    c('LYS','K'),c('LEU','L'),c('MET','M'),c('ASN','N'),
    c('PRO','P'),c('GLN','Q'),c('ARG','R'),c('SER','S'),
    c('THR','T'),c('VAL','V'),c('TRP','W'),c('TYR','Y'),
    c('SEL','C')   # this last one just in case???  DO NOT convert Selenocysteines to U, that will break stuff.  Retain as Cysteines.
)

## Custom functions


align.orthologs <- function() {
    
    ## Filter down a set of orthologous transcript sequences for the most-alignable subset
    ## Produce a multiple alignment for this subset and visualize it
    ## Classify the mutant position(s)
    
    tmp <- new.list(qw(ortho,aligns,exons,vars,vars2))
    tmp$vars <- this.vars    # will upgrade and return
    tmp$ortho <- new.list(names(orgs))
    
    ## Screen ortholog peptides for similarity to target transcript
    for (j in 1:O) {
        
        g.ortho <- this.ortho[[j]]$ortho
        g.data <- this.ortho[[j]]$data
        peps <- this.ortho[[j]]$peps
        
        if (nrow(peps)>0) {
            ## generate multiple alignment scores; this.tseq$aa is the reference
            pa <- sapply(1:nrow(peps), function(g) pairwiseAlignment(this.tseq$aa, peps$peptide[g], substitutionMatrix=BLOSUM100, gapOpening=0, gapExtension=-5) )
            peps$aln_score <- sapply(pa, function(x) attributes(x)$score )
            ## which have max score (allow ties)
            ms <- zapsmall(peps$aln_score)==max(zapsmall(peps$aln_score))
            ## get aligned sequences, as vectors
            qry <- lapply(pa, function(x) unlist(strsplit(as.character(attributes(x)$pattern),"")) )
            sbj <- lapply(pa, function(x) unlist(strsplit(as.character(attributes(x)$subject),"")) )
            ## which have sufficient identity%
            peps$idt_pct <- sapply(1:length(qry), function(k) sum(qry[[k]]==sbj[[k]])/sum(qry[[k]]!="-") )
            si <- peps$idt_pct >= min.ortho.idt
            ## which have both
            ok <- which(falsify(ms & si))
            if (length(ok)>0) {
                ## if > 1 have both, take first
                peps$acceptable[ok[1]] <- TRUE
                ## drop rest
                peps <- peps[peps$acceptable,]
            } else {
                ## retain one best ortholog for reference purposes, but do not flag as acceptable
                peps <- peps[which(ms)[1],]
            }
            
            ## Generate dataframe
            df <- if (j==1) {
                      cbind(peps[,1:3], g.data[match(peps[,3],g.data[,1]),2:ncol(g.data)], aa_len=0, peps[,4:7])
                  } else {
                      cbind(peps[,1:3], g.data[match(peps[,3],g.data[,1]),2:ncol(g.data)], g.ortho[match(peps[,3],g.ortho[,1]),!mgrepl(qw(chrom,transcript,peptide),colnames(g.ortho))], aa_len=0, peps[,4:7])
                  }
            ## fill peptide-length column
            df$aa_len <- nchar(df$peptide)
            ## store INITIAL org-wise dataframes -- rbound into one dataframe below
            ## further on, will upgrade & return
            tmp$ortho[[j]] <- df
        }
    }
    
    
    ## Want all org dfs rbound together, but non-source dfs have more columns than source df.
    ## This is because non-sources have ortholog data while source does not.
    ## Must reshape source df to match non-source (add some NAs, basically).
    ## Then rbind into final df.
    
    ## rbind all non-source dfs that actually have data (=found orthologs)
    ok <- setdiff(which(listLengths(tmp$ortho)>0), 1)  # remove source entry (#1)
    x <- do.call(rbind2,tmp$ortho[ok])
    ## convert source df to list and align to non-source col names
    y <- as.list(tmp$ortho[[1]])[match(colnames(x),colnames(tmp$ortho[[1]]))]
    ## rename source cols to whatever non-source is using
    names(y) <- colnames(x)
    ## fill missing positions with NAs
    y[listLengths(y)==0] <- NA
    ## add source as first df row
    tmp$ortho <- rbind2(as.data.frame(y),x)
    rownames(tmp$ortho)[1] <- names(orgs)[1]
    ## source row lacks rowname; it will be the one org missing from current rownames
    message("Lack orthologs: ",paste(setdiff(names(orgs),rownames(tmp$ortho)),collapse=" "))
    message("Have orthologs: ",paste(setdiff(rownames(tmp$ortho),names(orgs)[1]),collapse=" "))
    
    ## Send final sequnce set to multiple alignment, PNG generation, and position classification
    ## MAs courtesy of muscle
    ## PNG/classification courtesy of /home/apa/local/bin/colorMultialign
    
    allseq <- tmp$ortho$peptide
    names(allseq) <- rownames(tmp$ortho)
    f0 <- paste0(this.tmp.prefix,"ortho.all.fa")
    write.fasta(allseq, f0)
    
    okseq <- allseq[truthify(tmp$ortho$acceptable)]
    names(okseq) <- rownames(tmp$ortho)[tmp$ortho$acceptable]
    O2 <- length(okseq)  # may have lost some orgs
    message("Good orthologs: ",paste(setdiff(names(okseq),names(orgs)[1]),collapse=" "))
    
    if (O2<2) {
        
        message("Too few sequences to proceed...\n")
        return()
        
    } else {
        
        message(this.trans,": aligning orthologs...")
        ref.is.neg <- tmp$ortho$strand[1] == -1
        if (!is.na(multialn.ortho)) {
            tmp$aligns <- do.call(cbind, lapply(read.fasta(multialn.ortho), function(x) unlist(strsplit(x,"")) ))
            f2 <- multialn.ortho
        } else {
            f1 <- paste0(this.tmp.prefix,"ortho.muscle_in.fa")
            f2 <- paste0(this.tmp.prefix,"ortho.muscle_out.fa")
            ##save.image(RData)
            if (!file.exists(f1)) write.fasta(okseq, f1)
            if (file.exists(f2)) {
                tmp$aligns <- do.call(cbind, strsplit(read.fasta(f2),""))
            } else {
                tmp$aligns <- t(muscle(f1, f2, as.matrix=TRUE, fix.terminals="AA"))
            }
        }
        mark.row <- which(colnames(tmp$aligns)==names(orgs)[1])
        
        ## map exon positions for transcripts.  initial: matrices for (-) strand genes will have to be row-reversed
        message("Mapping exon structure...")
        tmp$exons <- this.exons
        
        ## Convert genomic position to AA position in source transcript
        ## Then map position in source transcript to MA position -> mark.pos
        
        if (ref.is.neg) tmp$exons <- tmp$exons[nrow(tmp$exons):1,]
        
        for (r in 1:nrow(tmp$vars)) {
            x <- cbind(tmp$vars$Pos[r]-tmp$exons[,3], tmp$exons[,4]-tmp$vars$Pos[r])
            wx <- which(x[,1]>0 & x[,2]>0)   # row of tmp$exons which carries mutation #r
            message(paste(r,this.trans,wx,tmp$ortho$strand[1],nrow(x)))
            if (length(wx)==0) stop("Variant #",r," at position ",tmp$vars$Pos[r]," is not found on the exons of transcript ",this.trans,"!")
            tpos.nt <- this.vars$TransBp[r]
            tpos.aa <- this.vars$TransAa[r]
            ma.pos <- which(tmp$aligns[,mark.row]!="-")[tpos.aa]
            nt.phase <- (tpos.nt%%3)-1
            if (nt.phase==-1) nt.phase <- 2
            message("NTP: ",tpos.nt," | PHA: ",nt.phase," | AAP: ",tpos.aa," | MAP: ",ma.pos)
            tmp$vars$Phase[r] <- nt.phase
            tmp$vars$MAPos[r] <- ma.pos
            tmp$vars$OrthoN[r] <- ncol(tmp$aligns)-1
            tmp$vars$OrthoR[r] <- sum(tmp$aligns[ma.pos,]!="-")-1
            codon.start <- tpos.nt-nt.phase
            codon.end <- codon.start+2
            codon.seq <- this.tseq2$nt[codon.start:codon.end]
            codon.seq1 <- paste(codon.seq,collapse="")
            message(paste("REF:",codon.seq1,codons[match(codon.seq1,codons[,1]),2],paste0("(from ",this.tseq2$aa[tpos.aa],")")))
            codon.seq[nt.phase+1] <- ifelse(ref.is.neg, tmp$vars$AltNT.T[r], tmp$vars$AltNT.G[r])
            message(paste("ALT:",codon.seq1,codons[match(codon.seq1,codons[,1]),2]))
        }
        tmp$vars$RefNT.T[tmp$vars$RefNT.T==0] <- ""  # annotation errors
        ok <- which(tmp$vars$RefNT.T != "")
        mark.dat <- unique(data.frame(Pos=tmp$vars$MAPos[ok], Alt=tmp$vars$AltAA[ok]))
        mark.pos <- paste(mark.dat$Pos,collapse=",")
        mark.alt <- paste(mark.dat$Alt,collapse=",")
        message(mark.row," : ",mark.pos," : ",mark.alt)
        
        this.tmp.prefix.f2 <- sub("out.fa$","",f2)
        cmd <- paste("/home/apa/local/bin/colorMultialign -f",f2,"-o",paste0(this.tmp.prefix.f2,"res.png"),"-c residues --as-cons -mr",mark.row,"-mp",mark.pos,"-ma",mark.alt,"-w 151 --max-flank 75")
        message(cmd)
        system(cmd)  # could capture STDERR, for something
        
        cmd <- paste("/home/apa/local/bin/colorMultialign -f",f2,"-o",paste0(this.tmp.prefix.f2,"cons.png"),"-c conservation --tabular -mr",mark.row,"-mp",mark.pos,"-ma",mark.alt,"-w 151 --max-flank 75")
        message(cmd)
        system(cmd)
        ## gthumb(paste0(this.tmp.prefix,"cons.png"))
        
        ## Finally, classify marked positions in MA -- have colorMultialign produce tabular output, rows = classified positions?
        ## Include AA substitution, if not splice effect
        classify <- read.delim(paste0(this.tmp.prefix.f2,"cons.png.table.txt"), as.is=TRUE)
        class2vars <- classify[tmp$vars$MAPos[ok],c(O2+c(2:6),2:(O2+1))]
        
        if (nrow(class2vars)<nrow(tmp$vars)) {
            ## some variant positions were not classifiable?  error in data entry?
            for (r in (nrow(class2vars)+1):length(ok)) class2vars[r,] <- NA
        }
        tmp$vars2 <- cbind(tmp$vars, class2vars)
        
        tmp
    }
}










process.pdb <- function(tmp) {
    
    ## if PDB file specified, and PDB(s) assigned to this transcript, analyze PDB for side-chain proximity to other residues
    
    p <- 1  # currently can only report on 1 PDB file per transcript, even if multiple given
    chain <- this.trans.pdb$Chain[p]
    pdbname <- paste0("PDB:",this.trans.pdb$ID[p])
    pdb <- this.trans.pdb.dat[[p]]
    tmp <- new.list(qw(vars2,report))
    tmp$report <- new.list(this.trans.pdb$ID[p])
    
    ## Extract Expected AA Sequence
    pdb.exp <- pdb[which(pdb[,1]=="SEQRES")[grepl(paste0("^ +[0-9]+ ",chain),pdb[pdb[,1]=="SEQRES",2])],2]
    pdb.exp.len <- as.numeric(unlist(strsplit(pdb.exp[1]," +"))[4])
    pdb.exp <- real(unlist(lapply(pdb.exp, function(x) translate(unlist(strsplit(x," +"))[5:17], TLAs[,1], TLAs[,2]) )))
    pdb.exp2 <- paste(pdb.exp, collapse="")
    if (length(pdb.exp)!=pdb.exp.len) message("WARNING: ",pdbname," reported sequence length ",pdb.exp.len," does not match resulting length ",nchar(pdb.exp),"!")
    
    ## Extract all atoms and het atoms 
    atoms <- rbind( pdb[which(pdb[,1]=="ATOM"),], pdb[which(pdb[,1]=="HETATM"),] )
    atoms <- data.frame( RECORD=atoms[,1], rownameless(t(sapply(atoms[,2], function(x) c(SERIAL=gsub(" ","",substr(x,1,5)), DATA=substr(x,6,80)) ))) )
    atoms[[2]] <- as.numeric(atoms[[2]])
    
    ## Get distance matrix of all atoms
    apos <- cbind(atoms[,1:2], rownameless(t(sapply(atoms[[3]], function(x) real(as.numeric(unlist(strsplit(substr(x,22,45)," +")))) ))))
    colnames(apos)[3:5] <- qw(X,Y,Z)
    adist <- as.matrix(dist(apos[,3:5]))
    
    ## Get residue identity of all non-het atoms
    ares <- cbind(atoms[atoms[,1]=="ATOM",1:2], rownameless(t(sapply(atoms[[3]][atoms[,1]=="ATOM"], function(x) c(ATOM=gsub(" ","",substr(x,2,5)), RES=substr(x,7,9), POS=as.numeric(substr(x,12,15)), CHAIN=substr(x,11,11)) ))))
    ares[,4] <- translate(ares[,4], TLAs[,1], TLAs[,2])
    for (j in c(2,5)) mode(ares[[j]]) <- "numeric"
    
    ## Extract Observed AA Sequence
    pdb.obs <- unique(ares[ares[,6]==chain,4:5])[,1]
    pdb.obs2 <- paste(pdb.obs, collapse="")
    pdb.obs.len <- length(pdb.obs)
    
    ## Align observed PDB residues into expected PDB residues
    fp <- paste0(this.tmp.prefix,"muscle.",pdbname,".in.fa")
    write.fasta(c(EXP=pdb.exp2, OBS=pdb.obs2), fp)
    ##pdb.oea <- t(muscle(fp, paste0(this.tmp.prefix,"muscle.",pdbname,".out.fa"), as.matrix=TRUE, fix.terminals="AA"))
    pdb.oea <- t(muscle(fp, paste0(this.tmp.prefix,"muscle.",pdbname,".out.fa"), as.matrix=TRUE, fix.terminals="AA"))
    wpdb1 <- which(pdb.oea[,2]!="-")[1]
    wares1 <- which(ares[,6]==chain)[1]
    if (ares[wares1,4] != pdb.oea[wpdb1,2]) message("WARNING: congruent PDB chains don't start with same AA!")  # sanity check
    pdb.num <- ares[wares1,5]-wpdb1+1  # pdb residue numbering start number of first residue in chain
    pdb.oea <- data.frame(pdb.oea, RESNUM=pdb.num:(pdb.num+pdb.exp.len-1))  # pdb residue numbering for expected and observed sequences
    
    ## Collect CONECTs
    connect <- lapply(pdb[which(pdb[,1]=="CONECT"),2], function(x) sapply(1:(nchar(x)/5), function(y) as.numeric(substr(x,(y-1)*5+1,y*5)) ) )
    
    ##LINK         OD2 ASP B 559                MG    MG B1301     1555   1555  2.10 6345 
    ##LINK        MG    MG A1303                 O   HOH F 213     1555   1555  2.65  
    ##LINK         OH  TYR A 821                 P    DT F   9     1555   1555  1.62 2526
    
    ##SSBOND   1 CYS A    3    CYS A   16                          1555   1555  2.04  
    ##SSBOND   2 CYS A    5    CYS A   14                          1555   1555  2.04  
    ##SSBOND   3 CYS A    7    CYS A   12                          1555   1555  2.04
    
    ## Select reference transcript (if none matching any org, use whole multiple alignment)
    ## Align chain with reference; identify given variant positions in chain
    wo <- which(colnames(tmp$aligns)==gsub(" ","_",this.trans.pdb$Organism[p]))
    
#################### FIXME: DON'T KNOW IF BRANCHES OF DECISION TREE BELOW PRODUCE EQUIVALENT PRODUCTS
    
    if (length(wo)>0) {
        
        ## PDB chain org is one of the ortholog organisms; take that ortholog as 'ref' and align only to this.
        pa.ref <- pa.ref2 <- tmp$aligns[,wo]
        ref.gaps <- find.runs((pa.ref=="-")+0)
        ref.gaps <- ref.gaps[listLengths(ref.gaps)>0]
        ref.gaps <- ref.gaps[names(ref.gaps)==1]
        if (!is.na(multialn.pdb)) {
            pa.aln <- pa.aln2 <- do.call(cbind, lapply(read.fasta(multialn.pdb), function(x) unlist(strsplit(x,"")) ))
            f4 <- multialn.pdb
        } else {
            ## Align using muscle -- alignments with pairwiseAlignment() are garbage
            f3 <- paste0(this.tmp.prefix,"muscle.",pdbname,"-MA.in.fa")
            f4 <- paste0(this.tmp.prefix,"muscle.",pdbname,"-MA.out.fa")
            write.fasta(c(REF=paste(pa.ref[pa.ref!="-"],collapse=""), PDB=pdb.exp2), f3)
            ##pa.aln <- pa.aln2 <- t(muscle(f3, f4, as.matrix=TRUE, fix.terminals="AA"))
            pa.aln <- pa.aln2 <- t(muscle(f3, f4, as.matrix=TRUE, fix.terminals="AA"))
        }
        ## gap testing            pa.aln2 <- rbind(pa.aln[1:999,], matrix(c(".","X"),nrow=1), pa.aln[1000:1429,], matrix(c(".","X"),nrow=1), pa.aln[1430:nrow(pa.aln),])
        pa.aln2[pa.aln2[,1]=="-",1] <- "."  # ref gaps in PDB alignment become "." not "-"
        ins <- which(pa.aln2[,1]==".")   # positions where PDB chain has an insertion relative to reference (BEFORE adding reference MA gaps)
        for (g in 1:length(ref.gaps)) pa.aln2 <- rbind(pa.aln2[1:(ref.gaps[[g]][1]-1),], matrix("-",diff(ref.gaps[[g]])+1,2), pa.aln2[ref.gaps[[g]][1]:nrow(pa.aln2),])
        ins2 <- which(pa.aln2[,1]==".")  # positions where PDB chain has an insertion relative to reference (AFTER adding reference MA gaps)
        if (length(ins2)>0) {
            for (n in 1:length(ins2)) pa.ref2 <- c(pa.ref2[1:(ins2[n]-1)], ".", pa.ref2[ins2[n]:length(pa.ref2)])
        }
        ref.pdb.aln <- pa.aln[,match(qw(PDB,REF),colnames(pa.aln))]
        
    } else {
        
        ## PDB chain org is NOT one of the ortholog organisms; use whole multiple alignment to position PDB chain.  Take 'ref' as the source-org transcript.
        pa.ref <- pa.ref2 <- tmp$aligns[,colnames(tmp$aligns)==names(orgs)[1]]
        if (!is.na(multialn.pdb)) {
            pa.aln <- pa.aln2 <- do.call(cbind, lapply(read.fasta(multialn.pdb), function(x) unlist(strsplit(x,"")) ))
            f4 <- multialn.pdb
        } else {
            f3  <- paste0(this.tmp.prefix,"muscle.",pdbname,"-MA.in.fa")
            f4  <- paste0(this.tmp.prefix,"muscle.",pdbname,"-MA.out.fa")
            f4f <- paste0(this.tmp.prefix,"muscle.",pdbname,"-MA.out.flat.txt")
            okseq <- tmp$ortho$peptide[tmp$ortho$acceptable]
            write.fasta(c(okseq,PDB=pdb.exp2), f3)
            ##pa.aln <- pa.aln2 <- t(muscle(f3, f4, f4f, "-maxiters 64", as.matrix=TRUE, fix.terminals="AA"))
            pa.aln <- pa.aln2 <- t(muscle(f3, f4, f4f, "-maxiters 64", as.matrix=TRUE, fix.terminals="AA"))
        }
        colnames(pa.aln) <- colnames(pa.aln2) <- gsub(" ","",colnames(pa.aln))
        wp <- colnames(pa.aln)=="PDB"
        for (r in 1:nrow(pa.aln)) {
            if (all(pa.aln[r,!wp]=="-") & pa.aln[r,wp]!="-") pa.aln[r,wp] <- "."  # gaps unique to PDB chain
        }
        ins <- which(pa.aln2[,wp]==".")   # positions where PDB chain has an insertion relative to rest of MA
        ref.pdb.aln <- pa.aln[,match(c("PDB",names(orgs)[1]),colnames(pa.aln))]
        
    }
    ok <- which(tmp$vars$RefNT.T != "")
    mark.pos <- paste(tmp$vars$MAPos[ok],collapse=",")
    mark.alt <- paste(tmp$vars$AltAA[ok],collapse=",")
    this.tmp.prefix.f4 <- sub("out.fa$","",f4)
    system(paste("/home/apa/local/bin/colorMultialign -f",f4,"-o",paste0(this.tmp.prefix.f4,"cons.png"),"-c conservation --as-cons -mp",mark.pos,"-ma",mark.alt))
    
    ## must do double-gapping with MA position vectors, not just residues -- IS THIS DONE ALREADY?
    
    ## Make coord-conversion table for raw PDB chain -> chain / ref alignment -> ortholog MA
    chain2MA <- rep(NA,length(pdb.exp))  # one ortholog MA position for each PDB chain residue
    pdb.pos <- pa.pos <- c2m.pos <- 1
    for (ma.pos in 1:length(pa.ref)) {  # ortholog MA for reference transcript (same sequence as PDB-alignment reference transcript, ref.pdb.aln[,2], but gaps probably vary)
        if (pa.ref[ma.pos]=="-") {
            ## ortho MA position is a gap; skip
            message(ma.pos," Gap skip")
        } else {
            ## ortho MA position is a residue, continue
            if (ref.pdb.aln[pa.pos,2]=="-") {
                ## pdb reference pos is a gap; skip forward until residue corresponding to ortho MA position is found
                while (ref.pdb.aln[pa.pos,2]!=pa.ref[ma.pos]) pa.pos <- pa.pos + 1
                message(ma.pos," Skipped forward to ",pa.pos)
            }
            ## at this point, ref.pdb.aln reference residue matches ortho MA residue; neither are gaps
            if (ref.pdb.aln[pa.pos,1]=="-") {
                ## chain half of PDB alignment is a gap; skip
                message(ma.pos," PDB null")
            } else if (ref.pdb.aln[pa.pos,1]==".") {
                ## chain half of PDB alignment is a gap introduced by the PDB chain; do not advance pa.pos
                message(ma.pos," PDB insert")
                chain2MA[c2m.pos] <- NA
                pdb.pos <- pdb.pos + 1
                c2m.pos <- c2m.pos + 1
                next
            } else {
                ## chain half of PDB alignment is a residue
                if (ref.pdb.aln[pa.pos,1]==pdb.exp[pdb.pos]) {
                    ## PDB alignment residue agrees with pending PDB residue
                    chain2MA[c2m.pos] <- ma.pos
                    message(paste(ma.pos,"PDB positions congruent: ",pa.pos,pdb.pos,ref.pdb.aln[pa.pos,1],pdb.exp[pdb.pos]))
                } else {
                    ## residues do not agree...
                    message(paste(ma.pos,"WARNING: PDB positions not congruent! ",pa.pos,pdb.pos,ref.pdb.aln[pa.pos,1],pdb.exp[pdb.pos]))
                }
                pdb.pos <- pdb.pos + 1
                c2m.pos <- c2m.pos + 1
            }
            pa.pos <- pa.pos + 1
        }
    }
    ## QC
    ## x=rbind(pdb.exp, pa.ref[chain2MA])
    ## x
    ## sum(x[1,]==x[2,])/ncol(x)  # percent identity between pdb chain and reference ortholog
    
    ## Extract links/ssbonds/connections/proximal non-self atom distances for variant positions
    tmp$vars2 <- cbind(tmp$vars2, PDB.File="", PDB.Pos=0, PDB.Res="", PDB.Rep=0, PDB.Int=0, PDB.Comments="")  # in future, comments will have results for autodetection of known structural interactions
    ok2 <- which(!is.na(tmp$vars2$MAPos))
    pdb.w <- match(tmp$vars2$MAPos[ok2], chain2MA)
    tmp$vars2$PDB.File[ok2] <- this.trans.pdb$ID[1]
    tmp$vars2$PDB.Pos[ok2] <- pdb.oea$RESNUM[pdb.w]
    tmp$vars2$PDB.Res[ok2] <- pdb.oea$OBS[pdb.w]
    pdb.usevar <- ok2[which(!is.na(tmp$vars2$PDB.Pos[ok2]) & (tmp$vars2$PDB.Res[ok2] != "-"))]  # which variant residues are represented in PDB
    tmp$report[[p]] <- new.list(apply(tmp$vars2[pdb.usevar,],1,function(x) paste0(x[9],":",x[7],x[12],x[8]) ))
    
    for (j in 1:length(pdb.usevar)) {
        w.ares <- which(ares[,5] %in% tmp$vars2$PDB.Pos[pdb.usevar[j]] & ares[,6] %in% chain)  # atom serial numbers corresponding variant residue
        v <- as.matrix(adist)[w.ares,]  # distances in Angstroms to all other atoms in PDB
        v[,colnames(v) %in% rownames(v)] <- Inf   # ignore all within-residue distances
        v[,which(ares[,5] == ares[as.numeric(rownames(v)[1]),5]-1 & ares[,6] == chain)] <- Inf  # ignore all distances to Nt-neighbor residue's atoms
        v[,which(ares[,5] == ares[as.numeric(rownames(v)[1]),5]+1 & ares[,6] == chain)] <- Inf  # ignore all distances to Ct-neighbor residue's atoms
        v <- v[which(!(ares[w.ares,3] %in% qw(N,CA,C,O))),]  # remove backbone atoms; should be no useful structural interactions with these
        if (any(v<=angstroms[1])) {
            ## some pair of non-neighboring atoms are within reporting distance of each other
            nearby <- named.list(lapply(1:nrow(v), function(y){ z=v[y,]<=angstroms[1]; cbind(DIST=v[y,z],SERIAL=which(z),ares[match(which(z),ares[,2]),3:6]) }), rownames(v))
            nearby <- lapply(nearby, function(y){ z=y[!(y[,3] %in% qw(N,CA,C,O)),]; z[order(z[,1]),] })  # remove distances to other residues' backbone atoms; sort by distance
            nearby <- nearby[sapply(nearby,nrow)>0]
            tmp$report[[p]][[j]] <- do.call(rbind, lapply(1:length(nearby), function(n) suppressWarnings(cbind(ares[as.numeric(names(nearby)[n]),2:6],nearby[[n]])) ))  # report for reportable distances for this variant residue
            
            tmp$vars2$PDB.Rep[[pdb.usevar[j]]] <- luniq(unlist(slice.list(nearby,2)))
            if (any(unlist(slice.list(nearby,1))<=angstroms[2])) {
                ## some pair of non-neighboring non-backbone atoms are within interacting distance of each other
                tmp$vars2$PDB.Int[[pdb.usevar[j]]] <- luniq(unlist(sapply(nearby, function(x) x[x[,1]<=angstroms[2],2] )))
            }
        }
    }
    
    tmp
}


