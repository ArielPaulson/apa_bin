

get.orthologs <- function(i) {
    
    ## Given a transcript, and identity cutoff, and a list of organisms, construct a table of Ensembl ortholog data
    
    message("\n",genes[i,1],": processing orthologs...")
    ortho.tmp <- new.list(names(orgs))
    
    for (j in 1:O) {
        ## Select which symbol attribute to use, depending on Mart (basically human = "hgnc_symbol", nonhuman = first "*_symbol" attrib if any, otherwise "wikigene_name")
        symb <- grep("_symbol$",la[[j]][,1],value=TRUE)
        symb <- ifelse (names(orgs)[j]=="Homo_sapiens","hgnc_symbol",setdiff(symb,"hgnc_symbol")[1])
        if (is.na(symb)) symb <- "wikigene_name"
        if (!symb %in% la[[j]][,1]) symb <- "external_gene_name"
        ## Acquire (ortho|para|homo)log genes, if any
        g.ortho <- getBM(attributes=grep("log_",grep(orgs.ens[j],la[[1]][,1],value=T),value=T), filters="ensembl_transcript_id", values=genes$ensembl_transcript_id[i], mart=marts[[1]])  # "log_" targeting "(ortho|para|homo)log_" attributes
        ##message("\n",orgs.ens[j]); print(g.ortho); message("")
        message(paste0(sum(!is.na(g.ortho[[1]]))," ",orgs.ens[j],"_gene_ensembl"))
        if (j>1 & all(is.na(g.ortho[[1]]))) next  # skip non-source orgs with no orthologs
        ## strip sci name from these col names, else future rbinds will have issues
        colnames(g.ortho) <- sub(paste0(orgs.ens[j],"_"),"",colnames(g.ortho))
        ## if source org, use gene for query transcript (to eliminate dealing with inparalogs), otherwise use resulting ortholog
        gids <- ternary(j==1, genes$ensembl_gene_id[i], g.ortho[,1])
        ## get gene data for ortholog(s)
        g.data <- getBM(attributes=c("ensembl_gene_id",symb,qw(chromosome_name,strand,start_position,end_position,gene_biotype,description)), filters="ensembl_gene_id", values=gids, mart=marts[[j]])   # once also included: source, status
        ## standardize this col name, else future rbinds will have issues
        colnames(g.data)[2] <- "symbol"
        ## find all peptides for ortholog(s)
        peps <- getBM(qw(ensembl_peptide_id,ensembl_transcript_id,ensembl_gene_id), "ensembl_gene_id", gids, marts[[j]])
        ## remove transcripts with no peptide (=noncoding)
        peps <- peps[peps[,1]!="",]
        ## if source org gene has > 1 peptide, select the one corresponding to the query transcript
        if (j==1) peps <- peps[peps[,2]==genes$ensembl_transcript_id[i],]
        ## get sequence for each 'peps' peptide
        p.seq <- getSequence(id=peps[,1], type="ensembl_peptide_id", seqType="peptide", mart=marts[[j]])
        ## add sequences to 'peps'
        peps <- cbind(peps, aln_score=0, idt_pct=0, acceptable=FALSE, p.seq[match(peps[,1],p.seq[,2]),1,drop=FALSE])
        
        ## Screen peptide(s) for similarity to query
        
        ## generate multiple alignment scores; tseq$aa[[i]] is the reference
        pa <- sapply(1:nrow(peps), function(g) pairwiseAlignment(tseq$aa[[i]], peps$peptide[g], substitutionMatrix=BLOSUM100, gapOpening=0, gapExtension=-5) )
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
        ortho.tmp[[j]] <- df
    }
    ortho.tmp
}



align.orthologs <- function() {
    
    ## this.trans, ortho, vars, tseq2, genes, multialn.ortho, codons, 
    
    ## Produce a multiple alignment for a set of orthologous transcript sequences
    
    tmp <- new.list(qw(ortho,aligns,exons,vars,vars2))
    tmp$ortho <- this.ortho  # will upgrade and return
    tmp$vars <- this.vars    # will upgrade and return
    
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
    f0 <- paste0(prefix,".ortho.all.fa")
    write.fasta(allseq, f0)
    
    okseq <- allseq[tmp$ortho$acceptable]
    names(okseq) <- rownames(tmp$ortho)[tmp$ortho$acceptable]
    O2 <- length(okseq)  # may have lost some orgs
    message("Good orthologs: ",paste(setdiff(names(okseq),names(orgs)[1]),collapse=" "))

    if (O2<2) {
        
        message("Too few sequences to proceed...\n")
        return()
        
    } else {
        
        message(this.trans,": aligning orthologs...")
        if (!is.na(multialn.orth)) {
            tmp$aligns <- do.call(cbind, lapply(read.fasta(multialn.orth), function(x) unlist(strsplit(x,"")) ))
            f2 <- multialn.orth
        } else {
            f1 <- paste0(prefix,".ortho.muscle_in.fa")
            f2 <- paste0(prefix,".ortho.muscle_out.fa")
            save.image(RData)
            write.fasta(okseq, f1)
            ##tmp$aligns <- t(muscle(f1, f2, as.matrix=TRUE, fix.terminals="AA"))
            tmp$aligns <- t(muscle(f1, f2, as.matrix=TRUE, fix.terminals="AA"))
        }
        mark.row <- which(colnames(tmp$aligns)==names(orgs)[1])
        
        ## map exon positions for transcripts.  initial: matrices for (-) strand genes will have to be row-reversed
        message("Mapping exon structure...")
        x <- getBM(attributes=qw(ensembl_exon_id,exon_chrom_start,exon_chrom_end,cds_start,cds_end,cds_length), filters="ensembl_transcript_id", values=this.trans, mart=marts[[1]])
        x <- cbind(x[,1,drop=FALSE],x[,2:3],exon_length=x[,3]-x[,2]+1,x[,4:5],cds_length=x[,5]-x[,4]+1)
        tmp$exons <- x[order(x[,2]),]  # sort by start pos
        
        ## Convert genomic position to AA position in source transcript
        ## Then map position in source transcript to MA position -> mark.pos

        if (tmp$ortho$strand[1] == -1) tmp$exons <- tmp$exons[nrow(tmp$exons):1,]
        
        for (r in 1:nrow(tmp$vars)) {
            x <- cbind(tmp$vars$Pos[r]-tmp$exons[,2], tmp$exons[,3]-tmp$vars$Pos[r])
            wx <- which(x[,1]>0 & x[,2]>0)
            message(paste(i,r,wx,tmp$ortho$strand[1] == -1,nrow(x)))
            if (length(wx)==0) {
                message("Variant at position ",tmp$vars$Pos[r]," is not within boundaries of transcript ",this.trans)
                next
            }
            if (tmp$ortho$strand[1] == -1) {
                tpos.nt <- sum(tmp$exons[1:wx,4])-x[wx,1]
                tmp$vars$AltNT.T[r] <- translate(tmp$vars$AltNT.G[r],qw(A,C,G,T,a,c,g,t),qw(T,G,C,A,t,g,c,a))
            } else {
                tpos.nt <- sum(tmp$exons[1:wx,4])-x[wx,2]
                tmp$vars$AltNT.T[r] <- tmp$vars$AltNT.G[r]
            }
            tmp$vars$RefNT.T[r] <- tseq2$nt[[this.trans]][tpos.nt]
            w1 <- which(tmp$exons[,5]==1)
            utr5 <- sum(tmp$exons[1:w1,4])-tmp$exons[w1,7]
            tpos.aa <- ceiling((tpos.nt-utr5)/3)
            ma.pos <- which(tmp$aligns[,mark.row]!="-")[tpos.aa]
            nt.phase <- (tpos.nt%%3)-1
            message(paste(tpos.nt,nt.phase,tpos.aa,ma.pos))
            if (nt.phase==-1) nt.phase <- 2
            tmp$vars$TransBp[r] <- tpos.nt
            tmp$vars$Phase[r] <- nt.phase
            tmp$vars$TransAa[r] <- tpos.aa
            tmp$vars$MAPos[r] <- ma.pos
            tmp$vars$RefAA[r] <- tseq2$aa[[this.trans]][tpos.aa]
            codon.start <- tpos.nt-nt.phase
            codon.end <- codon.start+2
            codon.seq <- tseq2$nt[[this.trans]][codon.start:codon.end]
            message(paste("REF:",codon.seq,codons[match(paste(codon.seq,collapse=""),codons[,1]),2],paste0("(",tseq2$aa[[this.trans]][tpos.aa],")")))
            codon.seq[nt.phase+1] <- ifelse(tmp$ortho$strand[1] == -1, tmp$vars$AltNT.T[r], tmp$vars$AltNT.G[r])
            message(paste("ALT:",codon.seq,codons[match(paste(codon.seq,collapse=""),codons[,1]),2]))
            tmp$vars$AltAA[r] <- codons[match(paste(codon.seq,collapse=""),codons[,1]),2]
        }
        tmp$vars$RefNT.T[tmp$vars$RefNT.T==0] <- ""  # annotation errors
        ok <- which(tmp$vars$RefNT.T != "")
        mark.dat <- unique(data.frame(Pos=tmp$vars$MAPos[ok], Alt=tmp$vars$AltAA[ok]))
        mark.pos <- paste(mark.dat$Pos,collapse=",")
        mark.alt <- paste(mark.dat$Alt,collapse=",")
        
        prefix.f2 <- sub(".out.fa$","",f2)
        cmd <- paste("/home/apa/local/bin/colorMultialign -f",f2,"-o",paste0(prefix.f2,".res.png"),"-c residues --as-cons -mr",mark.row,"-mp",mark.pos,"-ma",mark.alt)
        message(cmd)
        system(cmd)  # could capture STDERR, for something
        cmd <- paste("/home/apa/local/bin/colorMultialign -f",f2,"-o",paste0(prefix.f2,".cons.png"),"-c conservation --tabular -mr",mark.row,"-mp",mark.pos,"-ma",mark.alt)
        message(cmd)
        system(cmd)
        ## gthumb(paste0(prefix,".cons.png"))
        
        ## Finally, classify marked positions in MA -- have colorMultialign produce tabular output, rows = classified positions?
        ## Include AA substitution, if not splice effect
        classify <- read.delim(paste0(prefix.f2,".cons.png.table.txt"), as.is=TRUE)
        class2vars <- classify[tmp$vars$MAPos[ok],c(O2+c(2:6),2:(O2+1))]
        
        if (nrow(class2vars)<nrow(tmp$vars)) {
            ## some variant positions were not classifiable?  error in data entry?
            for (r in (nrow(class2vars)+1):length(ok)) class2vars[r,] <- NA
        }
        tmp$vars2[[i]] <- cbind(tmp$vars, class2vars)
        
        tmp
    }
}



process.pdb <- function() {
    
    ## WHERE DID 'pdb' COME FROM ????  (now "this.pdb")
    
    ## trans.pdb, trans.pdb.dat, aligns, ortho, vars, vars2
    
    ## if PDB file specified, and PDB(s) assigned to this transcript, analyze PDB for side-chain proximity to other residues
    
    p <- 1  # CURRENTLY: can only report on 1 PDB file per transcript, even if multiple given
    chain <- this.trans.pdb$Chain[p]
    pdbname <- paste0("PDB:",this.trans.pdb$ID[p])
    pdb <- this.trans.pdb.dat[[p]]
    tmp <- new.list(qw(vars2,report))
    tmp$report <- new.list(this.trans.pdb$ID[p])
    
    ## Extract Expected AA Sequence
    pdb.exp <- this.pdb[which(this.pdb[,1]=="SEQRES")[grepl(paste0("^ +[0-9]+ ",chain),this.pdb[this.pdb[,1]=="SEQRES",2])],2]
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
    fp <- paste0(prefix,".muscle.",pdbname,".in.fa")
    write.fasta(c(EXP=pdb.exp2, OBS=pdb.obs2), fp)
    ##pdb.oea <- t(muscle(fp, paste0(prefix,".muscle.",pdbname,".out.fa"), as.matrix=TRUE, fix.terminals="AA"))
    pdb.oea <- t(muscle(fp, paste0(prefix,".muscle.",pdbname,".out.fa"), as.matrix=TRUE, fix.terminals="AA"))
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
    wo <- which(colnames(aligns)==gsub(" ","_",trans.pdb$Organism[p]))
    
#################### FIXME: DON'T KNOW IF BRANCHES OF DECISION TREE BELOW PRODUCE EQUIVALENT PRODUCTS
    
    if (length(wo)>0) {
        
        ## PDB chain org is one of the ortholog organisms; take that ortholog as 'ref' and align only to this.
        pa.ref <- pa.ref2 <- aligns[,wo]
        ref.gaps <- find.runs((pa.ref=="-")+0)
        ref.gaps <- ref.gaps[listLengths(ref.gaps)>0]
        ref.gaps <- ref.gaps[names(ref.gaps)==1]
        if (!is.na(multialn.pdb)) {
            pa.aln <- pa.aln2 <- do.call(cbind, lapply(read.fasta(multialn.pdb), function(x) unlist(strsplit(x,"")) ))
            f4 <- multialn.pdb
        } else {
            ## Align using muscle -- alignments with pairwiseAlignment() are garbage
            f3 <- paste0(prefix,".muscle.",pdbname,"-MA.in.fa")
            f4 <- paste0(prefix,".muscle.",pdbname,"-MA.out.fa")
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
        pa.ref <- pa.ref2 <- aligns[,colnames(aligns)==names(orgs)[1]]
        if (!is.na(multialn.pdb)) {
            pa.aln <- pa.aln2 <- do.call(cbind, lapply(read.fasta(multialn.pdb), function(x) unlist(strsplit(x,"")) ))
            f4 <- multialn.pdb
        } else {
            f3 <- paste0(prefix,".muscle.",pdbname,"-MA.in.fa")
            f4 <- paste0(prefix,".muscle.",pdbname,"-MA.out.fa")
            f4f <- paste0(prefix,".muscle.",pdbname,"-MA.out.flat.txt")
            okseq <- ortho$peptide[ortho$acceptable]
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
    ok <- which(vars$RefNT.T != "")
    mark.pos <- paste(vars$MAPos[ok],collapse=",")
    mark.alt <- paste(vars$AltAA[ok],collapse=",")
    prefix.f4 <- sub(".out.fa$","",f4)
    system(paste("/home/apa/local/bin/colorMultialign -f",f4,"-o",paste0(prefix.f4,".cons.png"),"-c conservation --as-cons -mp",mark.pos,"-ma",mark.alt))
    
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
    tmp$vars2 <- cbind(vars2, PDB.File="", PDB.Pos=0, PDB.Res="", PDB.Rep=0, PDB.Int=0, PDB.Comments="")  # in future, comments will have results for autodetection of known structural interactions
    ok2 <- which(!is.na(vars2$MAPos))
    pdb.w <- match(vars2$MAPos[ok2], chain2MA)
    tmp$vars2$PDB.File[ok2] <- trans.pdb$ID[1]
    tmp$vars2$PDB.Pos[ok2] <- pdb.oea$RESNUM[pdb.w]
    tmp$vars2$PDB.Res[ok2] <- pdb.oea$OBS[pdb.w]
    pdb.usevar <- ok2[which(!is.na(vars2$PDB.Pos[ok2]) & (vars2$PDB.Res[ok2] != "-"))]  # which variant residues are represented in PDB
    tmp$report[[p]] <- new.list(apply(vars2[pdb.usevar,],1,function(x) paste0(x[9],":",x[7],x[12],x[8]) ))
    
    for (j in 1:length(pdb.usevar)) {
        w.ares <- which(ares[,5] %in% vars2$PDB.Pos[pdb.usevar[j]] & ares[,6] %in% chain)  # atom serial numbers corresponding variant residue
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


