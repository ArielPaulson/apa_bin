#!/usr/bin/env perl


sub extract_variants_new {
    
    ## Extracts snpEff annotation data from a VCF INFO field.
    ## Expects newer snpEff versions, which use the ANN (not EFF) tag.
    
    ## http://snpeff.sourceforge.net/SnpEff_manual.html#input
    ## http://varnomen.hgvs.org/bg-material/numbering/
    
    my @fields = @{ $_[0] };
    
    #;ANN=A|missense_variant|MODERATE|hax1|ENSDARG00000036764|transcript|ENSDART00000167573.1|protein_coding|1/6|c.37G>A|p.Val13Ile|186/1370|37/939|13/312||,A|downstream_gene_variant|MODIFIER|ubap2l|ENSDARG00000063219|transcript|ENSDART00000151799.2|protein_coding||c.*2811G>A|||||2109|,A|downstream_gene_variant|MODIFIER|ubap2l|ENSDARG00000063219|transcript|ENSDART00000092112.7|protein_coding||n.*2811G>A|||||1251|,A|downstream_gene_variant|MODIFIER|ubap2l|ENSDARG00000063219|transcript|ENSDART00000092138.7|protein_coding||n.*2811G>A|||||1251|,A|downstream_gene_variant|MODIFIER|ubap2l|ENSDARG00000063219|transcript|ENSDART00000151400.2|protein_coding||c.*2811G>A|||||1258|,A|non_coding_transcript_exon_variant|MODIFIER|hax1|ENSDARG00000036764|transcript|ENSDART00000053380.6|protein_coding|1/7|n.37G>A||||||;
    
    #print $DET join("\t", qw/ TranscriptID NucPos NucRef NucAlt PepPos PepRef PepAlt TransBio GeneBio Effect /), "\n";
    
    my ($vchr, $vpos, $vid, $vref, $valt) = @fields[0..4];
    my $vtag = join(":", $vchr, $vpos, $vref, $valt);
    
    my ($ann) = ($fields[7] =~ /ANN=([^;]+)/);
    foreach my $transeff (split /,/, $ann) {
        $delins1++ if $transeff =~ /delins/;
        ## e.g.: A|missense_variant|MODERATE|hax1|ENSDARG00000036764|transcript|ENSDART00000167573.1|protein_coding|1/6|c.37G>A|p.Val13Ile|186/1370|37/939|13/312||
        my ($effect, $tid, $tbio, $nmut, $pmut, $nloc, $ploc) = (split /\|/, $transeff)[(1,6,7,9,10,11,13)];
        next if $effect =~ /(intergenic|upstream|downstream|intron)/;
        next if $effect eq 'synonymous_variant' || $effect eq 'stop_retained_variant';  # synonymous
        #next if $nmut =~ /\d\[+-]/;  # why?
        if ($nmut =~ /\d\[+-]/) {
            print "$nmut: next\n";
            next;
        }
        $tid =~ s/\.\d+$//;  # drop transcript version number; not necessary since all ops are tied to $anno version.
        my $strand = $tdat{$tid}{Strand};
        my $gsbio = $tdat{$tid}{GeneSimpleBiotype};
        die "$0: Transcript '$tid' not found in geno/anno datasets!  Wrong versions?\n  $ann\n  $transeff\n\n" unless exists $tdat{$tid};
        my ($ntype, $npos, $nref, $nalt) = &parse_HGVS_n($tid, $nmut, $vref, $nloc);
        my ($ptype, $ppos, $pref, $palt);
        ($ptype, $ppos, $pref, $palt) = &parse_HGVS_p($tid, $pmut, $ntype) if $pmut;
        $delins2++ if $ptype eq 'delins';
        print $TAB join("\t", $tid, $vtag, $nmut, $ntype, $npos, $nref, $nalt, $pmut, $ptype, $ppos, $pref, $palt, $tbio, $gsbio, $effect, $transeff), "\n";
    }
    
}



sub parse_HGVS_n {
    
    ## Input: transcript ID, HGVS string, VCF ref string, and snpEff nuc coords (if any):
    ## Output: variant type, transcript coords, ref string, and alt string.
    
    my ($TID, $NMUT, $VREF, $NLOC) = @_;
    
    $NLOC =~ s/\/.*//;  # keep only first number
    (my $VREFRC = $VREF) =~ tr/ACGTacgt/TGCAtgca/;
    
    my $strand = $tdat{$TID}{strand};
    my $tcdss = $tdat{$TID}{Trans_CDS_Start};
    my $tcdse = $tdat{$TID}{Trans_CDS_End};
    
    my ($htype, $hdir1, $hdir2, $nstart, $nend, $npos, $nseq1, $nref, $nalt);
    my $ntype = 'snp';   # unless replaced below

    ## Formula matching
    if ($NMUT =~ /([gncm])\.([\+\*-])?(\d+)_([\+\*-])?(\d+)(ins|del|dup)([A-Za-z]+)?/) {
        ($htype, $hdir1, $nstart, $hdir2, $nend, $ntype, $nseq1) = ($1, $2, $3, $4, $5, $6, $7);
    } elsif ($NMUT =~ /([gncm])\.([\+\*-])?(\d+)(ins|del|dup)([A-Za-z]+)?/) {
        ($htype, $hdir1, $npos, $ntype, $nseq1) = ($1, $2, $3, $4, $5);
    } elsif ($NMUT =~ /([gncm])\.([\+\*-])?(\d+)([A-Z])>([A-Z])/) {
        ($htype, $hdir1, $npos, $nref, $nalt) = ($1, $2, $3, $4, $5);
    } else {
        print "WARNING: failed to parse HGVS string '$NMUT'!\n";
    }
    
    ## Position correction
    if ($htype eq 'n' || $htype eq 'c') {
        if ($nstart) {
            ($nstart, $nend) = ($NLOC, $NLOC+$nend-$nstart) if $NLOC;
            if (!$hdir1 && !$hdir2) {
                $npos = "$nstart-$nend";
            } elsif ($hdir1 eq '-' && $hdir2 eq '-') {
                $npos = ($tcdss-$nstart).'-'.($tcdss-$nend);
            } elsif ($hdir1 eq '*' && $hdir2 eq '*') {
                $npos = ($tcdse+$nstart).'-'.($tcdse+$nend);
            } elsif ($hdir1 eq '+' && $hdir2 eq '+') {
                print "WARNING: introns excluded; we should not be seeing '+' modifiers! ($NMUT)\n";
            } else {
                print "WARNING: mixed modifiers! ($NMUT)\n";
            }
        } elsif ($npos) {
            $npos = $NLOC if $NLOC;
            if (!$hdir1) {
                ## do nothing;
            } elsif ($hdir1 eq '-') {
                $npos = $tcdss-$npos;
            } elsif ($hdir1 eq '*') {
                $npos = $tcdse+$npos;
            } elsif ($hdir1 eq '+') {
                print "WARNING: introns excluded; we should not be seeing '+' modifiers! ($NMUT)\n";
            } else {
                print "WARNING: mixed modifiers! ($NMUT)\n";
            }
        } else {
            print "WARNING: position parsing failure! ($NMUT)\n"
        }
    } else {
        print "WARNING: no handling for HGVS $htype-type coords! ($NMUT)\n"
    }
    
    ## Ref/Alt correction
    if ($ntype eq 'ins') {
        $nref = $strand eq '-' ? $VREFRC : $VREF;
        $nalt = $strand eq '-' ? "$nseq1$VREFRC" : "$VREF$nseq1";
    } elsif ($ntype eq 'del') {
        $nref = $nseq1;
        $nalt = '';
    } elsif ($ntype eq 'dup') {
        $nref = $nseq1;
        $nalt = "$nseq1$nseq1";
    } elsif ($ntype eq 'snp') {
        ## do nothing
    } else {
        print "WARNING: vartype parsing failure! ($NMUT)\n"
    }
    
    return ($ntype, $npos, $nref, $nalt);
    
    
    
    ## OLD CODE ##
    if (0) {
        my ($nmut, $nalt1, $tid, $vref, $nloc, $ann, $transeff);  # prevent use-strict fails
        if ($nmut =~ /(ins|del|dup)/) {
            if ($nmut =~ /([gncm])\.([\+\*-])?(\d+)_([\+\*-])?(\d+)(ins|del|dup)([A-Za-z]+)/) {
                ## multi-bp indel
                ($htype, $hdir1, $nstart, $hdir2, $nend, $ntype, $nalt1) = ($1, $2, $3, $4);
                $npos = "$nstart-$nend";
            } elsif ($nmut =~ /([gncm])\.([\+\*-])?(\d+)(ins|del|dup)([A-Za-z]+)/) {
                ## single-bp indel
                ($htype, $hdir1, $npos, $ntype, $nalt1) = ($1, $2, $3);
            } else {
                ## ??
                print "\nCannot interpret nt indel (1): $nmut\n\n";
            }
            if ($ntype eq 'ins') {
                $nref = $strand eq '-' ? $VREFRC : $VREF;
                $nalt = $strand eq '-' ? "$nalt1$VREFRC" : "$VREF$nalt1";
            } elsif ($ntype eq 'del') {
                $nref = $nalt1;
                $nalt = '';
            } elsif ($ntype eq 'dup') {
                $nref = $nalt1;
                $nalt = "$nalt1$nalt1";
            } else {
                print "\nCannot interpret nt indel (2): $nmut\n\n";
            }
        } else {
            ($htype, $hdir1, $npos, $nref, $nalt) = ($nmut =~ /([gncm])\.([\+\*-])?(\d+)([A-Z])>([A-Z])/);
        }
        if ($nmut && !$npos) {
            print "Cannot parse variant '$nmut':\n  $ann\n  $transeff\n\n";
        }
        if ($nloc) {
            ## if $nloc present, then it is correct and $npos is wrong
            if ($npos =~ /-/) {
                my ($start, $end) = split /-/, $npos;
                $npos = join('-', $nloc, $nloc+$end-$start);
            } else {
                $npos = $nloc;
            }
        } elsif ($htype eq 'c') {
            if ($nstart) {
                print "C-TYPE (1): $htype, $hdir1, $nstart, $hdir2, $nend, $ntype, $nalt1\n";
                ## not yet sure what to do here, need to see one
                $npos = undef;
            } elsif ($npos) {
                print "C-TYPE (2): $htype, $hdir1, $npos, $nref, $nalt\n";
                $npos = $hdir1 eq '-' ? $tcdss-$npos : $hdir1 eq '*' ? $tcdse+$npos : undef;
                
            }
            print "WARNING: Unparseable c-type HGVS annot '$nmut' !\n" unless $npos;
        }
    }
    
}



sub parse_HGVS_p {
    
    ## Input: transcript ID, HGVS string, and nucleotide variant type:
    ## Output: variant type, peptide coords, ref string, and alt string.
    
    my ($TID, $PMUT, $NTYPE) = @_;
    
    my ($pstart, $pend, $ppos, $pseq1, $aa1, $aa2, $pref, $palt);
    my $ptype = $NTYPE eq 'snp' ? 'sap' : '';   # unless replaced below
    
    ## Formula matching and variant-type modification
    if ($PMUT =~ /p\.(\w{3})(\d+)_(\w{3})(\d+)(delins|ins|del|dup)([A-Za-z]+)?/) {
        ($aa1, $pstart, $aa2, $pend, $ptype, $pseq1) = ($AAs{$1}, $2, $AAs{$3}, $4, $5, $6);
    } elsif ($PMUT =~ /p\.(\w{3})(\d+)(delins|ins|del|dup)([A-Za-z]+)?/) {
        ($pref, $ppos, $ptype, $pseq1) = ($AAs{$1}, $2, $3, $4);
    } elsif ($PMUT =~ /p\.(\w{3})(\d+)fs/) {
        ($pref, $ppos, $palt) = ($AAs{$1}, $2, '');
        $ptype = $ppos == 1 ? 'trunc' : 'shift';  # $ppos==1 means start loss; others, frameshift
    } elsif ($PMUT =~ /p\.(\w{3})(\d+)(\w{3})ext\*/) {
        ($pref, $ppos, $palt) = ($AAs{$1}, $2, "$AAs{$3}+");
        $ptype = 'ext';
    } elsif ($PMUT =~ /p\.(\w{3})(\d+)(\w{3}|\*)/) {
        ($pref, $ppos, $palt) = ($AAs{$1}, $2, $AAs{$3});
        $ptype = 'trunc' if $palt eq '*';
    } else {
        print "WARNING: failed to parse HGVS string '$PMUT'!\n";
    }
    
    ## Position correction
    if ($pstart) {
        $ppos = "$pstart-$pend";
    } elsif ($ppos) {
        ## do nothing;
    } else {
        print "WARNING: position parsing failure! ($PMUT)\n"
    }
    
    ## Ref/Alt correction
    if (exists $iddtypes{$ptype}) {  # ins, del, dup, delins
        if ($pstart) {
            $pref = $pstart+1==$pend ? "$aa1$aa2" : $aa1.('.'x($pend-$pstart-1)).$aa2;  # ref is flanking AAs; must likewise flank the alt sequence
        }
        if ($ptype =~ 'ins') {  # ins, delins
            $palt .= $AAs{$1} while $pseq1 =~ /(.{3})/g;
            #print "\nConverted: '$pseq1' -> '$palt'\n\n";
            $palt = "$aa1$palt$aa2";
        } elsif ($ptype eq 'dup') {
            $palt = $pref =~ /-/ ? "$pref+$pref" : "$pref$pref";  # one outcome for deletions, two for duplications
        } elsif ($ptype eq 'del') {
            $palt = '';
        } else {
            print "WARNING: vartype parsing failure (1)! ($PMUT)\n"
        }
    } elsif ($ptype =~ /^(sap|ext|trunc|shift)$/) {
        ## do nothing
    } else {
        print "WARNING: vartype parsing failure! ($PMUT)\n"
    }
    if ($pref =~ /[A-Z][a-z]{2}/) {
        print "\nFailed to convert ref AA '$pref'! ($PMUT)\n\n";
    }
    if ($palt =~ /[A-Z][a-z]{2}/) {
        print "\nFailed to convert alt AA '$palt'! ($PMUT)\n\n";
    }
    
    return ($ptype, $ppos, $pref, $palt);
    
    
    
    ## OLD CODE ##
    if (0) {
        my ($pmut, $palt1, $tid, $vref, $ploc, $ann, $transeff, $ntype);  # prevent use-strict fails
        if ($pmut =~ /(ins|del|dup)/) {
            if ($pmut =~ /p\.(\w{3})(\d+)_(\w{3})(\d+)(delins|ins|del|dup)([A-Za-z]+)?/) {
                ($aa1, $pstart, $aa2, $pend, $ptype, $palt1) = ($AAs{$1}, $2, $AAs{$3}, $4, $5, $6);
                $ppos = "$pstart-$pend";
                $pref = $pstart+1==$pend ? "$aa1$aa2" : $aa1.('.'x($pend-$pstart-1)).$aa2;
            } elsif ($pmut =~ /p\.(\w{3})(\d+)(delins|ins|del|dup)([A-Za-z]+)?/) {
                ($pref, $ppos, $ptype, $palt1) = ($AAs{$1}, $2, $3, $4);
            } else {
                print "\nCannot interpret aa indel (1): $pmut\n\n";
            }
            $delins2++ if $ptype eq 'delins';
            if ($ptype =~ 'ins') {  # ins or delins
                $palt .= $AAs{$1} while $palt1 =~ /(.{3})/g;
                print "\nConverted: '$palt1' -> '$palt'\n\n";
                $palt = "$aa1$palt$aa2";
            } elsif ($ptype eq 'dup') {
                $palt = $pref =~ /-/ ? "$pref+$pref" : "$pref$pref";  # one outcome for deletions, two for duplications
            }
        } elsif ($pmut =~ /p\.(\w{3})(\d+)fs/) {
            ($pref, $ppos, $palt) = ($AAs{$1}, $2, '');
            $ptype = $ppos == 1 ? 'trunc' : 'shift';  # $ppos==1 means start loss; others, frameshift
        } elsif ($pmut =~ /p\.(\w{3})(\d+)(\w{3})ext\*/) {
            ($pref, $ppos, $palt) = ($AAs{$1}, $2, "$AAs{$3}+");
            $ptype = 'ext';
        } elsif ($pmut =~ /p\.(\w{3})(\d+)(\w{3}|\*)/) {
            ($pref, $ppos, $palt) = ($AAs{$1}, $2, $AAs{$3});
            $ptype = 'trunc' if $palt eq '*';
        }
        if ($pmut && !$ppos) {
            print "\nCannot parse variant '$pmut':\n  $ann\n  $transeff\n\n";
        }
        if ($pref =~ /[A-Z][a-z]{2}/) {
            print "\nFailed to convert ref AA '$pref'! ($pmut)\n\n";
        }
        if ($ppos && !$ptype) {
            if ($ntype eq 'snp') {
                $ptype = 'sap';
            } else {
                print "\nFailing to ptype '$pmut'!\n\n";
            }
        }
    }
    
}



sub create_sequences_1 {
    
    ## This will be the new version, once multi-variant sequences are a thing
    
    my $tid = shift;
    my $strand = $tdat{$tid}{Strand};
    my $tseqn = $tdat{$tid}{SEQNT};
    my $tseqp = $tdat{$tid}{SEQAA};
    
    foreach my $i (sort {$a <=> $b} keys %{ $vdat{$tid} }) {
        my ($vtag, $nmut, $pmut, $effect) = @{ $vdat{$tid}{$i}{V} };
        next if $already{$tid}{$vtag};  # snpEff can output multiple effects for the SAME transcript/variant pair; e.g. see ENSDART00000130003 for cbio.heb.102 migr_F167 @ chr20:204921, A->G
        $already{$tid}{$vtag} = 1;
        
        ## NUCLEOTIDE
        
        my ($ntype, $npos, $nref, $nalt) = @{ $vdat{$tid}{$i}{N} };
        my @npos = split /-/, $npos;
        my $nrlen = length($nref);
        my $nalen = length($nalt);
        my $nrseq = substr($tseqn, $npos-1, $nrlen);
        my $nrseq2 = substr($tseqn, $npos-3, $nrlen+4);
        print "$tid $strand $vtag N: $ntype, $npos: $nref->$nalt";
        $total++;
        if ($nrseq ne $nref) {
            $bad++;
            print ' BAD';
        } else {
            $good++;
            print '    ';
        }
        print " ($nrseq, $nrseq2), $nmut\n";
        
        if ($ntype eq 'snp') {
            
        } elsif ($ntype eq 'ins') {
            
        } elsif ($ntype eq 'del') {
            
        } elsif ($ntype eq 'dup') {
            
        } else {
            die "$0: Don't know what to do with nucleotide variant type '$ntype'!\n";
        }
        
        ## PEPTIDE
        
        my ($ptype, $ppos, $pref, $palt) = @{ $vdat{$tid}{$i}{P} };
        next unless $ptype;
        my @ppos = split /-/, $ppos;
        my $prlen = length($pref);
        my $palen = length($palt);
        my $prseq = substr($tseqp, $ppos-1, $prlen);
        my $prseq2 = substr($tseqp, $ppos-2, $prlen+4);
        print "$tid $strand $vtag P: $ptype, $ppos: $pref->$palt";
        $prseq ne $pref ? print ' BAD' : print '    ';
        print " ($prseq, $prseq2), $pmut\n";
        
        if ($ptype eq 'sap') {
            
        } elsif ($ptype eq 'ins') {
            
        } elsif ($ptype eq 'del') {
            
        } elsif ($ptype eq 'dup') {
            
        } elsif ($ptype eq 'shift') {
            
        } elsif ($ptype eq 'trunc') {
            
        } elsif ($ptype eq 'ext') {
            
        } else {
            die "$0: Don't know what to do with peptide variant type '$ptype'!\n";
        }
        
    }
    
}



sub create_sequences_M {

    ## Will eventually produce multi-variant sequences.
    
    my $tid = shift;
    my $tseqn = $tdat{$tid}{SEQNT};
    my $tseqa = $tdat{$tid}{SEQAA};
    
    ## NUCLEOTIDE
    
    ## First, test change positions for overlaps
    my (%npos, %ppos);
    foreach my $pos (sort {$b <=> $a} keys %{ $vdat{$tid}{N} }) {
        
        
    }
    
    ## Second, create the single-variant set
    my (%npos, %ppos);
    foreach my $pos (sort {$b <=> $a} keys %{ $vdat{$tid}{N} }) {
        
        
    }
    
    ## Third, expand single-variant sequences to multi-variant IF no variants overlap
    ## NOT: can't do this for indel-incorporated sequences above
    my (%npos, %ppos);
    foreach my $pos (sort {$b <=> $a} keys %{ $vdat{$tid}{N} }) {   ## ALWAYS WORK BACKWARDS THROUGH SEQUENCE -- in case of indels
        
        
    }
    
    ## PEPTIDE 
    
    
}


