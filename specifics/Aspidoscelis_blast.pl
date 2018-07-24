#!/usr/bin/env perl

=pod

=head1 SYNOPSIS

This is an interface to run blast / parse results / output filtered results and alignment maps.  Databases 
are de novo RNAseq assemblies for various Aspidoscelis species; inputs are gene sequences from other species.

perl Aspidoscelis_blast.pl [-i input] [-s species] [-p program] [-c cores] [-o output_folder] [-mx min_subj_sim] [-mc min_subj_covg] [-mi min_subj_ident] [-ml min_align_length] --multi --nomap


=head1 OPTIONS

=over

=item B<-i <input_file>>

Input file (fasta)

=item B<-s <species>>

Species to blast against; default 'gularis'

=item B<-p <program>>

Blast program to run:

=over

=item B<blastn>

Nucleotide query -> nucleotide database

=item B<blastp>

Protein query -> protein database

=item B<blastx>

Nucleotide query -> translated nucleotide query -> protein database

=back

=item B<-c <cores>>

An integer (2 by default) specifying the number of CPUs to use for the blast.  Please check server load before allocating > 2 CPUs.

=item B<-o <output_folder>>

An output folder to write output to

=item B<-mx <decimal>>

Minimum alignment similarity index (identity% * subject-coverage%) to allow (on [0,1])

=item B<-mc <decimal>>

Minimum alignment coverage% (of subject) to allow (on [0,1])

=item B<-mi <decimal>>

Minimum alignment identity% to allow (on [0,1])

=item B<-ml <integer>>

Minimum alignment length to allow (in nt/aa)

=item B<--multi>

Allow a query to receive multiple hits from the same subject (default = takes first best alignment only)

=item B<--nomap>

Do not print coverage maps

=back


=head1 EXAMPLES

=item B<Nucleotide blast, no filters>

perl Aspidoscelis_blast.pl -i genes.fa -s gularis -p blastn

=item B<Protein blast, minimum query coverage 50%, minimum query identity 50%>

perl Aspidoscelis_blast.pl -i genes.fa -s gularis -p blastp -mc 0.5 -mi 0.5

=item B<Blastx, minimum alignment length 1kb, no map images>

perl Aspidoscelis_blast.pl -i genes.fa -s gularis -p blastx -ml 1000 --nomaps


=head1 OUTPUTS

=item B<blast_results.txt>

Raw blast output (-m 8 tabular), no headers

=item B<blast_filtered_results.txt>

Filtered blast output, with headers

=item B<filtered_results_README.txt>

A readme file explaining blast_filtered_results.txt

=item B<blast_mapdata.txt>

Alignment data used to create the blast map images

=item B<blast_maps.R>

R script which generates the blast map images

=item B<*_map.png>

The blast map images, one per query, showing alignments of Aspidoscelis contigs

=item B<blast_map_color_key.png>

The color key for the blast map images


=head1 VERSION

$Revision: 1.0$

=head1 AUTHOR

Ariel Paulson [apa@stowers.org]

=head1 DEPENDENCIES

Perl 5.8; NCBI blast; R; Getopt::Long; Pod::Usage; Cwd; 

=cut 


##############################################  BEGIN ACTUAL CODE  ##############################################


use Getopt::Long;
use Pod::Usage;
use Cwd;
use strict;

## Formatting a new nucleotide db: "formatdb -p F -a F -o T -i <db file> -l <log file>"
## Formatting a new protein db:    "formatdb -p T -a F -o T -i <db file> -l <log file>"

## add new blast programs, add or alter blast dbs / flags here
my %blastflags = (
    "blastn" => { "dbtype" => "nt", "flags" => "" },   # nucleotide
    "blastp" => { "dbtype" => "aa", "flags" => "" },   # protein
    "blastx" => { "dbtype" => "aa", "flags" => "" },   # protein db, translated nucleotide query
    "tblastn" => { "dbtype" => "nt", "flags" => "" }   # translated nucleotide db, protein query
);

# add new species here
my %knownspecies = map {($_=>1)} qw/ gularis /;
my $path = '/n/facility/Bioinformatics/analysis/Baumann/LizardAssembly/blast';

## get command line
my ($help, $man);
my ($prog, $input, $minsim, $mincov, $minidt, $minaln, $multi, $nomap);
my ($species, $cores, $outdir) = ('gularis', 2, 'results');           # defaults
GetOptions(
    "s=s" => \$species,      # Aspidoscelis species
    "p=s" => \$prog,         # blast program to run; must be in %blastflags above
    "c=i" => \$cores,        # number of CPUs to use for blast
    "i=s" => \$input,        # input fasta to blast
    "o=s" => \$outdir,       # output directory
    "mx=f" => \$minsim,      # minimum subject similarity to consider (ident% * covg%, on [0,1])
    "mc=f" => \$mincov,      # minimum subject coverage percent to consider (on [0,1])
    "mi=f" => \$minidt,      # minimum subject identity percent to consider (on [0,1])
    "ml=i" => \$minaln,      # minimum alignment length to consider (in nt/aa)
    "multi" => \$multi,      # allow a query to take > 1 alignment from the same subject (default no)
    "nomap" => \$nomap,      # do not generate blast hit coverage maps for queries? (maps print by default)

    "help|?" => \$help,
    "man!" => \$man
) or pod2usage(2);

pod2usage(-exitstatus => 0, -verbose => 3) if $help || $man;

## test command line
if (-e $outdir) {
    die "The folder '$outdir' already exists!\n";
} else {
    system "mkdir -p $outdir";
    sleep 1;
    die "Unable to create output directory '$outdir'!\n" unless -e $outdir;
}
my $knownprogs = join ', ', map {"'$_'"} sort keys %blastflags;
my $knownspecs = join ', ', map {"'$_'"} sort keys %knownspecies;
die "Input file '$input' not found!\n" unless -e $input;
die "Unknown program type '$prog'!  Must be one of $knownprogs\n" unless $blastflags{$prog};
die "Unknown species '$species'!  Must be one of $knownspecs\n" unless $knownspecies{$species};

# DBs include: 1 nucleotide seq and 6 translated frames.
# DB names are: <species>_<type>.fa, where <type> is either 'aa' or 'nt'.
# Each species need a key file, names <species>_key.txt:
# key file columns: header, length, status, tag vol, read vol, then 6 frames columns:
#  each frame entry must be the string "frame,orflength,maxorf", where:
#   "frame" is one of +0 +1 +2 -0 -1 -2; "orflen" is length of max orf in frame; "maxorf" = 0 or 1 if this orf is largest of all 6 frames.  Multiple frames can be maxorf.
# the headers MUST match the headers in the nucleotide database fasta.
# protein headers also MUST match headers in protein fasta, but:
#  unlike nucleotide headers, protein headers are reconstructed by the script, as "header|frame":
#  where "header" = nucleotide header and "frame" is one of the frame values specified above.
#  thus, the script will expect to find headers like "header|+0", "header|+1", "header|+2", etc in the protein db.

chomp(my $start = `date`);
my $intemp = "$outdir/temp_blast_input.fa";
my $blastout = "$outdir/${prog}_results.txt";
my $blastfilt = "$outdir/${prog}_filtered_results.txt";
my $readme = "$outdir/filtered_results_README.txt";
my $mapdata = "$outdir/${prog}_mapdata.txt";
my $mapscript = "$outdir/${prog}_maps.R";

my (%qdata, %sdata, %fasta, $header, $hcount, $hcount2, %scored);

## parse the input file
open FA, $input;
open FA2, "> $intemp";
while (<FA>) {
    $_ =~ s/[\n\r]+$//;
    if ($_ =~ /^>(.*)/) {
        $hcount++;
        $qdata{$hcount} = [$1, 0];    # 0 will be length
        print FA2 ">$hcount\n";
    } else {
        print FA2 "$_\n" if $_;
        $fasta{$hcount} .= $_;
    }
}
close FA;
close FA2;
$qdata{$_}->[1] = length($fasta{$_}) foreach keys %fasta;

## parse the key file
open KEY, "$path/${species}_key.txt";
while (<KEY>) {
    $_ =~ s/[\n\r]+$//;
    my ($header, $length, $status, @frames) = split /\t/, $_;
    if ($prog eq 'blastn') {
        $sdata{$header} = [$length, $status];
    } else {
        foreach (@frames) {
            my ($frame, $orflen, $maxorf) = split ',', $_;
            my ($fnum) = ($frame =~ /(\d)$/);  # get number w/o +,-
            my $aalen = int( ($length-$fnum) / 3 );
	    my $header2 = $prog eq 'blastp' ? "$header|$frame" : $header;
            $sdata{$header2} = [$aalen, $status, $orflen, $maxorf];
        }
    }
}
close KEY;

my %hashtype = ("CODING" => 0, "BORDERLINE" => 1, "NONCODING" => 2);
my %YN = (1 => 'Y', 0 => 'N');
my (%Bqueries, %Bsubjects, $linecount);

## run the blast
my $flags = $blastflags{$prog}{flags};
my $db = "$path/${species}_$blastflags{$prog}{dbtype}.fa";
system "blastall -p $prog $flags -i $intemp -d $db -a $cores -m 8 -F F $flags > $blastout";

## parse the blast
open IN, $blastout;
while (<IN>) {
    (my $line = $_) =~ s/[\n\r]+$//;
    my ($query, $line2) = split /\t/, $line, 2;
    my ($qname, $qrylen) = @{ $qdata{$query} };
    my ($subject, $matchidt, $mlen, $misses, $gaps, $qpos1, $qpos2, $spos1, $spos2, $Eval, $score) = split /\t/, $line2;
    $Bqueries{ALL}{$query} = 1;
    $Bsubjects{ALL}{$subject} = 1;
    unless ($sdata{$subject}) {
        print "Contig '$subject' not found in header file: skipping alignment!\n";
        next;
    }
    next if ($minaln && $mlen < $minaln);
    my ($qstart, $qend, $sstart, $send);
    my ($subjlen, $status, $orflen, $maxorf) = @{ $sdata{$subject} };
    ($qpos1 < $qpos2) ? ( ($qstart, $qend) = ($qpos1, $qpos2) ) : ( ($qstart, $qend) = ($qpos2, $qpos1) );
    ($spos1 < $spos2) ? ( ($sstart, $send) = ($spos1, $spos2) ) : ( ($sstart, $send) = ($spos2, $spos1) );
    my $identpct = $matchidt / 100;                   # as pct subject len
    next if ($minidt && $identpct < $minidt);         # insufficient identities
    if ($prog eq 'blastp' || $prog eq 'blastx') {
        $mlen = $subjlen if ($subjlen - $mlen <= 3); # allowing 1 codon flexibility for 100% coverage
    }
    my $scovpct = $mlen / $subjlen;
    next if ($mincov && $scovpct < $mincov);      # insufficient subject coverage
    my $sindex = $identpct * $scovpct;
    next if ($minsim && $sindex < $minsim);        # insufficient overall similarity index
    $identpct = sprintf("%0.2f", 100*$identpct);	 # clean up before printing
    $scovpct = sprintf("%0.2f", 100*$scovpct);	 # clean up before printing
    my $qcovpct = sprintf("%0.2f", 100*$mlen / $qrylen);
    $sindex = sprintf("%0.2f", 100*$sindex);		 # clean up before printing
    my $idtdec;
    if ($identpct =~ /^100/) {
        $idtdec = 10;
    } elsif ($identpct =~ /^(\d+)/) {
        $idtdec = int($1/10)+1;  # identity decile (drop decimal and divide by 10)
    } else {
        print "Failed to determine identity decile from '$identpct'!\n";
    }
    my $morefields = "$qrylen\t$qcovpct\t$subjlen\t$scovpct\t$sindex\t$status\t$orflen\t$YN{$maxorf}";
    $Bqueries{PASS}{$query} = 1;
    $Bsubjects{PASS}{$subject} = 1;
    $scored{$query}{$subject}{$score}{"$qname\t$line2\t$morefields"} = join "\t", ($qname, $qrylen, $qstart, $qend, $subjlen, $sstart, $send, $idtdec, $hashtype{$sdata{$subject}->[1]});
    $linecount++;
}
close IN;

my $allq = scalar keys %{ $Bqueries{ALL} };
my $passq = scalar keys %{ $Bqueries{PASS} };
my $missq = scalar (keys %fasta) - $allq;
my $alls = scalar keys %{ $Bsubjects{ALL} };
my $passs = scalar keys %{ $Bsubjects{PASS} };
print "$linecount alignments parsed\n";
print "$passq/$allq queries passed; $missq did not blast\n";
print "$passs/$alls aligned subjects passed.\n";

my ($prebest, $best);
open OUT, "> $blastfilt";
print OUT "Query\tSubject\tIdentity%\tAlignLen\tMismatches\tGaps\tQueryStart\tQueryEnd\tSubjStart\tSubjEnd\tE-Val\tScore\tQueryLen\tQueryCov%\tSubjLen\tSubjCov%\tSubjSim%\tSubjStatus\tORFLen\tMaxOrf?\n";
open MAP, "> $mapdata";
print MAP "Query\tQ_Len\tQ_Start\tQ_End\tS_Len\tS_Start\tS_End\tIdt_Dec\tCtg_Stat\n";
foreach my $query (sort keys %scored) {
    foreach my $subject (sort keys %{ $scored{$query} }) {
        $prebest += scalar keys %{ $scored{$query}{$subject} };
        if ($multi) {
            foreach my $score (keys %{ $scored{$query}{$subject} }) {
                foreach my $line (keys %{ $scored{$query}{$subject}{$score} }) {
                    $best++;
                    print OUT "$line\n";
                    print MAP "$scored{$query}{$subject}{$score}{$line}\n";
                }
            }
        } else {
            my $topscore = [sort {$b <=> $a} keys %{ $scored{$query}{$subject} }]->[0];
            foreach my $line (keys %{ $scored{$query}{$subject}{$topscore} }) {
                $best++;
                print OUT "$line\n";
                print MAP "$scored{$query}{$subject}{$topscore}{$line}\n";
                last;     # takes first best hit only
            }
        }
    }
}
close OUT;
close MAP;

open README, "> $readme";
print README "Columns in the filtered blast output file:\n";
print README "Query\tName of the query sequence\n";
print README "Subject\tName of the subject sequence (contig)\n";
print README "Identity%\tPercent of identities in the alignment\n";
print README "AlignLen\tLength of the alignment (bp or aa)\n";
print README "Mismatches\tNumber of mismatches in the alignment\n";
print README "Gaps\tNumber of gaps in the alignment (NOT size of gaps in bp/aa)\n";
print README "QueryStart\tStart position of alignment on query sequence\n";
print README "QueryEnd\tEnd position of alignment on query sequence\n";
print README "SubjStart\tStart position of alignment on subject sequence\n";
print README "SubjEnd\tEnd position of alignment on subject sequence\n";
print README "E-Val\tExpected number of alignments of this score which could be found by chance in the subject database\n";
print README "Score\tScore of the alignment\n";
print README "QueryLen\tLength of query sequence\n";
print README "QueryCov%\tPercent of query covered by the alignment\n";
print README "SubjLen\tLength of subject sequence\n";
print README "SubjCov%\tPercent of subject covered by the alignment\n";
print README "SubjSim%\tSubject->query similarity index (subject identity % * subject coverage %)\n";
print README "SubjStatus\tAssumed coding status of subject based on ORF detection\n";
print README "ORFLen\tMaximum ORF length in this frame of the subject\n";
print README "MaxOrf?\tIs this ORF the largest in any frame of the subject? (may have > 1 largest-ORF frame)\n";
close README;
unlink $intemp;

my $OS = $^O;
my $src = $OS eq 'linux' ? '/n/projects/apa/R/apa_tools.R' : 'U:/apa/R/apa_tools.R';

## generate maps
my $literals = join ', ', map { "\"$_\"" } ('\"', '\\\|', ':', ';', '#', '@', '\\\$', '%', '&', '\\\*', '\\\{', '\\\}', '\\\[', '\\\]', ' ');
open R, "> $mapscript";
print R &mapscript;
close R;
system "nohup R --vanilla < $mapscript > ${mapscript}out";
chomp(my $end = `date`);
print "\nStart: $start\n  End: $end\nComplete!\n";
exit;


sub mapscript {

## map concept:
## black bar on x axis = query.
## layers of colored bars above x axis = subjects (similar to NCBI blast alignment chart): try to "pack" like UCSC.
## colors = alignment identity (in deciles, thus 10 colors from red=90-100% to purple=0-10%).
## hashes through bar indicate assumed status of contig:
##   none = coding, | = borderline, / = noncoding
## if $mincov < 1, light grey flanks may appear on colored alignment bars.  These are the unaligned contig flanks.

my $script = <<EOF;

source("$src")  # need .IM(), fs.legal()
mapdat <- read.delim("$mapdata", sep="\t", header=T)  # Query, Q_Len, Q_Start, Q_End, S_Len, S_Start, S_End, Idt_Dec, Ctg_Stat
idtcols <- rev(c("red","darkorange","yellow","lawngreen","green3","cyan","dodgerblue","blue","purple3","sienna4"))
qry <- unique(mapdat[,1:2])

for (i in 1:nrow(qry)) {
     name <- as.character(qry[i,1])
     .IM(name)
     if (length(grep("^gi", name)) > 0) {
          fields <- unlist(strsplit(name, "|", fixed=T))
          desc <- fields[length(fields)]
          name2 <- imgname <- paste(fields[2:length(fields)-1], collapse="|")
     } else {
          desc <- ''
          name2 <- imgname <- name
     }
     imgname <- fs.legal(imgname)
     x <- which(mapdat[,1] == name)
     x <- rev(x[order(mapdat[x,3], mapdat[x,4])])
     ht <- length(x) + 1; ht
     wd <- qry[i,2]; wd
     imght <- 300 + (ht - 10) * 4
     yht <- ifelse(ht < 10, 10, ht)
     message(i, ht, imght)
     png(paste("$outdir/", imgname, "_map.png", sep=""), height=imght, width=1000)
     par(cex=1.2)
     plot(0, 0, col=0, xlim=c(1,wd), ylim=c(-2,yht), xlab="", ylab="", yaxt="n", main=paste(name2,": ",length(x)," hits\n",desc,sep=""))
     rect(0, -2, wd, 0, col=8, bord=1)      # gene base: grey with black border
     for (j in 1:length(x)) {
          color <- idtcols[mapdat[x[j],8]]
          f5s <- mapdat[x[j],3] - mapdat[x[j],6]
          f5e <- mapdat[x[j],3]
          f3s <- mapdat[x[j],4]
          f3e <- mapdat[x[j],4] + (mapdat[x[j],5] - mapdat[x[j],7])
          rect(f5s, j-0.3, f5e, j+0.3, col=8, bord=8)            # 5' flank
          rect(f3s, j-0.3, f3e, j+0.3, col=8, bord=8)            # 3' flank
          rect(f5e, j-0.3, f3s, j+0.3, col=color, bord=color)    # alignment
          rect(f5e, -2, f3s, 0, col=1, bord=1)    # fill in aligned area of gene 
     }
     dev.off()
     
     png("$outdir/blast_map_color_key.png", height=600, width=600)
     par(cex=1.3)
     plot(1, 1, col=0, xlim=c(1,8), ylim=c(0.5,10), axes=F, xlab="", ylab="", main="Blast Hit Identity Colors")
     for (i in 1:10) { 
          rect(1, i-0.55, 5, i+0.15, col=idtcols[i], bord=idtcols[i])
          text(6, i, pos=4, labels=paste("> ", 10*(i-1), "% Identity", sep=""))
     }
     dev.off()
}

EOF
return $script;
}

