#!/usr/bin/env perl
use Excel::Writer::XLSX;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use strict;
no strict 'refs';

# alignToRNA_summary.pl [-t targets.txt] [-f fasta] [-o output_folder] --no-trim --verbose


## Inputs
my $targets;    # a 2-col file with col 1 = fastq path, col 2 = alias to display instead of filename
my $outdir;     # output directory, do not clobber
my $notrim;     # deactivate adapter trimming
my $verbose;    # echo system commands?


GetOptions("t=s" => \$targets, "f=s" => \$fasta, "o=s" => \$outdir, "verbose" => \$verbose, "no-trim" => \$notrim) || pod2usage(2);


## Globals
my @fastqs;    # sample fastq paths
my @samples;   # sample names
my @usamples;  # unique sample names
my @bands;     # sample bands
my $useBands;  # are bands being used in the analysis? (1/0)
my $N;         # N rows in $targets (0-based)
my $M;         # N unique samples (0-based)
my @ncRNAs;    # ncRNA names for this index
my %idxData;   # idxstats data
my %sumData;   # summarized idxstats data
my %rnaData;   # ncRNA data
my %errData;   # bowtie.err data, etc
#my %trmData;   # trimmomatic summary   ## NOT CURRENTLY IN USE
my %allnames;  # sample/replicate alias tracking
my @classes;   # RNA classes
my @stats;     # bowtie.err fields
my $fatal;     # fatal error encountered


## Parse targets file
## Determine if 'band' column is present -- substitute '1' if not, just as a placeholder
## Make targets-ordered arrays of fastq files, sample IDs, and bands
## Make unique-sample-IDs array
my %already;
$targets = abs_path($targets);
open my $IN1, '<', $targets or die "$0: Cannot open targets file '$targets' for reading: $!\n";
while (<$IN1>) {
    next if $. == 1;  ##### EXPECTS HEADER
    $_ =~ s/[\n\r]+$//;
    my @data = split /\t/, $_;
    die "$0: Targets file must be 2 or 3 columns: fastq path, sample name, optional band name!\n" if $#data < 1;
    push @fastqs, $data[0];
    push @samples, $data[1];
    my $band = $data[2] || 1;
    push @bands, $band;
    push @usamples, $data[1] unless $already{$data[1]};
    $already{$data[1]} = 1;
}
$N = $.-2;  # -header & 0-based
$M = $#usamples;
close $IN1;


## Did the user specify bands or not?  If not, no band-wise output will be generated, and fastq aliases will lack the band field.
my $ubands = scalar keys %{ {map {($_=>1)} @bands} };
$useBands = $ubands == 1 && $bands[0] == 1 ? 0 : 1;


## Read ncRNA fasta to which reads will be aligned
## Group ncRNAs by class, e.g. rRNA, snoRNA, miRNA, etc: fasta headers are always "class|gene"
## Initialize reporting objects with zeroes to ensure printable values
## Calculate ncRNA lengths
open my $IN3, '<', $fasta or die "$0: Failed to open ncRNA fasta '$fasta': $!\n";
my ($ncRNA, $seq);
while (<$IN3>) {
    $_ =~ s/[\n\r]+$//;
    if ($_ =~ /^>(.*)/) {
	if ($seq) {
	    $rnaData{$ncRNA}->[2] = length($seq);
	    $seq = '';
	}
	$ncRNA = $1;
	push @ncRNAs, $ncRNA;
	my ($class, $gene) = (split /\|/, $ncRNA, 2);
	$rnaData{$ncRNA} = [$class, $gene, 0];
	foreach my $i (0..$N) {
	    my $samp = $samples[$i];
	    my $band = $bands[$i];
	    $idxData{$ncRNA}{$samp}{$band} = [0, 0];  # ensure printable values
	    $sumData{$class}{$samp}{$band} = [0, 0];
	    $idxData{$ncRNA}{$samp}{ALL} = [0, 0];
	    $sumData{$class}{$samp}{ALL} = [0, 0];
	}
    } else {
	$seq .= $_;
    }
}
$rnaData{$ncRNA}->[2] = length($seq);
close $IN3;


## Unique classes; stat field names from bowtie-error files
@classes = sort keys %sumData;
@stats = qw/ READS UNALIGNED ALL_ALIGNED MONO_ALIGNED MULTI_ALIGNED /;


## For each fastq file in the targets.txt, run the following block
foreach my $i (0..$N) {
    
    my $path = $fastqs[$i];
    my $samp = $samples[$i];
    my $band = $bands[$i];
    
    ## alias: include band data or not?
    my ($alias, $alias, $blurb);
    if ($useBands) {
	$alias = "$samp.$band";
	$blurb = "'$samp $band'";
    } else {
	$alias = $samp;
	$blurb = "'$samp'";
    }
    my $aliasdir = "$outdir/$alias";
    my $prefix = "$aliasdir/$alias";
    my $genodir = "$prefix.genome";
    
    ## trimmomatic stats files
    my $allhisto = "$prefix.trim.histo.txt";
    my $alnhisto = "$prefix.align.histo.txt";
    my $unahisto = "$prefix.unaln.histo.txt";
    my $trimlog = "$prefix.trimmomatic.log";
    my $trimtxt = "$prefix.trimmomatic.txt";
    
    ## alignment output files etc.
    my $logfile = "$prefix.log";
    my $alignlog = "$prefix.align_summary.txt";
    my $galignlog = "$genodir/align_summary.txt";
    my $idxfile = "$prefix.idxstats.txt";
    my $gidxfile = "$genodir/accepted_hits.idxstats.txt";
    my $pidxfile = "$genodir/accepted_hits_primary.idxstats.txt";
    my $genobam = "$genodir/accepted_hits.bam";
    my $primbam = "$genodir/accepted_hits_primary.bam";
    my $mapfq = "$prefix.mapped.fq.gz";
    my $unmapfq = "$prefix.unmapped.fq.gz";
    my $primbw = "$genodir/accepted_hits_primary.APM.bw";
    
    unless ($notrim) {
	my @trimstats = split / /, `grep ^Input $trimtxt`;  # e.g. "Input Reads: 5636013 Surviving: 5061126 (89.80%) Dropped: 574887 (10.20%)"
	$trimstats[$_] =~ s/[\(\)%]//g foreach (5,8);
#	$trmData{INPUT}{$samp}{$band} = [$trimstats[2], 1];
#	$trmData{KEEP}{$samp}{$band} = [$trimstats[4], $trimstats[5]/100];
#	$trmData{DROP}{$samp}{$band} = [$trimstats[7], $trimstats[8]/100];
	
## knocking this block out for now (deprecated)
##	if (0 && !$notrim) {
##	    &execute("gunzip -c $mapfq | sed -n '2~4p' | awk '{print length(\$1)}' | /home/apa/local/bin/colTally - | sort -n -k2,2 > $alnhisto");
##	    &execute("gunzip -c $unmapfq | sed -n '2~4p' | awk '{print length(\$1)}' | /home/apa/local/bin/colTally - | sort -n -k2,2 > $unahisto");
##	    my %tmp;
##	    $tmp{TOTAL} = histoparse($allhisto);
##	    $tmp{ALIGN} = histoparse($alnhisto);
##	    $tmp{UNALN} = histoparse($unahisto);
##	    my ($Lmin, $Lmax) = (sort {$a <=> $b} keys %{ $tmp{TOTAL} })[0,-1];
##	    open my $OUT1, '>', $allhisto or die "$0: Cannot open histogram file '$allhisto' for writing: $!\n";   # OVERWRITE
##	    print $OUT1 "ReadLen\tTRIMMED\tALIGNED\tUNALIGNED\n";
##	    print $OUT1 "$_\t". ($tmp{TOTAL}{$_}||0) ."\t". ($tmp{ALIGN}{$_}||0) ."\t". ($tmp{UNALN}{$_}||0) ."\n" foreach ($Lmin..$Lmax);
##	    close $OUT1;
##	    &execute("rm -f $alnhisto $unahisto");
##	}
	if (0) {   ## %trmData is not used by anything at this point
	    open my $IN4, '<', $allhisto or die "$0: Cannot open aligned read-length histogram file '$allhisto': $!\n";
	    while (<$IN4>) {
		next if $. == 1;
		$_ =~ s/[\n\r]+$//;
		my ($readlen, $class, $multi) = split /\t/, $_;
		$trmData{BYLEN}{$readlen}{$samp}{$band}{TOTAL}{ALL}++;
		if ($class eq '*') {  # not aligned
		    $trmData{BYLEN}{$readlen}{$samp}{$band}{UNALIGN}{ALL}++;
		} else {
		    $trmData{BYLEN}{$readlen}{$samp}{$band}{ALIGN}{ALL}++;
		    $trmData{BYLEN}{$readlen}{$samp}{$band}{ALIGN}{$class}++;
		    if ($multi) {
			$trmData{BYLEN}{$readlen}{$samp}{$band}{MULTI}{ALL}++;
			$trmData{BYLEN}{$readlen}{$samp}{$band}{MULTI}{$class}++;
		    } else {
			$trmData{BYLEN}{$readlen}{$samp}{$band}{MONO}{ALL}++;
			$trmData{BYLEN}{$readlen}{$samp}{$band}{MONO}{$class}++;
		    }
		}
	    }
	    close $IN4;
	}
    }
    
    ## parse bowtie-err file for overall alignment stats
    ## record as [N reads, % Total]
    open my $IN5, '<', $alignlog or die "$0: Cannot open bowtie error file '$alignlog': $!\n";
    while (<$IN5>) {
	$_ =~ s/%//;
	if ($_ =~ /^(\d+) reads/) {
	    $errData{READS}{$samp}{$band} = [$1, 1];    ### KEY NAMES MUST MATCH @stats ENTRIES
	} elsif ($_ =~ /^\s+(\d+) \(([\d.]+)\) were unpaired/) {
##	    $errData{UNPAIRED}{$samp}{$band} = [$1, $2/100];  # all treated as unpaired...
	} elsif ($_ =~ /^\s+(\d+) \(([\d.]+)\) aligned 0/) {
	    $errData{UNALIGNED}{$samp}{$band} = [$1, $2/100];
	} elsif ($_ =~ /^\s+(\d+) \(([\d.]+)\) aligned exactly/) {
	    $errData{MONO_ALIGNED}{$samp}{$band} = [$1, $2/100];
	} elsif ($_ =~ /^\s+(\d+) \(([\d.]+)\) aligned >1/) {
	    $errData{MULTI_ALIGNED}{$samp}{$band} = [$1, $2/100];
	} elsif ($_ =~ /^([\d.]+) overall/) {
	    $errData{ALL_ALIGNED}{$samp}{$band} = [$errData{MONO_ALIGNED}{$samp}{$band}->[0]+$errData{MULTI_ALIGNED}{$samp}{$band}->[0], $1/100];
	}
    }
    close $IN5;
    
    ## ncRNAs quantitated by samtools idxstats; parse this result
    ## record as [N aligns, RPKM]
    ## sumData just racking up counts
    open my $IN6, '<', $idxfile or die "$0: Cannot open idxstats output '$idxfile' for reading: $!\n";
    while (<$IN6>) {
	$_ =~ s/[\n\r]+$//;
	my ($ncRNA, $len, $aligns, $unaligns) = (split /\t/, $_);
	unless ($ncRNA eq '*') {
	    $idxData{$ncRNA}{$samp}{$band} = [$aligns, (1E9*$aligns)/($rnaData{$ncRNA}->[2]*$errData{ALL_ALIGNED}{$samp}{$band}->[0])];
	    $sumData{ $rnaData{$ncRNA}->[0] }{$samp}{$band}->[0] += $aligns;
	}
    }
    close $IN6;
    
    ## add RPKMs to sumData, after all counting finished
    $sumData{$_}{$samp}{$band}->[1] = $sumData{$_}{$samp}{$band}->[0] / $errData{ALL_ALIGNED}{$samp}{$band}->[0] foreach keys %sumData;
    
    ## whole-sample counts; above were for sample+band
    $idxData{$_}{$samp}{ALL}->[0] += $idxData{$_}{$samp}{$band}->[0] foreach keys %idxData;
    $errData{$_}{$samp}{ALL}->[0] += $errData{$_}{$samp}{$band}->[0] foreach keys %errData;
    $sumData{$_}{$samp}{ALL}->[0] += $sumData{$_}{$samp}{$band}->[0] foreach keys %sumData;
    
}


## Whole-sample percents and RPKMs; above were for sample+band
foreach my $samp (@usamples) {
    $idxData{$_}{$samp}{ALL}->[1] = (1E9*$idxData{$_}{$samp}{ALL}->[0])/($rnaData{$_}->[2]*$errData{ALL_ALIGNED}{$samp}{ALL}->[0]) foreach keys %idxData;
    $errData{$_}{$samp}{ALL}->[1] = $errData{$_}{$samp}{ALL}->[0] / $errData{READS}{$samp}{ALL}->[0] foreach keys %errData;
    $sumData{$_}{$samp}{ALL}->[1] = $sumData{$_}{$samp}{ALL}->[0] / $errData{ALL_ALIGNED}{$samp}{ALL}->[0] foreach keys %sumData;
}



########## OUTPUT



## Excel layout prep: file names, header sizes, column widths, levels to summarize
my %levelData = ('band' => ["$results/ByBand.xlsx", 2, $N], 'sample' => ["$results/BySample.xlsx", 1, $M], 'trim' => ["$results/TrimStats.xlsx", 1, $M]);  # [output file, N header lines, N data sets]
my ($sampleWidth, $spacer, $summColWidth, $classWidth, $ncrnaWidth, $lenWidth) = (14, 6, 20, 10, 70, 10);
my (%colwidths, %headers);
my @useLevels = $useBands ? qw/ band sample / : qw/ sample /;


## Produce one file per level (per sample, and per band if using bands)
foreach my $level (@useLevels) {
    
    ## Set up 
    my ($h, $n) = @{ $levelData{$level} }[1,2];  # N header lines, N data sets
    my @cg = ('Class','Gene','Length');
    my @blank = ('','','');
    my @cg1 = $level eq 'sample' ? @cg : @blank;
    my @levelSamples = $level eq 'sample' ? @usamples : @samples;
    
    ## Different header structures if using global+local or local-only
    $colwidths{sum}  = [$summColWidth, (map {$sampleWidth} (0..$n))];
    $headers{sum}{1} = ['', @levelSamples];
    $headers{sum}{2} = ['', @bands];  # only used for 'band' level
    $colwidths{rna}  = [$classWidth, $ncrnaWidth, $lenWidth, (map {$sampleWidth} (0..$n))];
    $headers{rna}{1} = [@cg1, @levelSamples];
    $headers{rna}{2} = [@cg, @bands];  # only used for 'band' level
    
    ## Initialize workbook
    ## Create formats
    unlink $levelData{$level}->[0];   # ABSOLUTELY ENSURE WORKBOOK IS CLOSED, at least from the filesystem's perspective
    my $workbook = Excel::Writer::XLSX->new($levelData{$level}->[0]);
    my $header_format = $workbook->add_format();
    my $count_format = $workbook->add_format();
    my $percent_format = $workbook->add_format();
    my $float_format = $workbook->add_format();
    $header_format->set_bold(1);
    $header_format->set_align('center');
    $count_format->set_num_format(3);
    $percent_format->set_num_format(10);
    $float_format->set_num_format('0.0000');
    
    ## Initialize 'Summary' sheet
    my $worksheet = $workbook->add_worksheet("Summary");
    $worksheet->set_column($_, $_, $colwidths{sum}->[$_]) foreach (0..$#{ $colwidths{sum} });
    $worksheet->freeze_panes($h, 0);
    my $row = -1;
    
    ## Write 'Summary' header lines
    foreach my $i (1..$h) {
	$worksheet->write_row(++$row, 0, $headers{sum}{$i}, $header_format);
    }
    
    ## Write 'Summary' data set 1: bowtie err stats, as counts
    foreach my $stat (@stats) {
	my @dat = @{ &excelify($errData{$stat}, 0, $level) };
	$worksheet->write_row(++$row, 0, [$stat, @dat], $count_format);
    }
    
    ## Write 'Summary' data set 2: ncRNA class-wise stats, as counts
    $row++;  # spacer
    foreach my $class (@classes) {
	my @dat = @{ &excelify($sumData{$class}, 0, $level) };
	$worksheet->write_row(++$row, 0, [$class, @dat], $count_format);
    }
    
    ## Write 'Summary' data set 3: bowtie err stats, as percents
    $row += 2;  # double spacer
    foreach my $stat (@stats) {
	my @dat = @{ &excelify($errData{$stat}, 1, $level) };
	$worksheet->write_row(++$row, 0, [$stat, @dat], $percent_format);
    }
    
    ## Write 'Summary' data set 4: ncRNA class-wise stats, as percents
    $row++;  # spacer
    foreach my $class (@classes) {
	my @dat = @{ &excelify($sumData{$class}, 1, $level) };
	$worksheet->write_row(++$row, 0, [$class, @dat], $percent_format);
    }
    $worksheet->write(++$row, 0, '(as % ALL_ALIGNED)');
    
    
    ## WRITE NCRNA-COUNTS, NCRNA-FPKMS SHEETS
    
    ## ncRNA sheet names, formats
    my @rnaSheets = (['ncRNA_Counts', 0, $count_format], ['ncRNA_RPKMs', 1, $float_format]);
    
    ## Two sheets: '0' = counts, '1' = FPKMs
    foreach my $i (0,1) {
	
	## Initialize sheet
	my $worksheet = $workbook->add_worksheet($rnaSheets[$i][0]);
	$worksheet->set_column($_, $_, $colwidths{rna}->[$_]) foreach (0..$#{ $colwidths{rna} });
	$worksheet->freeze_panes($h, 0);
	my $row = -1;
	
	## Write header
	foreach my $i (1..$h) {
	    $worksheet->write_row(++$row, 0, $headers{rna}{$i}, $header_format);
	}
	
	## Write counts/FPKMs data per ncRNA
	foreach my $ncRNA (@ncRNAs) {
	    my @dat = @{ &excelify($idxData{$ncRNA}, $rnaSheets[$i][1], $level) };
	    $worksheet->write_row(++$row, 0, [@{ $rnaData{$ncRNA} }, @dat], $rnaSheets[$i][2]);
	}
    }
    
    ## CLOSE ENTIRE WORKBOOK
    $workbook->close();
}


unless ($notrim) {
#    		$trmData{INPUT}{$samp}{$band} = [$trimstats[2], 1];
#		$trmData{KEEP}{$samp}{$band} = [$trimstats[4], $trimstats[5]/100];
#		$trmData{DROP}{$samp}{$band} = [$trimstats[7], $trimstats[8]/100];
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{TOTAL}{ALL}++;
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{UNALIGN}{ALL}++;
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{ALIGN}{ALL}++;
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{ALIGN}{$class}++;
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{MULTI}{ALL}++;
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{MULTI}{$class}++;
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{MONO}{ALL}++;
#		$trmData{BYLEN}{$readlen}{$samp}{$band}{MONO}{$class}++;

    if (0) {
   
    ## Set up
    my $level = 'trim';
    my ($h, $n) = @{ $levelData{$level} }[1,2];  # N header lines, N data sets
    my @cg = ('Class','Gene','Length');
    my @blank = ('','','');
    my @cg1 = $level eq 'sample' ? @cg : @blank;
    my @levelSamples = $level eq 'sample' ? @usamples : @samples;
    
    ## Different header structures if using global+local or local-only
    # use $colwidths{sum}, $headers{sum} as above
    $colwidths{rna}  = [$classWidth, $ncrnaWidth, $lenWidth, (map {$sampleWidth} (0..$n))];
    $headers{rna}{1} = [@cg1, @levelSamples];
    $headers{rna}{2} = [@cg, @bands];  # only used for 'band' level
    
    ## Initialize workbook
    ## Create formats
    unlink $levelData{$level}->[0];   # ABSOLUTELY ENSURE WORKBOOK IS CLOSED, at least from the filesystem's perspective
    my $workbook = Excel::Writer::XLSX->new($levelData{$level}->[0]);
    my $header_format = $workbook->add_format();
    my $count_format = $workbook->add_format();
    my $percent_format = $workbook->add_format();
    my $float_format = $workbook->add_format();
    $header_format->set_bold(1);
    $header_format->set_align('center');
    $count_format->set_num_format(3);
    $percent_format->set_num_format(10);
    $float_format->set_num_format('0.0000');
    
    ## Initialize 'Summary' sheet
    my $worksheet = $workbook->add_worksheet("Summary");
    $worksheet->set_column($_, $_, $colwidths{sum}->[$_]) foreach (0..$#{ $colwidths{sum} });
    $worksheet->freeze_panes($h, 0);
    my $row = -1;
    
    ## Write 'Summary' header lines
    foreach my $i (1..$h) {
	$worksheet->write_row(++$row, 0, $headers{sum}{$i}, $header_format);
    }
    
    ## Write 'Summary' data set 1: bowtie err stats, as counts
    foreach my $stat (@stats) {
	my @dat = @{ &excelify($errData{$stat}, 0, $level) };
	$worksheet->write_row(++$row, 0, [$stat, @dat], $count_format);
    }
    
    ## Write 'Summary' data set 2: ncRNA class-wise stats, as counts
    $row++;  # spacer
    foreach my $class (@classes) {
	my @dat = @{ &excelify($sumData{$class}, 0, $level) };
	$worksheet->write_row(++$row, 0, [$class, @dat], $count_format);
    }
    
    ## Write 'Summary' data set 3: bowtie err stats, as percents
    $row += 2;  # double spacer
    foreach my $stat (@stats) {
	my @dat = @{ &excelify($errData{$stat}, 1, $level) };
	$worksheet->write_row(++$row, 0, [$stat, @dat], $percent_format);
    }
    
    ## Write 'Summary' data set 4: ncRNA class-wise stats, as percents
    $row++;  # spacer
    foreach my $class (@classes) {
	my @dat = @{ &excelify($sumData{$class}, 1, $level) };
	$worksheet->write_row(++$row, 0, [$class, @dat], $percent_format);
    }
    $worksheet->write(++$row, 0, '(as % ALL_ALIGNED)');
    
    
    ## WRITE NCRNA-COUNTS, NCRNA-FPKMS SHEETS
    
    ## ncRNA sheet names, formats
    my @rnaSheets = (['ncRNA_Counts', 0, $count_format], ['ncRNA_RPKMs', 1, $float_format]);
    
    ## Two sheets: '0' = counts, '1' = FPKMs
    foreach my $i (0,1) {
	
	## Initialize sheet
	my $worksheet = $workbook->add_worksheet($rnaSheets[$i][0]);
	$worksheet->set_column($_, $_, $colwidths{rna}->[$_]) foreach (0..$#{ $colwidths{rna} });
	$worksheet->freeze_panes($h, 0);
	my $row = -1;
	
	## Write header
	foreach my $i (1..$h) {
	    $worksheet->write_row(++$row, 0, $headers{rna}{$i}, $header_format);
	}
	
	## Write counts/FPKMs data per ncRNA
	foreach my $ncRNA (@ncRNAs) {
	    my @dat = @{ &excelify($idxData{$ncRNA}, $rnaSheets[$i][1], $level) };
	    $worksheet->write_row(++$row, 0, [@{ $rnaData{$ncRNA} }, @dat], $rnaSheets[$i][2]);
	}
    }
    
    ## CLOSE ENTIRE WORKBOOK
    $workbook->close();
    
    }  # END IF 0
    
} # END UNLESS $NOTRIM

exit;





#################################################################################################################
#################################################################################################################
#################################################               #################################################
#################################################  SUBROUTINES  #################################################
#################################################               #################################################
#################################################################################################################
#################################################################################################################





sub excelify {
    
    # turn specific data slices into hashes of Excel-printable arrays
    
    my ($HASH, $J, $LEVEL) = @_;
    my @DAT;
    if ($LEVEL eq 'band') {
	push @DAT, $$HASH{ $samples[$_] }{ $bands[$_] }->[$J] foreach (0..$N);
    } elsif ($LEVEL eq 'sample') {
	push @DAT, $$HASH{$_}{ALL}->[$J] foreach @usamples;
    } else {
	die "$0: unknown level type '$LEVEL'!\n";
    }
    return \@DAT;
}


sub histoparse {
    
    # hashify a trimmomatic read-length histogram file
    
    my $HIST = shift;
    my %HASH;
    
    open my $IN7, '<', $HIST or die "$0: Cannot open histogram file '$HIST': $!\n";
    while (<$IN7>) {
	$_ =~ s/[\n\r]+$//;
	my ($N, $L) = split /\t/, $_;
	$HASH{$L} = $N;
    }
    close $IN7;
    
    return \%HASH;
}
