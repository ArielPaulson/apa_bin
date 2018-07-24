#!/usr/bin/env perl
use strict;
use Roman;
use Data::Dumper;


# #!/n/site/inst/Linux-x86_64/sys/bin/perl

### SUMMARY OF FUNCTIONS:
## timestamp    : three different types of platform-independent timestamps
## monthconv    : convert month names to numbers and vice versa
## benchmark    : platform-independent system resource stamps
## logreport    : send a message to screen & logfile simultaneously
## revcomp      : reverse complement for DNA/RNA, with full degeneracy/masking support
## word_exp     : get the expected frequency of a DNA/RNA/AA word, given background character frequencies
## numjust      : pad a digit with leading zeroes
## blockify     : convert a sequence to a fasta block
## unique       : return unique values from an array, in order of appearance
## chrsort      : sort a list of chromosomes the way I prefer.  Special handling for roman chrnums (yeast) and drosophila chrs
## get_memedata : get a standard hash of data from a MEME run
## UTR_weld     : take a hash of CDS entries + hash of UTR entries and return a hash of exons

### UNDER CONSTRUCTION: may have some functionality at this time, or may not
## runtime      : takes 2 'FULL' timestamps and returns elapsed time, in various scales
## getfile      : read a tabular text file into an array-of-arrays or hash-of-arrays
## writefile    : write an array, array-of-arrays, or hash-of-arrays back to a tabular text file
## getfasta     : read a fasta file into a hash
## writefasta   : write a fasta file from a hash
## get_fimodata : get a standard hash of data from a FIMO run
## get_tomdata  : get a standard hash of data from a Tomtom run



sub timestamp {
    
    ## returns date, time, or date+time timestamps
    
    my $timetype = shift;	# FULL = date + time | DATE = date | TIME = time
    my $timestamp;
    
    my %daynames = (0,'Sun', 1,'Mon', 2,'Tue', 3,'Wed', 4,'Thu', 5,'Fri', 6,'Sat');
    my %monthnames = (0,'Jan', 1,'Feb', 2,'Mar', 3,'Apr', 4,'May', 5,'Jun', 6,'Jul', 7,'Aug', 8,'Sep', 9,'Oct', 10,'Nov', 11,'Dec');
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    
    if ($timetype eq 'FULL') {
	$timestamp = "$daynames{$wday} $monthnames{$mon} ".sprintf("%02d %4d %02d:%02d:%02d", $mday,$year+1900,$hour,$min,$sec);
    } elsif ($timetype eq 'DATE') {
	$timestamp = "$daynames{$wday} $monthnames{$mon} ".sprintf("%02d %4d", $mday,$year+1900);
    } elsif ($timetype eq 'TIME') {
	$timestamp = sprintf("%02d:%02d:%02d", $hour,$min,$sec);
    }
    
    return $timestamp;
}



sub runtime {   # UNDER CONSTRUCTION
    
    ## takes 2 'FULL' timestamps and returns the difference, in the specified timescale
    
    ### requires 2 timestamp('FULL') entries
    my (%timestamps, %timeparts, %timediff, $scale);
    ($timestamps{1}, $timestamps{2}, $scale) = @_;  # timestamps with format e.g. 'Thu Jun 23 2011 12:57:08'; for $scale see %scales below
    my %scales = map {($_=>1)} qw/ FULL YEAR MONTH DAY HOUR MINUTE SECOND /;  # time scale to return results in; 'FULL' = full breakout
    # note that scale 'MONTH' will return month breakdown as months elapsed + fraction of end-month elapsed
    my %monthdays = (1,31, 2,28, 3,31, 4,30, 5,31, 6,30, 7,31, 8,31, 9,30, 10,31, 11,30, 12,31);
    my %monthnums = (1,'Jan', 2,'Feb', 3,'Mar', 4,'Apr', 5,'May', 6,'Jun', 7,'Jul', 8,'Aug', 9,'Sep', 10,'Oct', 11,'Nov', 12,'Dec');
    my %timescales = ('YEAR',31536000, 'DAY',86400, 'HOUR',3600, 'MIN',60);  # times in seconds
    foreach my $t (1,2) {
	my ($dname, $mon, $date, $yr, $time) = split /\s+/, $timestamps{$t};
	my ($hr, $min, $sec) = split /:/, $time;
	$timeparts{$t}{YEAR} = $yr;
	$timeparts{$t}{MONTH} = monthconv($mon,2);
	$timeparts{$t}{DAY} = $date;
	$timeparts{$t}{HOUR} = $hr;
	$timeparts{$t}{MIN} = $min;
	$timeparts{$t}{SEC} = $sec;
    }
    $timediff{$_} = $timeparts{2}{$_} - $timeparts{1}{$_} foreach qw/ YEAR MON DAY HOUR MIN SEC /;
    $timediff{DAY}-- if ($timediff{DAY} > 0 && $timediff{HOUR} != 0);     # whole days elapsed
    $timediff{ALLDAYS} = $timediff{DAYS};   # initialize with days elapsed during end month




    if ($timediff{MONTH} > 0) {     # a month or more
	$timediff{ALLDAYS} += $monthdays{$_} foreach ($timeparts{1}{MONTH}..$timeparts{2}{MONTH}-1);  # add all days for all months
	if ($timediff{DAY} != 0) {  # remainder days exist -- adjust day count
	    $timediff{ALLDAYS} += $monthdays{ $timeparts{1}{MONTH} } - $timeparts{1}{DAY};  # add days elapsed during start month
	    $timediff{MONTH}--;     # whole months elapsed
	}
    }
    if ($timediff{YEAR} > 0 && $timediff{MONTH} != 0) {
	$timediff{YEAR}--;  # whole years elapsed
    }
    
    my $seconds = $timediff{YEAR}*31536000 + $timediff{ALLDAYS}*86400 + $timediff{HOUR}*3600 + $timediff{MIN}*60 + $timediff{SEC};
    my @elapsed;
    push @elapsed, "Elapsed | $seconds seconds |";
    push @elapsed, "$timediff{YEAR} years" if $timediff{YEAR};
    push @elapsed, "$timediff{MONTH} months" if $timediff{MONTH};
    push @elapsed, "$timediff{DAY} days" if $timediff{DAY};
    push @elapsed, "$timediff{HOUR} hours" if $timediff{HOUR};
    push @elapsed, "$timediff{MIN} minutes" if $timediff{MIN};
    push @elapsed, "$timediff{SEC} seconds" if $timediff{SEC};
    my $elapsed = join ', ', @elapsed;
    return "$elapsed\n";
}



sub monthconv {
    
    ## converting month names to numbers or vice versa
    
    my ($inmonth, $dir) = @_;   # $dir = 1 for num -> name | $dir = 2 for name -> num
    my %monthnames1 = (1,'Jan', 2,'Feb', 3,'Mar', 4,'Apr', 5,'May', 6,'Jun', 7,'Jul', 8,'Aug', 9,'Sep', 10,'Oct', 11,'Nov', 12,'Dec');
    my %monthnames2 = reverse %monthnames1;
    
    if ($dir == 1) {
	$inmonth =~ s/^0+//;  # drop any leading zeroes
	return $monthnames1{$inmonth};
    } elsif ($dir == 2) {
	my $outmonth = $monthnames2{$inmonth};
	$outmonth = "0$outmonth" if $outmonth < 10;
	return $outmonth;
    } else {
	die "Conversion direction '$dir' must be 1 (num->name) or 2 (name->num)!\n";
    }
}



sub benchmark {
    
    ## captures 'top' info on a process (Linux -OR- Windows) and returns it, annotated (or writes it to file)
    
    my ($pid, $platform, $file, $msg, $script) = @_;	# $msg, $script, $file can be null
#    print "Received benchmark request: '$pid', '$msg', '$script', '$platform', '$file'\n";
    my $string;
    if ($platform eq 'LINUX') {
	my @output = split /\n/, `top -b -p $pid -n 1`;
	my @stuff = split /\s+/, $output[7];
	$string = join "\t", ($msg, $script, @stuff);
    } elsif ($platform eq 'WIN2K') {
	my @output = split /\n/, `tasklist /V /FI "PID eq $pid"`;
	my @stuff = split /\s+/, $output[3];
	$string = join "\t", ($msg, $script, @stuff);
    } else {
	$string = "\nBenchmark subroutine doesn't know what to do on $platform operating systems!\nNo RAM usage data will be available.";
	warn $string;
    }
    if ($file) {
	open OUT, ">> $file" or warn "\n&benchmark cannot append to $file: $!\n";
	print OUT "$string\n";
	close OUT;
    } else {
	return $string;
    }
}



sub logreport {
    
    ## prints a message to STDOUT and to a given logfile 
    
    my ($msg, $logfile, $quiet) = @_;
    $msg =~ s/[\n\r]+$//;	# so you don't have to remember if it needs newlines
    if ($logfile) {
	open OUT, ">> $logfile" or warn "\nlogreport cannot append to $logfile: $!\n";
	print OUT "$msg\n";
	close OUT;
    }
    print "$msg\n" unless $quiet;
}



sub word_exp {
    
    ## expected frequency of a DNA/RNA/AA word given known background frequencies
    ## divide search space size by return value to get expected number of words
    
    my ($word, $alpha, $ref) = @_;   # $word = motif, $ref = bkg hash ref (%s or counts), $alpha = DNA, RNA, AA, or FREE (last = use bkgfreq as-is; don't test)
    
    my %alphabets = (
	'DNA' => { 'required' => { map {($_=>1)} qw/ A C G T / }, 
		   'allowed' => { map {($_=>1)} qw/ A C G T R Y S W K M H D V B N / } 
	},
	'RNA' => { 'required' => { map {($_=>1)} qw/ A C G U / }, 
		   'allowed' => { map {($_=>1)} qw/ A C G U R Y S W K M H D V B N / } 
	},
	'AA'  => { 'required' => { map {($_=>1)} qw/ A C D E F G H I K L M N P Q R S T V W Y / }, 
		   'allowed' => { map {($_=>1)} qw/ A C D E F G H I K L M N P Q R S T V W Y B Z J X U O / } 
	},
	'FREE' => 1
    );
    my %convdegen = (
	'DNA' => {   'R' => [qw/ A G /],
		     'Y' => [qw/ C T /],
		     'S' => [qw/ C G /],
		     'W' => [qw/ A T /],
		     'K' => [qw/ G T /],
		     'M' => [qw/ A C /],
		     'B' => [qw/ C G T /],
		     'D' => [qw/ A G T /],
		     'H' => [qw/ A C T /],
		     'V' => [qw/ A C G /],
		     'N' => [qw/ A C G T /]
	},
	'RNA' => {   'R' => [qw/ A G /],
		     'Y' => [qw/ C U /],
		     'S' => [qw/ C G /],
		     'W' => [qw/ A U /],
		     'K' => [qw/ G U /],
		     'M' => [qw/ A C /],
		     'B' => [qw/ C G U /],
		     'D' => [qw/ A G U /],
		     'H' => [qw/ A C U /],
		     'V' => [qw/ A C G /],
		     'N' => [qw/ A C G U /]
	},
	'AA' => {    'B' => [qw/ N D /],
		     'Z' => [qw/ E Q /],
		     'J' => [qw/ I L /],
		     'X' => [qw/ A C D E F G H I K L M N P Q R S T V W Y /]
	}
    );
    
    $word = "\u$word";
    $alpha = "\U$alpha";
    die "Unknown alphabet '$alpha'!  Must be one of 'DNA', 'RNA', 'AA', or 'FREE'.\n" unless $alphabets{$alpha};

    ### SET UP BKG FREQ HASH
    my ($sum, %freqs);   # background letter frequencies
    if ($ref) {  # known bkg; but are they frequencies?
	my ($sub1, $sum1);
	foreach (keys %$ref) {
	    $sub1++ if $$ref{$_} < 1;
	    $sum1 += $$ref{$_};
	    die "Cannot have negative background frequencies!  &word_exp halting.\n" if $_ < 0;
	}
	if ($sub1 == scalar (keys %$ref)) {   # all frequencies; ok
	    $freqs{ "\U$_" } = $$ref{$_} foreach keys %$ref;        # copy from $ref; ensure capitalization
	} elsif ($sub1 == 0) {  # all counts; ok
	    $freqs{ "\U$_" } = $$ref{$_} / $sum1 foreach keys %$ref;  # copy from $ref; divide by sum; ensure capitalization
	} else {   # mixed? bad
	    die "Cannot mix frequencies and counts in background hash!  &word_exp halting.\n";
	}
    } elsif ($alpha eq 'FREE') {     # calculate bkg as uniform distrib of extant letters?
	my %letters = map {($_=>1)} (split //, $word);
	my $N = scalar keys %letters;
	$freqs{$_} = 1 / $N foreach keys %letters;
    } else {     # assume random bkg frequencies for known alphabet
	my $N = scalar keys %{ $alphabets{$alpha}{required} };
	$N += 2 if ($word =~ /[UO]/ && $alpha eq 'AA');   # add these too
	$freqs{$_} = 1 / $N foreach keys %{ $alphabets{$alpha}{required} };
    }
    $sum += $freqs{$_} foreach keys %freqs;
    die "Background frequencies sum to $sum, not 1!  &word_exp halting.\n" if sprintf("%0.8f",$sum) != 1;
    
    ### QC INCOMING DATA | FIX DEGENERACIES
    my (%lost, %wrong, %degen, %nondeg, %nonstd);
    if ($alpha eq 'FREE') {
	my @lost;
	foreach my $letter (split //, $word) {
	    push @lost, $letter unless $freqs{$letter};   # freestyle still can't have letters without bkg freqs
	}
	my $lost = join ',', @lost;
	die "Background frequency hash missing values for some given letters ($lost)!  &word_exp halting.\n" if @lost;   # freestyle still can't have letters without bkg freqs
	$nondeg{W}{$_}++ foreach (split //, $word);
    } else {
	## test for all required letters
	foreach my $req (keys %{ $alphabets{$alpha}{required} }) {
	    $lost{$req} unless exists $freqs{$req};       # must have values for all required letters
	}
	my $lost = join ', ', (sort keys %lost);
	die "Background frequency hash missing values for required letters: $lost!  &word_exp halting.\n" if $lost;
	## test for unknown or degenerate letters in word
	foreach my $letter (split //, $word) {
	    if ($alphabets{$alpha}{required}{$letter}) {
		$nondeg{W}{$letter}++;   # non-degenerate
	    } elsif ($letter =~ /^[UO]$/ && $alpha eq 'AA') {
		$nonstd{W}{$letter}++;   # non-required non-degenerate amino acids
	    } elsif ($alphabets{$alpha}{allowed}{$letter}) {
		$degen{W}{$letter}++;    # allowed & not non-degenerate = degenerate
	    } else {
		$wrong{W}{$letter} = 1;    # cannot have unknown letters
	    }
	}
	my $wrongletters1 = join ', ', (sort keys %{ $wrong{W} });  # hopefully the string gets 'undef'
	die "Word contains non-$alpha letters: $wrongletters1!  &word_exp halting.\n" if $wrongletters1;
	## test for unknown or degenerate letters in bkg freq hash
	foreach my $letter (keys %freqs) {
	    if ($alphabets{$alpha}{required}{$letter}) {
		$nondeg{B}{$letter} = 1;   # non-degenerate
	    } elsif ($letter =~ /^[UO]$/ && $alpha eq 'AA') {
		$nonstd{B}{$letter} = 1;   # non-required non-degenerate amino acids
	    } elsif ($alphabets{$alpha}{allowed}{$letter}) {
		$degen{B}{$letter} = 1;    # allowed & not non-degenerate = degenerate
	    } else {
		$wrong{B}{$letter} = 1;    # cannot have unknown letters
	    }
	}
	my $wrongletters2 = join ', ', (sort keys %{ $wrong{B} });  # hopefully the string gets 'undef'
	die "Background frequency hash contains non-$alpha letters: $wrongletters2!  &word_exp halting.\n" if $wrongletters2;
	## distribute degenerates in bkg freq hash to non-degenerates
	foreach my $deg (keys %{ $degen{B} }) {  # nothing happens unless %degen2 has data
	    my @members = @{ $convdegen{$alpha}{$deg} };
	    my ($val, $N) = ($freqs{$deg}, scalar @members);
	    delete $freqs{$deg};   # remove degenerate entry
	    $freqs{$_} += $val / $N foreach @members;  # distribute value evenly among possible real letters
	}
	## test for allowable nonstd chars which aren't in bkg hash (DO NOT test degenerates)
	my %nslost;
	foreach my $letter (keys %{ $nonstd{W} }) {
	    $nslost{$letter} = 1 unless $nonstd{B}{$letter};
	}
	my $nslost = join ', ', (sort keys %nslost);  # hopefully the string gets 'undef'
	die "Background frequency hash missing the following nonstandard letters in word: $nslost!  &word_exp halting.\n" if $nslost;
    }

    ### Calculate weighted expectations for the word
    my $stdfreq = my $nonfreq = my $degfreq = 1;
    foreach my $letter (keys %{ $nondeg{W} }) {
	$stdfreq *= $freqs{$letter} foreach (1..$nondeg{W}{$letter});	# weighted product of frequencies for non-degenerates * number of occurrances
    }
    foreach my $letter (keys %{ $nonstd{W} }) {
	$nonfreq *= $freqs{$letter} foreach (1..$nonstd{W}{$letter});	# weighted product of frequencies for non-standards * number of occurrances
    }
    foreach my $letter (keys %{ $degen{W} }) {  # foreach degenerate base
	foreach (1..$degen{W}{$letter}) {       # foreach instance of degenerate base
	    my $dfreq;
	    $dfreq += $freqs{$_} foreach @{ $convdegen{$alpha}{$letter} };  # total frequency for degenerate letter * number of occurrances
	    $degfreq *= $dfreq;
	}
    }
    my $expfreq = $stdfreq * $nonfreq * $degfreq;
    if ($expfreq) {
	my $hmean = 1 / $expfreq;
	return $hmean;   # where we expect to see 1 occurrence every $hmean bp.
    } else {
	print "Expected frequency is zero: did you specify a bkg frequency of zero for any letters in the word?\n";
	return 0;
    }
}



sub revcomp {
    
    ## reverse-complement DNA/RNA will full degeneracy/masking support
    
    my $SEQ = shift;
    $SEQ = $$SEQ if $SEQ =~ /SCALAR/;   # convert references
    ($SEQ = reverse $SEQ) =~ tr/ACGTURYSWKMHDVBNacgturyswkmhdvbn/TGCAAYRSWMKDHBVNtgcaayrswmkdhbvn/;
    return \$SEQ;   # return reference
}



sub mean {
    
    ## returns the mean of an array
    
    my $ref = shift;
    my $sum;
    $sum += $_ foreach @$ref;
    $sum /= scalar @$ref;
    return $sum;
}



sub median {
    
    ## returns the median of an array
    
    my $ref = shift;
    my $N = scalar @$ref;
    my $med;
    if ($N % 2 == 1) {
	my $mid = ($N - 1) / 2;
	$med = (sort {$a <=> $b} @$ref)[$mid];  # direct median
    } else {
	my $midB = $N / 2;
	my $midA = $midB - 1;
	my @meds = (sort {$a <=> $b} @$ref)[$midA,$midB];  # straddle median
	$med = ($meds[0] + $meds[1]) / 2;
    }
    return $med;
}



sub stdev {
    
    ## returns the standard deviation of an array
    
    my $ref = shift;
    my $N = scalar @$ref;
    my $avg = mean($ref);
    my $sumsqr;
    $sumsqr += ($_-$avg)**2 foreach @$ref;
    my $sd = sqrt( $sumsqr/$N );
    return $sd;
}




sub numjust {
    
    ## justifies a number by adding leading zeroes
    
    my ($num, $width) = @_;
    my $len = length($num);
    if ($len < $width) {
	my $spacer = 0 x ($width - $len);
	$num = "$spacer$num";
    }
    return $num;
}



sub blockify {
    
    ## breaks a sequence into lines of length $N; e.g. for fastas
    
    my @data = @_;
    my $SEQ = $data[0];
    my $WIDTH = $data[1] ? $data[1] : 50;  # default width
    my (@lines, $start);
    my $blocks = length($SEQ) / $WIDTH;
    $blocks++ if (length($SEQ) % $WIDTH != 0);
    foreach (1..$blocks) {
	push @lines, substr($SEQ, $start, $WIDTH);
	$start += $WIDTH;
    }
    my $seqblock = join "\n", @lines;
    return \$seqblock;
}



sub unique {
    
    ## uniques an array
    
    my @array = @{ $_[0] };
    my (@unique, %already);
    foreach (@array) {
	push @unique, $_ unless $already{$_};  # maintains input order w/o using a sort step
	$already{$_} = 1;
    }
    return \@unique;
}



sub chrsort {
    
    ## sorts a list of chromosome names the way I prefer to sort them.  
    ## Specifically designed for UCSC chromosome naming conventions; will also work with some others.
    ## 1. numeric first (increasing numeric sort), 
    ## 2. alpha second (increasing alpha sort), 
    ## 3. randoms third (same order as their non-random counterparts) unless $interleave = 1 (see below)
    ## 4. scaffolds fourth (increasing numeric sort ignoring prefix); multiple sets are sorted increasing by prefix.
    ##    Tests for prefixes are of 2 types: "prefix_else" and "prefix####", the latter ignoring the prefix "chr"
    ##    It takes 5 entries with the same prefix to flag that prefix as indicating a scaffold (or something like it)
    
    # $dataref = chr names array ref
    # $interleave: 1 = sort random chrs next to their nonrandom counterparts; 0 = sort them after the canonical chromosomes
    # $scaflimit: threshold (>=) for prefix prevalence to trigger prefix=scaffold decision
    my @data = @_;
    my ($dataref, $interleave) = @data[0,1];
    my $scaflimit = $data[2] ? $data[2] : 5;
    my (%sets, %prefixes, %isscaf, %romanchrs, @final, $dscount, $drosophila);
    my @dros_strings = qw/ 2L 2R 3L 3R Het /;
    my @drosorder1 = qw/ 2 3 4 mitochondrion M U X Y /;   # fixed order for drosophila chr name part 1 (mito takes precedence over M)
    my @drosorder2 = qw/ L R /;             # fixed order for drosophila chr name part 2
#    my $Bclass = '\d._-';
    my $troubleshoot = 0;

    ## scaffold prefix, Roman, Drosophila tests
    foreach my $chr (@$dataref) {
	my $chrflag = $chr =~ /^chr/ ? 1 : 0;
	(my $ext = $chr) =~ s/^chr//;
	if ($ext ne 'M' && isroman($ext)) {  # roman numeral, but NOT 'M' THIS IS MITO
	    $romanchrs{$chr} = arabic($ext);
	    next;
	} else {
	    foreach my $dstring (@dros_strings) {
		$dscount++ if $chr =~ /$dstring/i;
#		$dscount += 3 if $chr =~ /_mitochondrion_genome/;  # flybase chrM string; guaranteed drosophila
	    }
	}
#	$prefixes{$1}++ if $ext =~ /^(\D+)[$Bclass]+$/;
	$prefixes{$1}++ if $ext =~ /^(\D+)\d/;
    }
    if ($dscount >= 3) {
	$drosophila = 1;   # match any 3 drosophila strings, you are considered drosophila
    } else {
	%romanchrs = () if (scalar (keys %romanchrs) <= 2);  # have to have > 2 romans to have any; otherwise just getting chrX and/or chrM
    }
    foreach my $prefix (keys %prefixes) {
	$isscaf{$prefix} = 1 if $prefixes{$prefix} >= $scaflimit;
    }

    ## chr breakout
    if ($drosophila) {
	
	my (%patterns, @search, @display, %drosordered);     # chrom arrangement templates
	if ($interleave) {
	    %patterns = ('numeric' => '^\d+[LR]?', 'alpha' => '^[MUXY]');
	    @search = @display = qw/ numeric alpha /;
	} else {
	    %patterns = ('numeric' => '^\d+[LR]?', 'numHet' => '^\d+[LR]?H', 'alpha' => '^[MUXY]', 'alphaHet' => '^[MUXY][LR]?H');
	    @search = qw/ numHet numeric alphaHet alpha /;     # Hets search first, else all go non-het
#	    @display = qw/ numeric numHet alpha alphaHet /;    # Hets display second in types
	    @display = qw/ numeric alpha numHet alphaHet /;    # Hets display last
	}
	$drosordered{$_} = [] foreach keys %patterns;
	
	foreach my $chr (@$dataref) {
	    (my $ext = $chr) =~ s/^chr//;
	    my $matched;
	    if ($chr =~ /mitochondrion_genome/) {   # flybase chrM
		$matched = 'alpha';
	    } else {
		foreach my $group (@search) {
		    next if $matched;  # already defined
		    $matched = $group if $ext =~ /$patterns{$group}/i;
		}
	    }
	    if ($matched) {
		push @{ $drosordered{$matched} }, $chr;
	    } else {
		print "Failed to classify '$chr'!\n";
	    }
	}
	
	foreach my $group (@display) {
	    my $flag = $troubleshoot ? "\t$group" : "";
	    push @final, "$_$flag" foreach sort @{ $drosordered{$group} };
	}
	
    } else {

	foreach my $chr (@$dataref) {
	    my $chrflag = $chr =~ /^chr/ ? 1 : 0;
	    (my $ext = $chr) =~ s/^chr//;
	    my $scaffold;
	    foreach my $prefix (keys %isscaf) {
		$scaffold = $prefix if $ext =~ s/^$prefix//;
	    }
	    if ($scaffold) {
		$sets{scaf}{$scaffold}{$ext}{$chr} = 1;
	    } elsif ($romanchrs{$chr}) {
		# do nothing; already stored in %romanchrs
	    } elsif ($chrflag && $ext =~ /^(\d+)(\D+)/) {  
		# chromosome, but extension starts with digits then nondigits?  Add it the numeric randoms.
		# TEST BEFORE "random" since this takes precedence!!!
		my ($ext2, $etc) = ($1, $2);
		$interleave ? ($sets{std}{numeric}{$ext2}{$chr} = 1) : ($sets{rand}{numeric}{$ext2}{$chr} = 1);
	    } elsif ($ext =~ s/_random$//) {
		if ($ext =~ /\D/) {
		    $interleave ? ($sets{std}{alpha}{$ext}{$chr} = 1) : ($sets{rand}{alpha}{$ext}{$chr} = 1);
		} else {
		    $interleave ? ($sets{std}{numeric}{$ext}{$chr} = 1) : ($sets{rand}{numeric}{$ext}{$chr} = 1);
		}
	    } elsif ($ext =~ /^(\d+)(\D+)/) {  
		# not a chrom or scaffold, but extension starts with digits then nondigits?  Add it the 'odds': list after standards but before randoms.
		# TEST AFTER "random" since "random" takes precedence!!!
		$sets{odd}{alpha}{$ext}{$chr} = 1;   # alpha only; no odd numeric
	    } elsif ($ext =~ /\D/) {
		$sets{std}{alpha}{$ext}{$chr} = 1;
	    } else {
		$sets{std}{numeric}{$ext}{$chr} = 1;
	    }
	}

        ## sort & compile
#       @final = ('ORIGINAL', @$dataref, "\n") if $troubleshoot;
	my $flag = $troubleshoot ? "\troman" : "";
	push @final, "$_$flag" foreach (sort { $romanchrs{$a} <=> $romanchrs{$b} } keys %romanchrs);
	foreach my $type (qw/ std odd rand /) {
	    if ($sets{$type}{numeric}) {
		my $flag = $troubleshoot ? "\t$type.numeric" : "";
		foreach my $ext (sort {$a <=> $b} (keys %{ $sets{$type}{numeric} })) {
		    push @final, "$_$flag" foreach (sort keys %{ $sets{$type}{numeric}{$ext} });
		}
	    }
	    if ($sets{$type}{alpha}) {
		my $flag = $troubleshoot ? "\t$type.alpha" : "";
		foreach my $ext (sort (keys %{ $sets{$type}{alpha} })) {
		    push @final, "$_$flag" foreach (sort keys %{ $sets{$type}{alpha}{$ext} });
		}
	    }
	}
	foreach my $prefix (sort keys %isscaf) {
	    if ($sets{scaf}{$prefix}) {
		my $flag = $troubleshoot ? "\tscaf.$prefix" : "";
		foreach my $ext (sort {$a <=> $b} keys %{ $sets{scaf}{$prefix} }) {
		    push @final, "$_$flag" foreach (sort keys %{ $sets{scaf}{$prefix}{$ext} });
		}
	    }
	}
    }
    
    return \@final;
}




#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################
#########################################################    FILE HANDLING SUBROUTINES    #########################################################




sub getfile {
    
    ## reads a text file into a standard format (array-of-arrays or hash-of-arrays) with various options
    
    # $filename = filename | $delim is the file delimiter ("\t", ",", etc) | $chomp = [01] remove terminal [\r\n]+ 
    # $header = integer: capture $header first lines as header | $skip = integer: skip $skip first lines
    # $nfields = integer: return row broken into $nfields fields (as "split /$delim/, $_, $nfields")
    # $keycol = integer: use data in column $keycol (USE 0-BASED!!!) as hash key for storing row (changes output from @-of-@ to %-of-@)
    my ($filename, $delim, $chomp, $header, $skip, $nfields, $keycol) = @ARGV;
    $chomp = 1 if !defined $chomp;
    my (@header, @data, %data);
    open MY_GETFILE_EXCLUSIVE_FILEHANDLE, $filename or die "Can't open '$filename': $!\n";
    while (<MY_GETFILE_EXCLUSIVE_FILEHANDLE>) {
	next if ($skip && $. <= $skip);
	$_ =~ s/[\n\r]+$// if $chomp;
	if ($header && $. <= $header) {
	    push @header, $_;
	    next;
	}
	if ($keycol) {
	    my @temp = split /$delim/, $_;  # gets key from absolute column, not $nfields-split column
	    warn "Key '$temp[$keycol]' already seen!  Overwriting...\n" if exists $data{$temp[$keycol]};
	    $data{$temp[$keycol]} = $nfields ? [split /$delim/, $_, $nfields] : [split /$delim/, $_];
	} else {
	    ($nfields) ? (push @data, [split /$delim/, $_, $nfields]) : (push @data, [split /$delim/, $_]);
	}
    }
    close MY_GETFILE_EXCLUSIVE_FILEHANDLE;
    if ($keycol) {
	return(\%data, \@header);
    } else {
	return(\@data, \@header);
    }
}



sub writefile {
    
    ## writes a text file from a contents block (must be single string, so join "\n" any arrays)
    
    # $filename = filename | $contents = string to write | $mode = '>' or '>>'
    my ($filename, $contents, $mode) = @ARGV;
    $mode = '>' unless $mode;
    open MY_WRITEFILE_EXCLUSIVE_FILEHANDLE, "$mode $filename" or die "Can't write to '$filename': $!\n";
    print MY_WRITEFILE_EXCLUSIVE_FILEHANDLE $contents;
    close MY_WRITEFILE_EXCLUSIVE_FILEHANDLE;
}



sub getfasta {
    
    ## reads a fasta into a hash (keys=headers, values=sequence)
    
    # $filename = filename | $chomp = [01] remove terminal [\r\n]+ (affects sequence block only)
    my ($filename, $chomp) = @ARGV;
    $chomp = 1 if !defined $chomp;
    my ($header, @headers, %data);
    open MY_GETFASTA_EXCLUSIVE_FILEHANDLE, $filename or die "Can't open '$filename': $!\n";
    while (<MY_GETFASTA_EXCLUSIVE_FILEHANDLE>) {
	next if $_ =~ /^#/;
	$_ =~ s/[\n\r]+$// if $chomp;
	if ($_ =~ /^>(.*)/) {
	    $header = $1;
	    push @headers, $header;  # record of header order
	} else {
	    $data{$header} .= $_ if $_ =~ /\S/;  # remove any lines with tabs or spaces
	}
    }
    close MY_GETFASTA_EXCLUSIVE_FILEHANDLE;
    return(\%data, \@headers);
}



sub writefasta {
    
    ## writes a sequence hash (keys=headers, values=sequence) to a fasta file
    
    # $filename = filename | $contents = string to write | $mode = '>' or '>>' | $width = line width
    my ($filename, $dataref, $mode, $width) = @ARGV;
    $mode = '>' unless $mode;
    $width = 50 unless $width;
    open MY_WRITEFASTA_EXCLUSIVE_FILEHANDLE, "$mode $filename" or die "Can't write to '$filename': $!\n";
    print MY_WRITEFASTA_EXCLUSIVE_FILEHANDLE ">$_\n", blockify($$dataref{$_}, $width), "\n" foreach keys %$dataref;
    close MY_WRITEFASTA_EXCLUSIVE_FILEHANDLE;
}



sub get_memedata {
    
    ## takes a meme directory and reads the meme.txt and meme.html files, or returns error
    ## if files exist, returns a hash with a variety of data per motif, current keys:
    ##  WIDTH, NSEQS, LOGLR, EVAL, INFO, ENTRO, CONSENSUS, PSSM (base => array), POS (seq => count), PVAL (pval => count).  
    ##   PVAL is per sequence.  PSSM base-arrays are the width of the motif.
    
    my $memedir = shift;
    my $memetxt = "$memedir/meme.txt";
    my $memehtml = "$memedir/meme.html";
    die "'$memedir/meme.txt' is not readable!\n" unless -e $memetxt;
    die "'$memedir/meme.html' is not readable!\n" unless -e $memehtml;
    
    my %degen = ('AG'=>'R', 'CT'=>'Y', 'CG'=>'S', 'AT'=>'W', 'GT'=>'K', 'AC'=>'M', 'CGT'=>'B', 'AGT'=>'C', 'ACT'=>'H', 'ACG'=>'V', 'ACGT'=>'N');

    my (%mdata, $motif, $reflag, $mbdflag, $mbdflag2, $pssmflag, $dashed, $line, %temp);
    open IN, $memetxt;
    while (<IN>) {
	$_ =~ s/[\n\r]//g;
	if ($reflag && $_ !~ /^-/) {  # given as regular expression; want denegerate consensus
	    my $consensus = my $regex = $_;
	    my @dchars;
	    while ($consensus =~ /\[([A-Z]+)\]/g) {
		my %chars = map {($_=>1)} (split //, $1);
		my $patt = join '', sort keys %chars;
		push @dchars, $degen{$patt};
	    }
	    my $start = my $idx = 0;
	    while ($start != -1) {
		$start = index($consensus, '[', 0);		# keep restarting from 0 because matches get eliminated
		last if $start == -1;
		my $end = index($consensus, ']', $start);
		substr($consensus, $start, $end-$start+1, $dchars[$idx]);
		$idx++;
	    }
	    $mdata{$reflag}{CONSENSUS} = $consensus;
	    $reflag = undef;
	} elsif ($_ =~ /^\s+Motif (\d+) regular expression/) {
	    $reflag = $1;
	} elsif ($_ =~ /\s+Motif (\d+) block diagrams/) {
	    $mbdflag = $1;	# prep for capture
	} elsif ($mbdflag && $_ !~ /^-/) {
	    next if $_ =~ /^SEQUENCE NAME\s/;
	    $mbdflag2 = 1;	# begin capture
	    my ($seq, $pval, $etc) = split /\s+/, $_, 3;
	    $mdata{$mbdflag}{POS}{$seq}++;
	    $mdata{$mbdflag}{PVAL}{$pval}++;
	} elsif ($mbdflag2 && $_ =~ /^-/) {
	    $mbdflag = $mbdflag2 = undef;	# end capture
	} elsif ($_ =~ /^MOTIF\s+(\d+)\s+width =\s+(\d+)\s+sites =\s+(\d+)\s+llr =\s+(\d+)\s+E-value =\s+(\S+)/) {
	    $motif = $1;
	    $mdata{$motif}{WIDTH} = $2;
	    $mdata{$motif}{NSEQS} = $3;
	    $mdata{$motif}{LOGLR} = $4;
	    $mdata{$motif}{EVAL} = $5;
	} elsif ($_ =~ /^\s+Motif (\d+) position-specific probability matrix/) {   # begin matrix capture
	    ($pssmflag, $dashed, $line) = ($1, 0, 0);
	    $mdata{$pssmflag}{PSSM} = {'A'=>[], 'C'=>[], 'G'=>[], 'T'=>[]};  # initialize
	} elsif ($dashed == 2) {	# second dashed line: terminate PSSM capture
	    $dashed = $pssmflag = undef;
	} elsif ($pssmflag) {		# capturing lines
	    if ($_ =~ /^-/) {
		$dashed++;
	    } elsif ($_ =~ /^\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+$/) {
		$mdata{$pssmflag}{PSSM}{A}->[$line] = $1;
		$mdata{$pssmflag}{PSSM}{C}->[$line] = $2;
		$mdata{$pssmflag}{PSSM}{G}->[$line] = $3;
		$mdata{$pssmflag}{PSSM}{T}->[$line] = $4;
		$line++;
	    }
	}
    }
    close IN;
    
    foreach my $motif (sort {$a <=> $b} keys %mdata) {
	$temp{$motif}{obs_seq} = scalar keys %{ $mdata{$motif}{POS} };
	if ($mdata{$motif}{NSEQS} != $temp{$motif}{obs_seq}) {
	    print "get_memedata error: motif $motif: stated sequence count '$mdata{$motif}{NSEQS}' != observed sequence count '$temp{$motif}{obs_seq}'!\n";
	}
    }
    
    my ($inflag, $enflag);
    open IN, $memehtml;
    while (<IN>) {
	$_ =~ s/[\n\r]//g;
	if ($_ =~ /^<big><b><a href=\"#summary_doc\">MOTIF (\d+)</) {
	    $motif = $1;
	} elsif ($_ =~ /^<tr><th><a href=\"#ic_doc\">Information Content</) {
	    $inflag = 1;
	} elsif ($_ =~ /^<tr><th><a href=\"#re_doc\">Relative Entropy</) {
	    $enflag = 1;
	} elsif ($_ =~ /^<tr><td align=\"center\"><b>(\S+) \(bits\)/) {
	    if ($inflag) {
		$mdata{$motif}{INFO} = $1;
		$inflag = undef;
	    } elsif ($enflag) {
		$mdata{$motif}{ENTRO} = $1;
		$enflag = undef;
	    }
	}
    }
    close IN;
    
    return \%mdata;
}



sub get_fimodata {
}



sub get_tomdata {
}



#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################
#########################################################    COORDINATE HANDLING SUBROUTINES    #########################################################




sub UTR_weld {
    
    ## takes 2 hash refs: 1 = CDS {coord => ID}, 2 = UTR {coord => ID}
    ##  optional 3rd arg = stop codon coord, if separated from the terminal CDS.
    ## ALL coords must have format "start\tend".
    ## *** scope should be confined to a SINGLE TRANSCRIPT ***.  This will prevent incorrect UTR-CDS associations.

    ## returned is a hash of exons, including welded UTR-CDS exons and standalone CDS or UTR exons.
    ## return hash tracks exon components: coord => { 'CDS' => { coord => ID }, 'UTR' => { coord => ID } }
    ##  this way if annotations were attached to CDS, UTR components you can drill down into an exon coord and find the original components
    ## ID hash values are not required (can all be 1 or whatever)
    
    my ($CDS, $UTR) = @_[0,1];
    my ($fail, %exons, %already, $STOP);
    $STOP = $_[2] if $_[2];
    
    ## test data completeness
    if (scalar (keys %$CDS) == 0) {
	$fail = 1;
	print "UTR_weld error: no CDS coords supplied!  Nothing to do.\n";
    } elsif (scalar (keys %$UTR) == 0) {
#	$fail = 1;
#	print "UTR_weld error: no UTR coords supplied!  Nothing to do.\n";
    }
    
    ## if $STOP, weld it to appropriate CDS first
    if ($STOP) {
	my $stopadd;
	my ($sstart, $send) = split /\t/, $STOP;
	foreach my $ccoord (keys %$CDS) {
	    my ($cstart, $cend) = split /\t/, $ccoord;
	    if ($cend+1 == $sstart) {        # + strand stop codon
		if ($stopadd) {  # can't attach to > 1 CDS!  Bad CDS list or stop codon position
		    (my $ccoord2 = $ccoord) =~ s/\t/-/;
		    $fail = 1;
		    print "Stop codon may attach to CDS $stopadd or $ccoord2!  Bad CDS list or stop position.\n";
		} else {
		    $stopadd = $ccoord;
		}
		my $newcoord = "$cstart\t$send";
		$$CDS{$newcoord} = $$CDS{$ccoord};
		delete $$CDS{$ccoord};
	    } elsif ($cstart-1 == $send) {   # - strand stop codon
		if ($stopadd) {  # can't attach to > 1 CDS!  Bad CDS list or stop codon position
		    (my $ccoord2 = $ccoord) =~ s/\t/-/;
		    $fail = 1;
		    print "Stop codon may attach to CDS $stopadd or $ccoord2!  Bad CDS list or stop position.\n";
		} else {
		    $stopadd = $ccoord;
		}
		$stopadd = $ccoord;
		my $newcoord = "$sstart\t$cend";
		$$CDS{$newcoord} = $$CDS{$ccoord};
		delete $$CDS{$ccoord};
	    }
	}
	(my $pstop = $STOP) =~ s/\t/-/;
	$stopadd =~ s/\t/-/;
	($stopadd) ? (print "Matched stop codon $pstop to CDS $stopadd.\n") : (print "Failed to add stop codon $pstop!\n");
    }
    
    ## test coords to ensure no overlapping entries anywhere
    my @test = ( (keys %$CDS), (keys %$UTR) );
    my @labels = ( (map {'CDS'} keys %$CDS), (map {'UTR'} keys %$UTR) );
    foreach my $i (0..$#test) {
	my ($start1, $end1) = split /\t/, $test[$i];
	foreach my $j (0..$#test) {
	    next if $i >= $j;
	    my ($start2, $end2) = split /\t/, $test[$j];
	    if ($start2 > $end1 || $start1 > $end2) {
		# no overlap; ok
	    } else {
		$fail = 1;
		(my $testI = $test[$i]) =~ s/\t/-/;
		(my $testJ = $test[$j]) =~ s/\t/-/;
		print "UTR_weld error: $labels[$i] $testI overlaps $labels[$j] $testJ!  Overlaps not allowed.\n";
	    }
	}
    }
    
    ## test to ensure no adjacent CDSs
    my @allCDS = keys %$CDS;
    foreach my $i (0..$#allCDS) {
	my ($start1, $end1) = split /\t/, $allCDS[$i];
	foreach my $j (0..$#allCDS) {
	    next if $i >= $j;
	    my ($start2, $end2) = split /\t/, $allCDS[$j];
	    if ($start2-1 == $end1) {
		$fail = 1;
		(my $coordI = $allCDS[$i]) =~ s/\t/-/;
		(my $coordJ = $allCDS[$j]) =~ s/\t/-/;
		print "UTR_weld error: CDS coords $coordI, $coordJ are adjacent!  Bad transcript structure.\n";
	    } elsif ($start1-1 == $end2) {
		$fail = 1;
		(my $coordI = $allCDS[$i]) =~ s/\t/-/;
		(my $coordJ = $allCDS[$j]) =~ s/\t/-/;
		print "UTR_weld error: CDS coords $coordJ, $coordI are adjacent!  Bad transcript structure.\n";
	    }
	}
    }
    
    unless ($fail) {
	## search for mergeable CDS/UTR sets
	foreach my $ccoord (keys %$CDS) {
	    my ($cstart, $cend) = split /\t/, $ccoord;
	    my (%welds, $ecoord);
	    foreach my $ucoord (keys %$UTR) {
		my ($ustart, $uend) = split /\t/, $ucoord;
		if ($uend+1 == $cstart) {        # 5' UTR junction
		    $welds{$ucoord} = [$ustart, $uend];
		} elsif ($ustart-1 == $cend) {   # 3' UTR junction
		    $welds{$ucoord} = [$ustart, $uend];
		}
	    }
	    $already{CDS}{$ccoord} = 1;
	    if (%welds) {  # CDS+UTR exon
		my @bounds = ($cstart, $cend);            # add CDS coords
		push @bounds, @$_ foreach values %welds;  # add UTR coords
		my $ecoord = join "\t", (sort {$a <=> $b} @bounds)[0,-1];  # terminals of all UTR+CDS coords
		$exons{$ecoord}{CDS}{$ccoord} = $CDS->{$ccoord};
		$exons{$ecoord}{UTR}{$_} = $UTR->{$_} foreach keys %welds;
		$already{UTR}{$_} = 1 foreach keys %welds;
	    } else {       # no appended UTRs -- standalone CDS exon
		$exons{$ccoord}{CDS}{$ccoord} = $CDS->{$ccoord};
	    }
	}
	## search for unmerged UTRs (standalone UTR exons)
	foreach my $ucoord (keys %$UTR) {
	    next if $already{UTR}{$ucoord};
	    $exons{$ucoord}{UTR}{$ucoord} = $UTR->{$ucoord};
	    $already{UTR}{$ucoord} = 1;
	}
	## double-check to see if all input coords are accounted for
	my $input = scalar @test;
	my $output = (scalar keys %{ $already{CDS} }) + (scalar keys %{ $already{UTR} });
	print "UTR_weld error: only $output of $input input coords accounted for!\n" if $input != $output;
    }
    return \%exons;
}



sub bedprep {
    
    ## takes a bed file and prepares it for list matching
    ## specifically, converts it into a hash like $HASH{CHR}{"START\tEND"} or $HASH{CHR}{STR}{"START\tEND"},
    ##  depending if strand-specific matching or not.
    ## really though, only the first 3 fields have to follow BED conventions (or 6 if $strmatch) and the others can be whatever data you want
    
    my ($bed, $strmatch) = @_;
    my ($line, $N, %struct);

    open IN, $bed || die "Cannot open file '$bed': $!\n";
    while (<IN>) {
	$line++;
	next if ($_ =~ /^#/ || $_ =~ /^track/);
	$_ =~ s/[\n\r]+$//;
	$N++;
	my ($chr, $start, $end, @more) = split /\t/, $_;  # @more could be fields up to ($id, $score, $strand) and beyond
	if ($strmatch && !$more[1]) {  # want strand, but no strand data
	    die "Strand matching requested, but file '$bed' line $line missing strand data!\n";
	} elsif ($strmatch) {
	    $struct{$chr}{$more[1]}{"$start\t$end"} = join "\t", @more;
	} else {
	    $struct{$chr}{"$start\t$end"} = join "\t", @more;
	}
    }
    close IN;
#    return (\%struct, $N);
    return \%struct;
}



sub gtfprep {
    
    ## takes a CUFFLINKS-STYLE gtf file and prepares it for list matching (hint: use annotation gtfs in my bowtie-index builds)
    ## specifically, converts it into a hash like $HASH{CHR}{"START\tEND"} or $HASH{CHR}{STR}{"START\tEND"},
    ##  depending if strand-specific matching or not.
    ## also tracks gene structure for exon-level matching.
    ## Use these %struct objects with GENE/EXON MATCHING ROUTINES ONLY.
    
    my ($gtf, $strmatch) = @_;
    my ($line, $N, %structG, %structE, %genecoord);

    open IN, $gtf || die "Cannot open file '$gtf': $!\n";
    while (<IN>) {
	$line++;
	next if ($_ =~ /^#/ || $_ =~ /^track/);
	$_ =~ s/[\n\r]+$//;
	my ($chr, $start, $end, $strand, $annots) = (split /\t/, $_)[0,3,4,6,8];
	my %annot = split /[\s";]+/, $annots;
	my $gene = $annot{gene_id};
	if ($strmatch && !$strand) {  # want strand, but no strand data
	    die "Strand matching requested, but file '$gtf' line $line missing strand data!\n";
	} else {
	    $structE{$gene}{"$start\t$end"} = 1;
	    $genecoord{$gene} = [$chr, $strand];
	}
    }
    close IN;
    
    # add whole-gene boundaries
    foreach my $gene (keys %genecoord) {
	$N++;
	my ($chr, $strand) = @{ $genecoord{$gene} };
	my $start = (split /\t/, (sort {$a <=> $b} keys %{ $structE{$gene} })[0])[0];
	my $end = (split /\t/, (sort {$b <=> $a} keys %{ $structE{$gene} })[0])[1];
	if ($strmatch) {
	    $structG{$chr}{$strand}{"$start\t$end"} = $gene;
	} else {
	    $structG{$chr}{"$start\t$end"} = $gene;
	}
    }
    return (\%structG, \%structE);
}



sub listmatch {
    
    ## takes 2 '%struct' objects and calls the matching routines on all available combinations to be compared
    
    my ($struct1, $struct2, $strmatch, $genematch, $neighbors, $log, $verbose) = @_;
    my (%struct, %alldata);
    my $quiet = $verbose ? 0 : 1;
    
    # Merge both %struct objects to get on the same page
    if ($strmatch) {
	foreach $chr (sort keys %$struct1) {  # this will be shared and %$struct1-specific stuff
	    foreach $str (sort keys %{ $$struct2{$chr} }) {
		$struct{$chr}{$str}{1} = $$struct1{$chr}{$str};
		$struct{$chr}{$str}{2} = $$struct2{$chr}{$str} if exists $$struct2{$chr}{$str};
	    }
	}
	foreach $chr (sort keys %$struct2) {  # this will be %$struct2-specific stuff
	    foreach $str (sort keys %{ $$struct2{$chr} }) {
		next if $struct{$chr}{$str};  # already got it above
		$struct{$chr}{$str}{2} = $$struct2{$chr}{$str};
	    }
	}
    } else {
	foreach $chr (sort keys %$struct1) {  # this will be shared and %$struct1-specific stuff
	    $struct{$chr}{1} = $$struct1{$chr};
	    $struct{$chr}{2} = $$struct2{$chr} if exists $$struct2{$chr};
	}
	foreach $chr (sort keys %$struct2) {  # this will be %$struct2-specific stuff
	    next if $struct{$chr};  # already got it above
	    $struct{$chr}{2} = $$struct2{$chr};
	}
    }
    
    # Now do the matching
    if ($strmatch) {
	print "Matching with strands...\n";
	foreach $chr (sort keys %struct) {
	    foreach $str (sort keys %{ $struct{$chr} }) {
		my $scalar1 += scalar keys %{ $struct{$chr}{$str}{1} };
		my $scalar2 += scalar keys %{ $struct{$chr}{$str}{2} };
		$alldata{N}{1} += $scalar1;
		$alldata{N}{2} += $scalar2;
		if ($scalar1 > 0 && $scalar2 > 0) {
		    my @criteria = ($chr, $str, $genematch, $neighbors);
		    scan_for_matches($struct{$chr}{$str}{1}, $struct{$chr}{$str}{2}, \@criteria, \%alldata);
		} elsif ($scalar2 > 0) {
		    &logreport(" Skipping chr $chr $str: $scalar1 (list 1) vs $scalar2 (list 2)\n", $log, $quiet);	# reporter
		    push @{ $alldata{incomparable}{2} }, "$chr\t$str\t$_\n" foreach sort {$a <=> $b} @{ $struct{$chr}{$str}{2} };
		    $alldata{fates}{skipped}{2} += scalar @list2;
		} else {
		    &logreport(" Skipping chr $chr $str: $scalar1 (list 1) vs $scalar2 (list 2)\n", $log, $quiet);	# reporter
		    push @{ $alldata{incomparable}{1} }, "$chr\t$str\t$_\n" foreach sort {$a <=> $b} @{ $struct{$chr}{$str}{1} };
		    $alldata{fates}{skipped}{1} += scalar @list1;
		}
	    }
	}
    } else {
	print "Matching without strands...\n";
	foreach $chr (sort keys %struct) {
	    my $scalar1 += scalar keys %{ $struct{$chr}{1} };
	    my $scalar2 += scalar keys %{ $struct{$chr}{2} };
	    $alldata{N}{1} += $scalar1;
	    $alldata{N}{2} += $scalar2;
	    if ($scalar1 > 0 && $scalar2 > 0) {
		my @criteria = ($chr, '', $genematch, $neighbors);
		scan_for_matches($struct{$chr}{1}, $struct{$chr}{2}, \@criteria, \%alldata);
	    } elsif ($scalar2 > 0) {
		&logreport(" Skipping chr $chr: $scalar1 (list 1) vs $scalar2 (list 2)\n", $log, $quiet);	# reporter
		push @{ $alldata{incomparable}{2} }, "$chr\t$_\n" foreach sort {$a <=> $b} @{ $struct{$chr}{2} };
		$alldata{fates}{skipped}{2} += scalar @list2;
	    } else {
		&logreport(" Skipping chr $chr: $scalar1 (list 1) vs $scalar2 (list 2)\n", $log, $quiet);	# reporter
		push @{ $alldata{incomparable}{1} }, "$chr\t$_\n" foreach sort {$a <=> $b} @{ $struct{$chr}{1} };
		$alldata{fates}{skipped}{1} += scalar @list1;
	    }
	}
    }
}



sub scan_for_matches {
    
    ## does the actual scanning through two lists to find matches
    
    
    my ($struct1, $struct2, $criteria, $alldata) = @_;
    my ($chr, $str, $genematch, $neighbors) = @$criteria;  # may not have all elements
    my $strmatch = $str ? 1 : 0;
    my $pstr = " $str" if (defined $str);
    my @list1 = sort {$a <=> $b} keys %$struct1;
    my @list2 = sort {$a <=> $b} keys %$struct2;
    my ($N1, $N2) = (scalar(@list1), scalar(@list2));
    &logreport(" Matching chr $chr$pstr: $N1 (list 1) against $N2 (list 2)\n");	# reporter
    ($matching) ? ($lead = "$chr\t$str") : ($lead = $chr);
    my (@X_starts, @X_ends, @X_else);
    foreach my $coord2 (@list2) {
	my ($start, $end) = split /\t/, $coord2, 2;
	push @X_starts, $start;
	push @X_ends, $end;
	push @X_else, $$struct2{$coord2};
    }
    my (%nnspacer, %matched, %cumatch);
    $nnspacer{1} = $nnspacer{2} = "\t";	# for "$distN" entry
    $nnspacer{1} .= "\t" while $$list1[0] =~ /\t/g;	# for remaining entries
    $nnspacer{2} .= "\t" while $$list2[0] =~ /\t/g;

    ## Match List 1 elements to List 2 elements
    $i = $first = 0;
    foreach $coord1 (@list1) {
#       ($matching) ? (print "$chr: $str: $coord1\t:\t", scalar @X_starts, "\n") : (print "$chr: $coord1\t:\t", scalar @X_starts, "\n");	# reporter
	my ($start1, $end1) = split /\t/, $coord1, 2;
	my $else1 = $$struct1{$coord11};
	my $hit;
	foreach $j ($first..$#X_starts) {
#	print "\t\t$j\n";			# reporter
	    my $start2 = $X_starts[$j];
	    my $end2 = $X_ends[$j];
	    if ($end2 < $start1) {		# not there yet
		unless ($matched{2}{$j}) {	# unmatched $j; get neighbors
		    my $dist3 = $end1 - $start2;	# 3' neighbor guaranteed
		                                        # REMEMBER: may have coterminal 3' neighbors -- need to record all (not yet implemented)
		    $$alldata{nn}{ $list2[$j] }{3} = "$dist3\t$list1[$i]";
		    if ($i == 0) {		# no 5' neighbor
			$$alldata{nn}{ $list2[$j] }{5} = $nnspacer{1};
		    } else {
			my $i2 = $i;
			{
			    $i2--;	# move backwards until nearest 5' neighbor is found, or list beginning is reached
			                # REMEMBER: may have coterminal 5' neighbors -- need to record all (not yet implemented)
			    my ($start1B, $end1B, $else1B) = split /\t/, $list1[$i2], 3;
			    if ($end1B < $start2) {
				my $dist5 = $start2 - $end1B;
				$$alldata{nn}{ $list2[$j] }{5} = "$dist5\t$list1[$i2]";
			    } else {
				($i2 > 0) ? (redo) : ($$alldata{nn}{ $list2[$j] }{5} = $nnspacer{1});	# continue 5', or quit
			    }
			}
		    }
		}
		if ($j == $#X_starts && !$hit) {	# if we are here, then list 2 ran out before list1 (5' neighbor only)
		    my $dist5 = $start1 - $end2;
		    $$alldata{nn}{ $list1[$i] }{5} = "$dist5\t$list2[$j]";
		    $$alldata{nn}{ $list1[$i] }{3} = $nnspacer{2};
		}
		next;
	    } elsif ($start2 > $end1) {		# gone too far
		unless ($hit) {			# unmatched $i; get neighbors
		    my $dist3 = $start2 - $end1;	# 3' neighbor guaranteed
		                                        # REMEMBER: may have coterminal 3' neighbors -- need to record all (not yet implemented)
		    $$alldata{nn}{ $list1[$i] }{3} = "$dist3\t$list2[$j]";
		    if ($j == 0) {		# no 5' neighbor
			$$alldata{nn}{ $list1[$i] }{5} = $nnspacer{2};
		    } else {
			my $j2 = $j;
			{
			    $j2--;	# move backwards until nearest 5' neighbor is found, or list beginning is reached
			                # REMEMBER: may have coterminal 5' neighbors -- need to record all (not yet implemented)
			    my ($start2B, $end2B, $else2B) = split /\t/, $list2[$j2], 3;
			    if ($end2B < $start1) {
				my $dist5 = $start1 - $end2B;
				$$alldata{nn}{ $list1[$i] }{5} = "$dist5\t$list2[$j2]";
			    } else {
				($j2 > 0) ? (redo) : ($$alldata{nn}{ $list1[$i] }{5} = $nnspacer{2});	# continue 5', or quit
			    }
			}
		    }
		}
		last;
	    } else {				# match of some kind
		$first = $j unless $hit;	# next time, start at the first match position
		$hit++;
#				print "\t$list2[$j]\n";	# reporter
		$matched{1}{$i}++;
		$matched{2}{$j}++;
		matchmaker($start1, $end1, $start2, $end2, $lead, $coord1, $list2[$j], $alldata, 0);
	    }
	}
	$i++;
    }

    ## Neighborize all trailing unmatched List 2 elements
    if ($neighbors) {
	$i--;	# was 1 too far
	foreach $j ($first..$#X_starts) {
	    unless ($matched{2}{$j}) {
		my ($start1, $end1, $else1) = split /\t/, $list1[$i], 3;
		my ($start2, $end2, $else2) = split /\t/, $list2[$j], 3;
		my $dist5 = $start2 - $end1;
		$$alldata{nn}{ $list2[$j] }{5} = "$dist5\t$list1[$i]";	# if unmatched $j, get 5' neighbors (since there are no 3's)
		                                                        # REMEMBER: may have coterminal 5' neighbors -- need to record all (not yet implemented)
	    }
	}
    }

    ## Tally stuff and store outputs
    $cumatch{$_} += scalar (keys %{ $matched{$_} }) foreach (1..2);
    if (scalar @X_temp > 0) {
	foreach (@X_temp) {
	    if ($exactonly) {
		($_ =~ /Exact$/) ? (push @{ $$alldata{output} }, "$_\n") : (push @{ $$alldata{mismatched} }, "$_\n");
	    } else {
		push @{ $$alldata{output} }, "$_\n";
	    }
	}
    }
    foreach $i (0..$#$list1) {
	unless (exists $matched{1}{$i}) {
	    my $string = "$lead\t$list1[$i]";
	    my $string2 = $neighbors ? "$string\t\|\t$$alldata{nn}{ $list1[$i] }{5}\t\|\t$$alldata{nn}{ $list1[$i] }{3}\n" : "$string\n";
	    push @{ $$alldata{unmatched}{1} }, $string2;    
	    $unmatched{1}++;
	}
    }
    foreach $j (0..$#$list2) {
	unless (exists $matched{2}{$j}) {
	    my $string = "$lead\t$list2[$j]";
	    my $string2 = $neighbors ? "$string\t\|\t$$alldata{nn}{ $list2[$j] }{5}\t\|\t$$alldata{nn}{ $list2[$j] }{3}\n" : "$string\n";
	    push @{ $$alldata{unmatched}{2} }, $string2;    
	    $unmatched{2}++;
	}
    }
}



sub matchmaker {
    ($start1, $end1, $start2, $end2, $lead, $element, $lr2j, $alldata, $tol) = @_;
    my ($dist, $flag);
    my ($mid1, $mid2) = (int($start1 + ($end1 - $start1) / 2), int($start2 + ($end2 - $start2) / 2));
    if ($start1 == $start2 && $end1 == $end2) {           # Exact match
	($dist, $flag) = ($end1 - $start1 + 1, 'Exact');
    } elsif ($start1 >= $start2 && $end1 <= $end2) {      # Eclipse 1 (eclipse of 1 by 2)
	$expansion = ($end1 - $start1 + 1) / ($end2 - $start2 + 1);
	($dist, $flag) = ($end1 - $start1 + 1, 'Eclipse 1');
    } elsif ($start1 <= $start2 && $end1 >= $end2) {      # Eclipse 2 (eclipse of 2 by 1)
	$expansion = ($end1 - $start1 + 1) / ($end2 - $start2 + 1);
	($dist, $flag) = ($end2 - $start2 + 1, 'Eclipse 2');
    } elsif ($end1 >= $start2 && $end1 <= $end2) {        # Overlap 1 (1 is 5'-most)
	$shift = ($end1 - $start2 + 1) / ($end2 - $start2 + 1);
	($dist, $flag) = ($end1 - $start2 + 1, 'Overlap 1');
    } elsif ($start1 >= $start2 && $start1 <= $end2) {    # Overlap 2 (2 is 5'-most)
	$shift = ($end2 - $start1 + 1) / ($end2 - $start2 + 1);
	($dist, $flag) = ($end2 - $start1 + 1, 'Overlap 2');
    } else {					# problem
	push @errors, "$lead\t$element\tERROR\t$lr2j\n";
    }
    push @X_temp, "$lead\t$element\tmatched\t$lr2j\t:\t$dist\t$flag";
    $mtypes{$flag}++;
}



sub genematch {
    
    my ($matchchr, $strand, $tcount) = @_;
    my (@gene_list, @tag_list, %genespace, %sites, $sign, %hits);
    
    &logreport("Sorting site lists...\n", $logfile) if $GLOBAL{PARAM}{VERBOSE};
    
    if ($strand == 1) {
	%genespace = %{ $GLOBAL{GTEDATA}{pgenespace}{$matchchr} } if $GLOBAL{GTEDATA}{pgenespace}{$matchchr};
	%sites = %{ $TEMP{psites}{$matchchr} } if $TEMP{psites}{$matchchr};
	$sign = '+';
    } else {
	%genespace = %{ $GLOBAL{GTEDATA}{ngenespace}{$matchchr} } if $GLOBAL{GTEDATA}{ngenespace}{$matchchr};
	%sites = %{ $TEMP{nsites}{$matchchr} } if $TEMP{nsites}{$matchchr};
	$sign = '-';
    }
    
    @gene_list = sort {$a <=> $b} (keys %genespace);
    @tag_list = sort {$a <=> $b} (keys %sites);
    &logreport("$sign strand: ".scalar @gene_list." genes, ".scalar @tag_list." tags.\nMatching...\n", $logfile) if $GLOBAL{PARAM}{VERBOSE};
    
    ### Assign all raw tags to features
    my ($startgene, $prev_tstart, $gcount, $icount, $loopcount) = (0, 0, 0, 0, 0);
    my $lastgene = scalar @gene_list - 1;
    my $first = undef;
    foreach my $site (@tag_list) {

	my ($tstart, $tend, $tag);
	unless (($tstart, $tend, $tag) = split /\t/, $site) {
	    &logreport("Cannot split $site!\n", $logfile);
	    die;
	}
	my ($matchsite);
	($strand == 1) ? ($matchsite = "$matchchr\t$tstart\t$strand") : ($matchsite = "$matchchr\t$tend\t$strand");
	push @{ $REPORTS{errors} }, "Sites out of order! $tstart $prev_tstart\n" if ($tstart < $prev_tstart);
	$prev_tstart = $tstart;

	### Locate all genes associated with given tag
	my ($matchcount, $i) = (0, 0);
	foreach my $i ($startgene..$lastgene) {

	    $loopcount++;			# this can be used as a reporter (gives number of loops required to match all tags against all genes; number varies on how you scope it)
	    my ($gstart, $gend, $gid);
	    unless (($gstart, $gend, $gid) = split /\t/, $gene_list[$i]) {
		&logreport("Cannot split $gene_list[$i]!\n", $logfile);
		die;
	    }
	    my @EMdata = ($gstart, $gend, $gid, $tstart, $tend, $tag, $matchsite, $site);	# special array with variables that need to be ferried out to &exonmatch

	    if ($gend < $tstart) {					# not there yet
		next;

	    } elsif ($gstart > $tend) {				# too far
		if ($matchcount > 0) {
		    $startgene = $first;			# begin next search at last search's first match...
		} else {
		    $hits{$site}++;
		    $icount++;
		    $startgene = $i;			# ...unless there was no match (then wait for tags to catch up)
		    $TAGDATA{$tag}{SITES}{I}{$matchsite}++;	# purely intergenic tag -- NOTE DIFFERENT HASH ORDER!
		}
		last;

	    } elsif ($tstart <= $gstart && $tend >= $gend) {	## eclipse (of gene by tag)
		$matchcount++;					# number of matches for this site
		$gcount++;
		$hits{$site}++;
		$first = $i if ($matchcount == 1);		# because only the first gets to be $first
		$TAGDATA{$tag}{SITES}{$matchsite}{$gid}{I}++ if ($tstart < $gstart || $tend > $gend);  # because this tag is partially intergenic
		&exonmatch('E', \@EMdata);

	    } elsif ($tstart <= $gstart) { 				## head overlap only
		$matchcount++;	
		$gcount++;
		$hits{$site}++;
		$first = $i if ($matchcount == 1);
		$TAGDATA{$tag}{SITES}{$matchsite}{$gid}{I}++ if ($tstart < $gstart);
		&exonmatch('H', \@EMdata);

	    } elsif ($tend >= $gend) {	 			## tail overlap only
		$matchcount++;	
		$gcount++;
		$hits{$site}++;
		$first = $i if ($matchcount == 1);
		$TAGDATA{$tag}{SITES}{$matchsite}{$gid}{I}++ if ($tend > $gend);
		&exonmatch('T', \@EMdata);

	    } elsif ($tstart >= $gstart && $tend <= $gend) { 	## containment (of tag by gene)
		$matchcount++;	
		$gcount++;
		$hits{$site}++;
		$first = $i if ($matchcount == 1);
		&exonmatch('C', \@EMdata);

	    } else {
		push @{ $REPORTS{errors} }, "Failed to place tag! $matchchr $strand $tag\n";
	    }
	}

	unless (exists $hits{$site}) {				# tags are left over after the last gene; = trailing intergenic
	    $hits{$site}++;
	    $TAGDATA{$tag}{SITES}{I}{$matchsite}++;		# purely intergenic tag -- NOTE DIFFERENT HASH STRUCTURE!
	    $icount++;
	}
    }

    &logreport("$gcount gene hits, $icount intergenic hits recorded for $tcount tag sites.\n", $logfile) if $GLOBAL{PARAM}{VERBOSE};
    if (scalar @tag_list > scalar (keys %hits)) {	# if for some reason some are still unaccounted for???
	push @{ $REPORTS{errors} }, scalar @tag_list." sites > ".scalar (keys %hits)." hits!\n";
    }

}





sub exonmatch {

    ## Set up exon-matching objects
    my ($geneoverlap, $EMdataref) = @_;	# $matchdata is ref to @EMdata in subroutine above (&genematch)
    my ($gstart1, $gend1, $gid1, $tstart1, $tend1, $tag1, $matchsite1, $site1) = @$EMdataref;
    my (%ENhits);
    my @exon_list = sort {$a <=> $b} (keys %{ $GLOBAL{GTEDATA}{gexonpool}{$gid1} });
##	print "Matching to ".scalar (@exon_list)." exons in gene $gid.\n" if $verbose;	# you don't want this... (normally)
    my ($startexon, $prev_estart, $eloopcount, $ecount, $intcount) = (0, 0, 0, 0, 0);
    my $lastexon = scalar (@exon_list) - 1;
    my $efirst = undef;

    ## Now match to exons within gene
    my ($ematchcount, $eloopcount, $j) = (0, 0, 0);
    foreach my $j ($startexon..$lastexon) {

	$eloopcount++;
	my ($estart, $eend, $eid);
	unless (($estart, $eend, $eid) = split /\t/, $exon_list[$j]) {
	    &logreport("Cannot split $exon_list[$j]!\n", $logfile);
	    die;
	}

	if ($eend < $tstart1) {				# not there yet
	    next;

	} elsif ($estart > $tend1) {			# too far
	    if ($ematchcount > 0) {
		$startexon = $efirst;		# begin next search at last search's first match...
	    } else {
		$ENhits{$site1}++;		# how many hits this site got in the routine -- whether to exons or introns
		$TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{N}++;	# purely intronic tag
		$GLOBAL{GTEDATA}{found_tags}{$gid1}{NC}++;
		$intcount++;
		$startexon = $j;		# ...unless there was no match (then wait for tags to catch up)
	    }
	    last;

	} elsif ($tstart1 <= $estart && $tend1 >= $eend) {	## eclipse (of exon by tag)
	    $ematchcount++;					# number of matches for this site
	    $ecount++;
	    $ENhits{$site1}++;
	    $efirst = $j if ($ematchcount == 1);		# because only the first gets to be $efirst
	    $TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{E}++;		# for the exonic component
	    $TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{N}++;		# for the intronic component
	    $GLOBAL{GTEDATA}{found_tags}{$gid1}{NC}++;

	} elsif ($tstart1 <= $estart) { 				## head overlap only
	    $ematchcount++;
	    $ecount++;
	    $ENhits{$site1}++;
	    $efirst = $j if ($ematchcount == 1);
	    $TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{E}++;
	    $TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{N}++ unless (exists $GLOBAL{GTEDATA}{t5exons}{$eid}); # head overlap on 5' exon = intergenic! (= already accounted)
	    $GLOBAL{GTEDATA}{found_tags}{$gid1}{NC}++;

	} elsif ($tend1 >= $eend) {	 			## tail overlap only
	    $ematchcount++;
	    $ecount++;
	    $ENhits{$site1}++;
	    $efirst = $j if ($ematchcount == 1);
	    $TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{E}++;
	    $TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{N}++ unless (exists $GLOBAL{GTEDATA}{t3exons}{$eid}); # tail overlap on 3' exon = intergenic! (= already accounted)
	    $GLOBAL{GTEDATA}{found_tags}{$gid1}{NC}++;

	} elsif ($tstart1 >= $gstart1 && $tend1 <= $gend1) { 	## containment (of tag by exon)
	    $ematchcount++;
	    $ecount++;
	    $ENhits{$site1}++;
	    $efirst = $j if ($ematchcount == 1);
	    $TAGDATA{$tag1}{SITES}{$matchsite1}{$gid1}{E}++;		# purely exonic tag
	    ## no $GLOBAL{GTEDATA}{found_tags}{$gid1}{C} here because all exons are assumed to have been seen already in transcriptome

	} else {
	    push @{ $REPORTS{errors} }, "Failed to place tag within gene! $gid1 $tag1\n";
	}
    }

    push @{ $REPORTS{errors} }, "0 exon matches reported for gene $gid1, site $site1, tag $tag1!\n" if (scalar (keys %ENhits) == 0);
}

1;
