#!/n/local/stage/perlbrew/perlbrew-0.43/perls/perl-5.16.0/bin/perl
use Data::Dumper;
use Roman;
use strict;

##!/n/site/inst/Linux-x86_64/sys/bin/perl

### GENERIC FUNCTIONS:
## timestamp        : three different types of platform-independent timestamps
## runtime          : takes 2 timestamps and returns elapsed time, in various scales
## monthconv        : convert month names to numbers and vice versa
## benchmark        : platform-independent system resource stamps
## logreport        : send a message to screen & logfile simultaneously
## numjust          : pad a digit with leading zeroes
## unique           : return unique values from an array, in order of appearance
## pairs2subgraphs  : given a 2-level hash describing a connectivity matrix (e.g. $hash{$key1}{$key2}, preferably reciprocal so that $hash{$key2}{$key1} also always exists), return all subgraphs.

### MATH FUNCTIONS:
## mean
## median
## stdev
## strip_NA

### BIOLOGICAL FUNCTIONS:
## revcomp          : reverse complement for DNA/RNA, with full degeneracy/masking support
## revcomp_regex    : reverse complement a DNA/RNA regular expression, with full degeneracy/masking support
## word_exp         : get the expected frequency of a DNA/RNA/AA word, given background character frequencies
## coding_potential : calculate coding potential for a DNA sequence, using one of several methods    ########## NOT WORKING PROPERLY ##########
## translate        : translate a DNA sequence
## blockify         : convert a sequence to a fasta block
## chrsort          : sort a list of chromosomes the way I prefer.  Special handling for roman chrnums (yeast) and drosophila chrs
## get_memedata     : get a standard hash of data from a MEME run (DOES NOT WORK FOR 4.8.1)
## UTR_weld         : take a hash of CDS entries + hash of UTR entries and return a hash of exons

### UNDER CONSTRUCTION: may have some functionality at this time, or may not
## getfile          : read a tabular text file into an array-of-arrays or hash-of-arrays
## writefile        : write an array, array-of-arrays, or hash-of-arrays back to a tabular text file
## getfasta         : read a fasta file into a hash
## writefasta       : write a fasta file from a hash
## get_fimodata     : get a standard hash of data from a FIMO run (DOES NOT WORK FOR 4.8.1)
## get_tomdata      : get a standard hash of data from a Tomtom run (DOES NOT WORK FOR 4.8.1)
## length_regex     : calculate length of a regular expression, also report if length is fixed or minimum


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
    
    ## takes 2 and returns the difference, in the specified timescale
    ## mode 1 = 'FULL' timestamp from above timestamp() function, e.g. "Fri May 4 2012 14:00:47"
    ## mode 2 = generic [YYYY-MM-DD HH:MM:SS] timestamp, e.g. "2012-05-04 14:00:47"
    ## mode 3 = Unix `date` style timestamp, e.g. "Fri May  4 14:00:47 CDT 2012"
    
    ### requires 2 timestamp('FULL') entries
    my (%timestamps, %timediff, %timeparts, %absolute, $scale, $mode);
    ($timestamps{1}, $timestamps{2}, $scale, $mode) = @_;  # timestamps with format e.g. 'Thu Jun 23 2011 12:57:08'; for $scale see %scales below
    my %scales = map {($_=>1)} qw/ FULL YEAR MONTH DAY HOUR MINUTE SECOND /;  # time scale to return results in; 'FULL' = full breakout
    # note that scale 'MONTH' will return month breakdown as months elapsed + fraction of end-month elapsed
    my %monthdays = (1,31, 2,28, 3,31, 4,30, 5,31, 6,30, 7,31, 8,31, 9,30, 10,31, 11,30, 12,31);
    my %monthsecs = map {($_=>$monthdays{$_}*86400)} keys %monthdays;
    my %monthnums = (1,'Jan', 2,'Feb', 3,'Mar', 4,'Apr', 5,'May', 6,'Jun', 7,'Jul', 8,'Aug', 9,'Sep', 10,'Oct', 11,'Nov', 12,'Dec');
    my %timescales = ('YEAR',31536000, 'DAY',86400, 'HOUR',3600, 'MIN',60, 'SEC',1);  # times in seconds | months variable, so not here
    foreach my $t (1,2) {
	my ($dname, $mname, $mon, $date, $yr, $zone, $hr, $min, $sec);
	if ($mode == 1) {
	    ($dname, $mon, $date, $yr, $hr, $min, $sec) = split /\s+/, $timestamps{$t};
	    $mon = monthconv($mname,2);
	} elsif ($mode == 2) {
	    ($yr, $mon, $date, $hr, $min, $sec) = split /[\s:-]+/, $timestamps{$t};
	} elsif ($mode == 3) {
	    ($dname, $mname, $date, $hr, $min, $sec, $zone, $yr) = split /[\s:]+/, $timestamps{$t};
	}
	$timeparts{$t}{YEAR} = $yr;
	$timeparts{$t}{MONTH} = $mon;
	$timeparts{$t}{DAY} = $date;
	$timeparts{$t}{HOUR} = $hr;
	$timeparts{$t}{MIN} = $min;
	$timeparts{$t}{SEC} = $sec;
	$absolute{$t} = ($yr-1)*31536000 + ($date-1)*86400 + $hr*3600 + $min*60 + $sec;  # fully elapsed years, days only
	$absolute{$t} += $monthsecs{$_} foreach 1..$timeparts{$t}{MONTH};  # add seconds per day of each fully elapsed month in YTD
    }
    my $seconds = $absolute{2} - $absolute{1};  # total seconds elapsed
    
    if ($scale eq 'FULL' || $scale eq 'YEAR') {
	if ($seconds > $timescales{YEAR}) {
	    $timediff{YEAR} = $seconds/$timescales{YEAR};
	    $seconds %= $timescales{YEAR};
	}
    }
    if ($scale eq 'FULL' || $scale eq 'MONTH') {
	my $thismonth = $timeparts{1}{MONTH};
	if ($seconds > $monthsecs{$thismonth}) {  # months begin elapsing from the start-time month
	    my $continue = 1;
	    while ($continue) {	
		$timediff{MONTH}++;
		$seconds -= $monthsecs{$thismonth};
		$thismonth++;
		$continue = 0 if $seconds < $monthsecs{$thismonth};  # stop if fewer seconds exist than do in the next month
	    }
	}
    }
    if ($scale eq 'FULL' || $scale eq 'DAY') {
	if ($seconds > $timescales{DAY}) {
	    $timediff{DAY} = $seconds/$timescales{DAY};
	    $seconds %= $timescales{DAY};
	}
    }
    if ($scale eq 'FULL' || $scale eq 'HOUR') {
#	if ($seconds > $timescales{HOUR}) {
	    $timediff{HOUR} = $seconds/$timescales{HOUR};
	    $seconds %= $timescales{HOUR};
#	}
    }
    if ($scale eq 'FULL' || $scale eq 'MIN') {
	if ($seconds > $timescales{MIN}) {
	    $timediff{MIN} = $seconds/$timescales{MIN};
	    $seconds %= $timescales{MIN};
	}
    }
	
    if ($scale eq 'FULL') {
	my @elapsed;
	push @elapsed, int($timediff{YEAR}), " years" if $timediff{YEAR};
	push @elapsed, int($timediff{MONTH}), " months" if $timediff{MONTH};
	push @elapsed, int($timediff{DAY}), " days" if $timediff{DAY};
	push @elapsed, int($timediff{HOUR}), " hours" if $timediff{HOUR};
	push @elapsed, int($timediff{MIN}), " minutes" if $timediff{MIN};
	push @elapsed, int($timediff{SEC}), " seconds" if $timediff{SEC};
	my $elapsed = join ', ', @elapsed;
	return "Elapsed: $elapsed\n";
    } else {
	return $timediff{$scale};
    }
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
    
    my ($msg, $logfile) = @_;
    $msg =~ s/[\n\r]+$//;	# so you don't have to remember if it needs newlines
    open OUT, ">> $logfile" or warn "\nlogreport cannot append to $logfile: $!\n";
    print OUT "$msg\n";
    close OUT;
    print "$msg\n";
}



sub revcomp {
    
    ## reverse-complement DNA/RNA with full degeneracy/masking support
    
    my $SEQ = shift;
    $SEQ = $$SEQ if $SEQ =~ /SCALAR/;   # convert references
    ($SEQ = reverse $SEQ) =~ tr/ACGTURYSWKMHDVBNacgturyswkmhdvbn/TGCAAYRSWMKDHBVNtgcaayrswmkdhbvn/;
    return \$SEQ;   # return reference
}



sub revcomp_regex {
    
    ## reverse-complement DNA/RNA regular expressions with full degeneracy/masking support
    ## DOES NOT SUPPORT ALL POSSIBLE REGULAR EXPRESSIONS.  Does not use regexp parse trees.  Designed for revcomping simple DNA/RNA motif regexps.
    ## Basic and non-greedy quantifiers, custom classes, and simple groupings are supported, or at least as far as I have tested.
    ## Anchors, escaped chars and class shortcuts, lookarounds, backreferences, logicals, no-memory groupings, and nested groupings are NOT supported.
    
    my $SEQ = shift;
    
    return('') if $SEQ =~ /\(\?[<!=:]/;  # no lookarounds or no-memory groupings
    return('') if $SEQ =~ /(?<!\\)\\(?!\\)/;  # no backrefs or escaped chars
    return('') if $SEQ =~ /\([^\(]*\(/;  # no nested groupings
    return('') if $SEQ =~ /^\^/ || $SEQ =~ /(?<!\\)\$$/;  # no anchors
   
    my $nonparen = '[^\(\)]+';
    my $nonbrack = '[^\[\]]+';
    my $nonenc = '[\\\]?[^\(\)\[\]]';  # aiming to capture escaped single chars, too -- but in general, escaped things not supported
    my $paren = '\('.$nonparen.'\)';
    my $brack = '\['.$nonbrack.'\]';
    my $quants = '(?<![\\\])[\*\+\?]{1,2}';  # NON-ESCAPED quantifiers -- but in general, escaped things not supported
    my $braces = '\{[\d,]+\}';
    my $anytarg = $paren.'|'.$brack.'|'.$nonenc;
    $SEQ = $$SEQ if $SEQ =~ /SCALAR/;   # convert references
    ($SEQ = reverse $SEQ) =~ tr/ACGTURYSWKMHDVBNacgturyswkmhdvbn\[\]\(\)\{\}/TGCAAYRSWMKDHBVNtgcaayrswmkdhbvn\]\[\)\(\}\{/;
    $SEQ =~ s/\?([+*])($anytarg)/$2$1?/g;  # un-reverse minimal quantifiers
    $SEQ =~ s/($quants|$braces)($anytarg)/$2$1/g;  # put quantifiers on the right side of targets
    $SEQ =~ s/\[($nonbrack)\^\]/[^$1]/g;  # put class negators on the left inside of brackets
    $SEQ =~ s/\{(\d*)(,?)(\d+)?\}/{$3$2$1}/g;  # un-reverse quantifier order inside braces
    return \$SEQ;   # return reference
}



sub length_regex {
    
    ##### UNDER CONSTRUCTION #####
    ## calculates the length of a regular expression
    ## also returns whether this length is a fixed length or a minimum length
    ## DOES NOT SUPPORT ALL POSSIBLE REGULAR EXPRESSIONS.  Does not use regexp parse trees.  Designed for revcomping simple DNA/RNA motif regexps.
    ## Basic and non-greedy quantifiers, custom classes, ^$-anchors, and simple groupings are supported, or at least as far as I have tested.
    ## Escaped chars, class shortcuts, and logicals are given SOME support.
    ## Non-^$ anchors, lookarounds, backreferences, no-memory groupings, and nested groupings are NOT supported.
    
    my ($RE, $view) = @_;
    my $nonparen = '[^\(\)]+';
    my $nonbrack = '[^\[\]]+';
    my $nonenc = '[\\\]?[^\(\)\[\]]';  # aiming to capture escaped chars, too
    my $paren = '\('.$nonparen.'\)';
    my $brack = '\['.$nonbrack.'\]';
    my $quants = '(?<!\\\)[\*\+\?]{1,2}';  # NON-ESCAPED quantifiers
    my $braces = '\{[\d,]+\}';
    my $anytarg = $paren.'|'.$brack.'|'.$nonenc;
    my ($minlen, $maxlen, $incr1, $incr2, $inf);
    
    return('NA','NA') if $RE =~ /\(\?[<!=:]/;  # no lookarounds or no-memory groupings
    return('NA','NA') if $RE =~ /\\\d/;  # no backrefs
    return('NA','NA') if $RE =~ /\([^\(]*\(/;  # no nested groupings
    
    $RE =~ s/^\///g;  # remove bouding slashes
    $RE =~ s/\/$//g;  # remove bouding slashes
    $RE =~ s/^\^//;   # remove anchors
    $RE =~ s/(?<!\\)\$$//;   # remove anchors
    
    while ($RE =~ s/($brack)($quants)//) {
	my ($m1, $q) = ($1, $2);
	if ($q =~ /\+/) {
	    $incr1 = $incr2 = 1;
	} else {
	    $incr1 = 0;
	    $incr2 = $q eq '?' ? 1 : 0;
	}
	$inf = 1 unless $q eq '?';
	$minlen += $incr1;
	$maxlen += $incr2;
	print "1. /$m1$q/ = ($incr1, $incr2)\n" if $view;
    }
    while ($RE =~ s/($brack)\{(\d+)(,?)(\d*)\}//) {
	$minlen += $2;
	$maxlen += $4;
	print "2. /$1\{$2$3$4\}/ = ($2, $4)\n" if $view;
    }
    while ($RE =~ s/($brack)//) {
	$minlen++;
	$maxlen++;
	print "3. /$1/ = 1\n" if $view;
    }
    while ($RE =~ s/($paren)($quants)//) {
	my ($m1, $q) = ($1, $2);
	(my $m2 = $m1) =~ s/[\(\)]//g;  # remove matched parens
	my ($min, $max) = (sort {$a <=> $b} map {length($_)} split /\|/, $m2)[0,-1];
	if ($q =~ /\+/) {
	    ($incr1, $incr2) = ($min, $max);
	} else {
	    $incr1 = 0;
	    $incr2 = $q eq '?' ? $max : 0;
	    $inf = 1 if $q =~ /\*/;
	}
	$inf = 1 unless $q eq '?';
	$minlen += $incr1;
	$maxlen += $incr2;
	print "4. /$m1$q/ = ($incr1, $incr2)\n" if $view;
    }
    while ($RE =~ s/($paren)\{(\d+)(,?)(\d*)\}//) {
	my ($m1, $b1, $b2, $b3) = ($1, $2, $3, $4);
	(my $m2 = $m1) =~ s/[\(\)]//g;  # remove matched parens
	my ($min, $max) = (sort {$a <=> $b} map {length($_)} split /\|/, $m2)[0,-1];
	$incr1 = $min*$b1;
	if ($b2) {
	    $incr2 = $max*$b3;   # zeroes out if infinite   # use of '|' inside parens could be dealt with here
	    $inf = 1 unless $b3;
	} else {
	    $incr2 = $incr1;
	}
	$minlen += $incr1;
	$maxlen += $incr2;
	print "5. /$m1\{$b1$b2$b3}/ = ($incr1, $incr2)\n" if $view;
    }
    while ($RE =~ s/($nonenc)($quants)//) {
	my ($m1, $q) = ($1, $2);
	if ($q =~ /\+/) {
	    $incr1 = $incr2 = 1;
	} else {
	    $incr1 = 0;
	    $incr2 = $q eq '?' ? 1 : 0;
	}
	$inf = 1 unless $q eq '?';
	$minlen += $incr1;
	$maxlen += $incr2;
	print "6. /$m1$q/ = ($incr1, $incr2)\n" if $view;
    }
    while ($RE =~ s/($nonenc)\{(\d+)(,?)(\d*)\}//) {
	$incr1 = $2;
	if ($3) {
	    $incr2 = $4;   # zeroes out if infinite
	    $inf = 1 unless $4;
	} else {
	    $incr2 = $incr1;
	}
	$minlen += $incr1;
	$maxlen += $incr2;
	print "7. /$1\{$2$3$4\}/ = ($incr1, $incr2)\n" if $view;
    }
    my $incr;
    $incr += length($_) foreach (split /\\/, $RE);  # add whatever loose chars remain -- split in escaped chars to get rid of backslashes
    $minlen += $incr;
    $maxlen += $incr;
    $maxlen = 'Inf' if $inf;
    print "8. /$RE/ = $incr\n" if $RE && $view;
    print "Final: $minlen, $maxlen\n" if $view;
    
    return($minlen, $maxlen);
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



sub coding_potential {
    
    ## calculates DNA sequence coding potential, using one of several published methods
    
    my %methods = (     # given in order of increasing sensitivity (and runtime)
	'odds' => 1,    # simple nucleotide odds ratio
	'PSF' => 2,     # position-specific nucleotide frequencies
	'CSF' => 3,     # codon structure factor | for methods 1-3 see Nikolaou & Almirantis, J Mol Evol 2004, PMID:15553086
	'CSF+' => 4,    # threshold-based RNY (allows for 'coding' calls, instead of just getting a value)
	'UFM' => 5      # Universal Feature Method | for methods 4-5 see Carels & Frias, Bioinfo Biol Insights 2009, PMID:20140062
	);
    
    my %revord = (0,2, 1,1, 2,0);   # inverted codon order
    my %stopcodons = map {($_=>1)} qw/ TAA TAG TGA /;  # ochre, amber, opal/umber
    
    my ($NTSEQ, $method, $tmode, $threshold) = @_;   # $threshold is 'coding' calling threshold for CSF+, UFM methods.  See defaults below.
    
    die "UFM method not yet fully implemented!  Try 'CSF' instead.\n" if $method eq 'UFM';
    my $mnames = join ', ', sort keys %methods;
    die "Unknown method '$method'!  Must be one of: '$mnames'\n" unless $methods{$method};
    if ($method eq 'UFM' && $tmode != 6) {
	warn "Method 'UFM' requires 6-frame translation: inceasing 'tmode' from $tmode to 6.\n";
	$tmode = 6;
    }
    $threshold = 1 if $method eq 'UFM' && !$threshold;    # see paper, section 'Optimization of classification thresholds'
    $threshold = 75 if $method eq 'CSF+' && !$threshold;  # ditto
    
    my (%FREQ, $sum);
    $FREQ{$_}++ foreach split //, $NTSEQ;
    my $N = length($NTSEQ);
    my $C = $N/3;
    foreach (keys %FREQ) {
	$FREQ{$_} /= $N;
	$sum += $FREQ{$_};
    }
    die "Background frequencies sum to $sum, not 1!  &coding_potential halting.\n" if sprintf("%0.8f",$sum) != 1;
    my $M = scalar keys %FREQ;
    
    my %CODONS = %{ &translate(\$NTSEQ, $tmode, 0, 0, 0, '', 1) };   # $NTSEQ, $tmode, $stats, $keepcase, $orfs, $stopchar, $ascodons
    
    my (%COM, %STOP);
    foreach my $frame (sort keys %CODONS) {
	my (%NUM, %PSF, %AG, %CGA);
	my ($NUMsum, %PSFsum, %AGsum, %CGAsum, @PSFfail);
#	print "$frame: @{ $CODONS{$frame} }\n";
#	print "$N, $C, ", (scalar @{ $CODONS{$frame} }), "\n";
	my $ncodons = scalar @{ $CODONS{$frame} };
	$NUM{$_} += 1/$ncodons foreach @{ $CODONS{$frame} };   # frequency per codon
	$NUMsum += $NUM{$_} foreach keys %NUM;
	die "Numerator codon frequencies sum to $NUMsum, not 1!  &coding_potential halting.\n" if sprintf("%0.8f",$NUMsum) != 1;
	my %DEN = map {($_=>1)} keys %NUM;  # initialize
	
	$STOP{$frame} = 0;  # ensure printable, if printing debugging messages
	unless ($method eq 'odds') {
	    foreach my $codon (@{ $CODONS{$frame} }) {
		my @bases = split //, $codon;
		$STOP{$frame}++ if $stopcodons{$codon};
		if ($method eq 'UFM') {
		    $AG{A} += 1/$ncodons if $bases[0] eq 'A';
		    $AG{G} += 1/$ncodons if $bases[0] eq 'G';
		    $CGA{C} += 1/$ncodons if $bases[0] eq 'C';
		    $CGA{G} += 1/$ncodons if $bases[1] eq 'G';
		    $CGA{A} += 1/$ncodons if $bases[2] eq 'A';
		} else {
		    $PSF{ $bases[$_] }{$_} += 1/$ncodons foreach (0..2);
		}
	    }
	    if ($method eq 'UFM') {
		$AGsum{$_} += $AG{$_} foreach qw/ A G /;
		$CGAsum{$_} += $CGA{$_} foreach qw/ C G A /;
		print "CGAsum: C=$CGAsum{C} G=$CGAsum{G} A=$CGAsum{A} | AGsum: A=$AGsum{A} G=$AGsum{G}\n";
	    } else {
		foreach my $i (0..2) {
		    $PSFsum{$i} += $PSF{$_}{$i} foreach qw/ A C G T N /;
		    push @PSFfail, " $i=$PSFsum{$i}" if sprintf("%0.8f",$PSFsum{$i}) != 1;
		}
		if (@PSFfail) {
		    my $msg = join "\n", ("Codon-position frequencies have failed to sum to 1 at the following positions:", @PSFfail, " &coding_potential halting.\n");
		    die $msg;
		}
	    }
	}
	
	if ($method eq 'UFM') {   # see paper, section 'Scoring the coding potential of ORFs with UFM'
	    my $W = 0.01;
	    $COM{$frame} = ($AG{A}*$AG{G}) / ($CGA{C}*$CGA{G}*$CGA{A}+$STOP{$frame}+$W);
#	    print "$frame: ($AG{A}*$AG{G}) / ($CGA{C}*$CGA{G}*$CGA{A}+$STOP{$frame}+$W) = $COM{$frame}\n";
	} else {
	    foreach my $codon (keys %NUM) {
		my @bases = split //, $codon;
		if ($method eq 'odds') {
		    $DEN{$codon} *= $FREQ{$_} foreach @bases;
		} else {
		    if ($method eq 'PSF') {
			$DEN{$codon} *= $PSF{ $bases[$_] }{$_} foreach (0..2);
		    } elsif ($method =~ /^CSF/) {
			$DEN{$codon} *= $PSF{ $bases[$_] }{ $revord{$_} } foreach (0..2);
		    }
		}
#		print "$frame: $codon: $NUM{$codon} / $DEN{$codon}\n";
#		print "$frame: $codon: $NUM{$codon} / $DEN{$codon} = ", ($NUM{$codon} / $DEN{$codon}), "\n";
		$COM{$frame} += $NUM{$codon} / $DEN{$codon};
	    }
#	    print "$frame: COM = $COM{$frame}\n";
	}
    }
    
    my $maxframe = (sort { $COM{$b} <=> $COM{$a} } keys %COM)[0];
    my $minframe = (sort { $COM{$a} <=> $COM{$b} } keys %COM)[0];
    my $maxcp = $COM{$maxframe};
    my $mincp = $COM{$minframe};
    my $iscoding = 0;
    if ($method eq 'UFM' || $method eq 'CSF+')  {
	$iscoding = $maxcp - $mincp > $threshold ? 1 : -1;
    } elsif ($method eq 'CSF+')  {
	$iscoding = $maxcp > $threshold ? 1 : -1;
    }
    return($maxframe, $maxcp, $minframe, $mincp, $maxcp-$mincp, $iscoding);   # best & worst frames + scores; $iscoding gives 1=TRUE | -1=FALSE | 0=UNTESTED
    
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
    ## RETURNS A SCALAR REFERENCE
    ## block does NOT end with a newline
    
    my @data = @_;   # (sequence, line width)
    my $SEQ = $data[0] =~ /SCALAR/ ? ${$data[0]} : $data[0];
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
    ## Specifically designed for UCSC chromosome naming conventions; probably works with some others.
    ## 1. numeric first (increasing numeric sort), 
    ## 2. alpha second (increasing alpha sort), 
    ## 3. randoms third (same order as their non-random counterparts) unless $interleave = 1 (see below)
    ## 4. scaffolds fourth (increasing numeric sort ignoring prefix); multiple sets are sorted increasing by prefix.
    ##    Tests for prefixes are of 2 types: "prefix_else" and "prefix####", the latter ignoring the prefix "chr"
    ##    It takes 5 entries with the same prefix to flag that prefix as indicating a scaffold (or something like it)
    
    # $dataref = chr names array ref
    # $interleave: 1 = sort random chrs next to their nonrandom counterparts; 0 = sort them after the canonical chromosomes
    # $scaflimit: threshold (>=) for prefix prevalence to trigger prefix=scaffold decision
    
    my ($dataref, $interleave) = @_[0,1];
    my $scaflimit = $_[2] ? $_[2] : 5;
    my (%sets, %prefixes, %isscaf, %romanchrs, @final, $dscount, $drosophila);
    my @dros_strings = qw/ 2L 2R 3L 3R Het /;
    my @drosorder1 = qw/ 2 3 4 mitochondrion M U X Y /;   # fixed order for drosophila chr name part 1 ('mitochondrion' takes precedence over M)
    my @drosorder2 = qw/ L R /;             # fixed order for drosophila chr name part 2
#    my $Bclass = '\d._-';
    my $troubleshoot = 0;

    ## scaffold prefix, Roman, Drosophila tests
    foreach my $chr (@$dataref) {
	my $chrflag = $chr =~ /^chr/ ? 1 : 0;
	(my $ext = $chr) =~ s/^chr//;
	if ($ext ne 'M' && isroman($ext)) {  # roman numeral, BUT NOT 'M' THAT IS MITO
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
		# TEST AFTER "random" since "random takes precedence!!!
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



sub UTR_weld {
    
    ## designed to return whole exons out of bookended CDS + UTR entries which are listed separately, e.g. from Flybase GFF.
    ## *** scope should be confined to a SINGLE TRANSCRIPT ***.  This will prevent incorrect UTR-CDS associations.
    ## takes 2 hash refs: 1 = CDS {coord => ID}, 2 = UTR {coord => ID}
    ##  optional 3rd arg = stop codon coord, if separated from the terminal CDS.
    ## ALL coords must have format "start\tend".
    
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



sub mean {
    
    ## returns the mean of an array
    
    my ($ref, $NArm) = @_;
    my @clean = @{ &strip_NA($ref, $NArm) };
    
    if ($clean[0] eq 'NA') {  # had NAs, and didn't specify removal
	return 'NA';
    } else {
	my $sum;
	$sum += $_ foreach @clean;
	$sum /= scalar @clean;
	return $sum;
    }
}



sub median {
    
    ## returns the median of an array
    
    my ($ref, $NArm) = @_;
    my @clean = @{ &strip_NA($ref, $NArm) };
    
    if ($clean[0] eq 'NA') {  # had NAs, and didn't specify removal
	return 'NA';
    } else {
	my $N = scalar @clean;
	my $med;
	if ($N % 2 == 1) {
	    my $mid = ($N - 1) / 2;
	    $med = (sort {$a <=> $b} @clean)[$mid];  # direct median
	} else {
	    my $midB = $N / 2;
	    my $midA = $midB - 1;
	    my @meds = (sort {$a <=> $b} @clean)[$midA,$midB];  # straddle median
	    $med = ($meds[0] + $meds[1]) / 2;
	}
	return $med;
    }
}



sub stdev {
    
    ## returns the standard deviation of an array
    
    my ($ref, $NArm) = @_;
    my @clean = @{ &strip_NA($ref, $NArm) };
    
    if ($clean[0] eq 'NA') {  # had NAs, and didn't specify removal
	return 'NA';
    } else {
	my $N = scalar @$ref;
	my $avg = mean($ref);
	my $sumsqr;
	$sumsqr += ($_-$avg)**2 foreach @$ref;
	my $sd = sqrt( $sumsqr/$N );
	return $sd;
    }
}



sub strip_NA {
    
    ## returns the input array with no blank elements
    ## $NArm specifies NA-removal status of calling function:
    ## - If '1', strip_NA returns 'NA' if any NAs.
    ## - Otherwise, returns stripped array.
    
    my ($ref, $NArm) = @_;

    my $N1 = scalar @$ref;
    my @clean;
    foreach (@$ref) {
	push @clean, $_ if (defined $_ && $_ ne '');
    }
    my $N2 = scalar @clean;
    if (!$NArm && $N1 > $N2) {  # NAs existed
	return ['NA'];
    } else {
	return \@clean;
    }
}



sub translate {
    
    ## translates a sequence, passed in as scalar ref
    ## options: 
    ## 1. choice of stop-codon character (default '*', also can print { # $ % } for { ochre amber opal })
    ## 2. choice of frames { 1 3 6 }, default 1
    ## 3. include stats with returned seq? { '' Frame Orf } (default '' = no); if 'Frame' generates frame-level stats; if 'Orf' generates ORF-level stats for each max orf per frame
    ## 4. leave AA in lowercase where codons had any lowercase letters (default no)
    ## 5. return codon arrays instead of amino acid sequence?  (intended for &coding_potential)
    
    my (%TRSEQ, %STATS);  # return objects
    
    my %stopcodons = ('TAA'=>'#', 'TAG'=>'$', 'TGA'=>'%');  ## ochre, amber, opal/umber

    my %aminoacids = (
	'TTT' => 'F',	'TTC' => 'F',	'TTA' => 'L',	'TTG' => 'L',
	'TCT' => 'S',	'TCC' => 'S',	'TCA' => 'S',	'TCG' => 'S',
	'TAT' => 'Y',	'TAC' => 'Y',	'TAA' => '*',	'TAG' => '*',
	'TGT' => 'C',	'TGC' => 'C',	'TGA' => '*',	'TGG' => 'W',
	'CTT' => 'L',	'CTC' => 'L',	'CTA' => 'L',	'CTG' => 'L',
	'CCT' => 'P',	'CCC' => 'P',	'CCA' => 'P',	'CCG' => 'P',
	'CAT' => 'H',	'CAC' => 'H',	'CAA' => 'Q',	'CAG' => 'Q',
	'CGT' => 'R',	'CGC' => 'R',	'CGA' => 'R',	'CGG' => 'R',
	'ATT' => 'I',	'ATC' => 'I',	'ATA' => 'I',	'ATG' => 'M',
	'ACT' => 'T',	'ACC' => 'T',	'ACA' => 'T',	'ACG' => 'T',
	'AAT' => 'N',	'AAC' => 'N',	'AAA' => 'K',	'AAG' => 'K',
	'AGT' => 'S',	'AGC' => 'S',	'AGA' => 'R',	'AGG' => 'R',
	'GTT' => 'V',	'GTC' => 'V',	'GTA' => 'V',	'GTG' => 'V',
	'GCT' => 'A',	'GCC' => 'A',	'GCA' => 'A',	'GCG' => 'A',
	'GAT' => 'D',	'GAC' => 'D',	'GAA' => 'E',	'GAG' => 'E',
	'GGT' => 'G',	'GGC' => 'G',	'GGA' => 'G',	'GGG' => 'G');
    
    my ($NTSEQ, $tmode, $stats, $keepcase, $orfs, $stopchar, $ascodons) = @_;
    
    $stats = 0 if $ascodons;   ##########  NO STATS IN THIS MODE  ##########
    
    $NTSEQ = $$NTSEQ if $NTSEQ =~ /SCALAR/;   # convert references
    $NTSEQ =~ s/^\s+//;	# clip leading whitespace
    $NTSEQ =~ s/\s+$//;	# clip trailing whitespace
    $NTSEQ =~ s/(.*)/\U$1/ unless $keepcase;	# all caps
    $NTSEQ =~ s/U/T/g;	# RNA->DNA
    warn "Sequence contains non-ACGTN characters!\n" if $NTSEQ =~ /[^ACGTNacgtn]/;   # nonstandard letters exist
    
    if ($stopchar eq 'names') {
	$aminoacids{$_} = $stopcodons{$_} foreach keys %stopcodons;
    } elsif ($stopchar && $stopchar ne '*') {
	$aminoacids{$_} = $stopchar foreach keys %stopcodons;
    }
    
    if ($stopchar && $stopchar ne 'names') {   # custom-defined stop char names: validate!
	my %AAnonStop = my %AAonlyStop = %aminoacids;
	delete $AAnonStop{$_} foreach keys %stopcodons;  ## amino acids only
	delete $AAonlyStop{$_} foreach keys %AAnonStop;  ## final stop codon names only
	my %AAchar = map {($_=>1)} values %AAnonStop;
	my %crossover;
	foreach (keys %stopcodons) {
	    $crossover{$_} = 1 if exists $AAchar{ $AAonlyStop{$_} };  ## we should NOT have any crossover between stop codon names and amino acid residues
	}
	if ($keepcase) {
	    my %AAcharLC = map { "\L$_" } keys %AAchar;
	    foreach (keys %stopcodons) {
		$crossover{$_} = 1 if exists $AAcharLC{ $AAonlyStop{$_} };  ## (might use lowercase AA characters) + (custom stops could be lowercase) = (ensure no overlap)
	    }
	}
	if (%crossover) {  # bad choice of stop codon names
	    my $fail = join ', ', sort keys %crossover;
	    warn "With the given settings, the following stop characters overlap with legal AA residue characters: '$fail'.\n Switching stop characters to { # $ % }...\n";
	    $aminoacids{$_} = $stopcodons{$_} foreach keys %stopcodons;
	}
    }
    my $stopclass = '[' . (join '', map { $aminoacids{$_} } keys %stopcodons) . ']';   # regex
    
    my (@framelist, %frameseq);
    if ($tmode == 1) {
	@framelist = ('+0');
    } elsif ($tmode == 3) {
	@framelist = ('+0','+1','+2');
    } elsif ($tmode == 6) {
	@framelist = ('+0','+1','+2','-0','-1','-2');
    }
    
    $frameseq{'+0'} = $NTSEQ;
    if ($tmode > 1) {
	my @jsplit = split //, $NTSEQ;
	$frameseq{'+1'} = join '', @jsplit[1..$#jsplit];
	$frameseq{'+2'} = join '', @jsplit[2..$#jsplit];
	if ($tmode == 6) {
	    $frameseq{'-0'} = ${ &revcomp(\$NTSEQ) };
	    my $temp1 = join '', @jsplit[0..$#jsplit-1];
	    my $temp2 = join '', @jsplit[0..$#jsplit-2];
	    $frameseq{'-1'} = ${ &revcomp(\$temp1) };
	    $frameseq{'-2'} = ${ &revcomp(\$temp2) };
	}
    }
    
    foreach my $frame (@framelist) {
	my (@translation, $stops, $ambigs, $much);
	my $seqlen = length($frameseq{$frame});
	my $startmax = int($seqlen/3);
	foreach my $start (0..$startmax) {  		# one for each codon
	    $start *= 3;
	    if ($start == $seqlen) {
		# done
	    } elsif (abs($seqlen-$start) <= 2) {	# reached end of translatable sequence
		my $off = abs($seqlen-$start);
		($seqlen - $start < 0) ? ($much = 'short') : ($much = 'long');
		my $msg = "Sequence was $off base";
		$msg .= "s" if ($off > 1);
		$msg .= " too $much for a full in-frame translation";
		$STATS{$frame}{MSG} = $msg;
	    } else {
		my ($codon, $residue);
		if ($keepcase) {
		    my $Rcodon = substr($frameseq{$frame}, $start, 3);
		    my $Ucodon = $codon = "\U$Rcodon";
		    $residue = $aminoacids{$Ucodon};
		    $stops++ if $stopcodons{$Ucodon};
		    $residue = "\L$residue" if $Rcodon ne $Ucodon;  # raw had lowercase letters
		} else {
		    $codon = substr($frameseq{$frame}, $start, 3);
		    $residue = $aminoacids{$codon};
		    $stops++ if $stopcodons{$codon};
		}
		if ($ascodons) {
		    push @translation, $codon if $residue;  # don't keep codons for untranslatable material
		    print "FAIL: '$codon' = '$residue'\n" unless $residue;
		} elsif ($residue) {
		    push @translation, $residue;
		} else {		# codon failed to translate: probably due to Ns or use of extended IUPAC alphabet
#		    push @errors, "header $HEADER | residue $residue | codon $codon | start $start | frame $frame | len $seqlen | seq $frameseq{$i}\n";
		    push @translation, 'X';
		    $ambigs++;
		}
	    }
	}
	
	if ($ascodons) {
	    $TRSEQ{$frame} = \@translation;
	    next;
	}
	
	my $AASEQ = join '', @translation;
	$TRSEQ{$frame} = \$AASEQ;
	
	if ($stats) {
	    $STATS{$frame}{NTLEN} = $seqlen;
	    $STATS{$frame}{AALEN} = $startmax;
	    $STATS{$frame}{STOPS} = $stops || 0;
	    $STATS{$frame}{AMBIGS} = $ambigs || 0;
	    $STATS{$frame}{AASEQ} = $AASEQ;
	    
	    my @orftest = split /$stopclass+/, $AASEQ;	# split on stops
	    my (%orflens, $maxorfcount);
	    $orflens{ length($_) }{$_} = 1 foreach @orftest;
	    my $orfmax = (sort {$b <=> $a} keys %orflens)[0];
	    $STATS{$frame}{MAXORFAA} = $orfmax;
	    $STATS{$frame}{MAXORFNT} = $orfmax*3;
	    $STATS{$frame}{NMAXORF} = scalar keys %{ $orflens{$orfmax} };
	    
	    if ($stats eq 'Orf') {
		foreach my $orfseq (keys %{ $orflens{$orfmax} }) {	# allows for > 1 longest ORF
		    $maxorfcount++;
		    $STATS{$frame}{ORF}{$maxorfcount}{AALEN} = $STATS{$frame}{MAXORFAA};
		    $STATS{$frame}{ORF}{$maxorfcount}{NTLEN} = $STATS{$frame}{MAXORFNT};
		    $STATS{$frame}{ORF}{$maxorfcount}{T5} = $orfseq eq $orftest[0] ? 1 : 0;    # ORF is 5'-terminal?
		    $STATS{$frame}{ORF}{$maxorfcount}{T3} = $orfseq eq $orftest[-1] ? 1 : 0;   # ORF is 3'-terminal?
		    $STATS{$frame}{ORF}{$maxorfcount}{M5} = $orfseq =~ /^M/i ? 1 : 0;          # ORF starts with an M?
		    $STATS{$frame}{ORF}{$maxorfcount}{S3} = $orfseq =~ /$stopclass$/i ? 1 : 0; # ORF ends with a stop?
		}
	    }
	}
    }
    
    if ($stats) {
	return( \%TRSEQ, \%STATS );
    } else {
	return( \%TRSEQ );
    }
}



sub pairs2subgraphs {
    
    ## extracts subgraphs from a connectivity hash
    ## connectivity hash must be complete (full matrix w/ diagonal, not just one triangle)
    ##  e.g. for any connected pair ($key1,$key2) must contain $hash{$key1}{$key2}, $hash{$key2}{$key1}, $hash{$key1}{$key1}, and $hash{$key2}{$key2}.
    ## subgraphs are returned in an ordered hash-of-arrays (i.e. @{ $subg{$i} } = @keys), where $i is rank, in order of decreasing size.
    
    my ($ref, $verbose) = @_;
    my %conn = %$ref;
    
    # First, expand connectivities for each key until exhausted
    my $itotal = my $iters = 0;
    {
	$iters++;
	my $prev_total = $itotal;
	$itotal = 0;
	if ($verbose) {
	    chomp(my $now = `date`);
	    print " iter $iters: $now | current $prev_total\n";
	}
	foreach my $key1 (keys %conn) {
	    foreach my $key2 (keys %{ $conn{$key1} }) {
		$conn{$key2}{$key1} = 1;
		$conn{$key1}{$_} = 1 foreach keys %{ $conn{$key2} };
	    }
	}
	$itotal += scalar keys %{ $conn{$_} } foreach keys %conn;
	redo if $itotal > $prev_total;
    }
    print "Network converged in $iters iterations: ", `date` if $verbose;
    
    # Second, collect unique key sets
    my (%unique, %subgraphs, $i);
    foreach my $key (keys %conn) {
	my @set = (sort keys %{ $conn{$key} });
	my $label = join ';', @set;
	$unique{$label} = \@set;
    }
    foreach my $label (sort { scalar(@{ $unique{$b} }) <=> scalar(@{ $unique{$a} }) } keys %unique) {
	$i++;
	$subgraphs{$i} = $unique{$label};
    }
    return \%subgraphs;
}

1;
