#!/usr/bin/perl



######  SMART PROCESSING NOT YET READY



### This script is a generalized version of the Solexa interpreter algorithm, which can take two ranked lists of coordinates and return matches, with or without overlaps.
### Element matching is accelerated by position tracking within the arrays, so that distant upstream/downstream elements are not needlessly searched.
### Both lists can have redundant and/or overlapping entries without creating problems.
### Column structure of the lists is based on track file data, i.e. {CHR#\tSTART\tEND\tDATA} or {CHR#\tSTART\tEND\tSTRAND\tDATA}.
### Generalized, these become {GROUP\tSTART\tEND\tDATA} and {GROUP\tSTART\tEND\tSUBGROUP\tDATA}, respectively.
### Use of subgroup matching is optional and can be used if the [strand / subgroup] column is present.
### Overlaps may be reported as such (option 'o'), and return # overlapping bp, or they may be included with matches (option 'm'), or simply ignored (option 'i').
### 
### COMMAND-LINE SUMMARY: > perl Generic_Redundant_List_Matcher.pl [List1 (first list)] [List2 (second list)] [y/n (subgroup matching)] [y/n (smart processing)] [o/m/i (overlap reporting)]
### E.g. 'perl Generic_Redundant_List_Matcher.pl peaks.txt yeast_genes.txt n y o'
### All five options must be used (and used correctly) for command-line execution, else script will die.  If no options given, script will prompt for input.

########################################     SET PARAMETERS     ########################################

if (scalar @ARGV > 0) {
	(scalar @ARGV == 5) ? (($list1, $list2, $match1, $smart1, $over1) = @ARGV) : (die "Inputs not equal to 5!  Cannot proceed.  5 elements required:\nList1 List2 SubgroupMatch?(y/n) SmartProcessing?(y/n) OverlapReport?(o/m/i)\n");
	if ($match1 =~ /^y$/i) {
		$matching = 1;
	} elsif ($match1 =~ /^n$/i) {
	} else {
		die "Subgroup matching choice must be [yY] or [nN] only!\n[yY]: match subgroups within groups\n[nN]: ignore subgroups during matching\n";
		redo;
	}
	if ($smart1 =~ /^y$/i) {
		$smart = 1;
	} elsif ($smart1 =~ /^n$/i) {
	} else {
		die "Smart processing choice must be [yY] or [nN] only!\n[yY]: Smart processing (faster for very large lists)\n[nN]: Simple processing (faster for small lists)\n";
		redo;
	}
	if ($over1 =~ /^o$/i) {
		$overlap = 1;
	} elsif ($over1 =~ /^m$/i) {
	} elsif ($over1 =~ /^i$/i) {
		$overlap = 0;
	} else {
		die "Overlap reporting choice must be [oO], [mM], or [iI] only!\n[oO]: report overlaps as such\n[mM]: treat overlaps as matches\n[iI]: ignore matches that overlap\n";
		redo;
	}
} else {
	print "\nLists must be in 'GROUP\\tSTART\\tEND\\tSUBGROUP\\tELSE\\n' format.\n\nSubgroup matching is optional; field may be omitted ('.' = ignore subgrouping)\n";
	print "\nEnter the first list:\n";
	chomp($list1 = <STDIN>);
	print "\nEnter the second list:\n";
	chomp($list2 = <STDIN>);
	print "\nMatch to subgroup? (y/n)\n";
	{
		chomp($match1 = <STDIN>);
		if ($match1 =~ /^y$/i) {
			$matching = 1;
		} elsif ($match1 =~ /^n$/i) {
		} else {
			print "Enter [yY] or [nN] only!\n";
			redo;
		}
	}
	print "\nUse smart processing to expedite processing of very large lists? (y/n)\n";
	{
		chomp($smart1 = <STDIN>);
		if ($smart1 =~ /^y$/i) {
			$smart = 1;
		} elsif ($smart1 =~ /^n$/i) {
		} else {
			die "Smart processing choice must be [yY] or [nN] only!\n";
			redo;
		}
	}
	print "\nReport overlaps as overlaps, as matches, or ignore? (o/m/i)\n";
	{
		chomp($over1 = <STDIN>);
		if ($over1 =~ /^o$/i) {
			$overlap = 1;
		} elsif ($over1 =~ /^m$/i) {
		} elsif ($over1 =~ /^i$/i) {
			$overlap = 0;
		} else {
			print "Enter [oO], [mM], or [iI] only!\n";
			redo;
		}
	}
}

########################################     PREPARE FOR PROCESSING     ########################################

my (%pools, @output, %mtypes, $lname1, $lname2);
$mtypes{$_} = 0 foreach ('Exact', 'Eclipse 1', 'Eclipse 2', 'Overlap 1', 'Overlap 2');	# ensure printable values
($list1 =~ /\/?([^\/]+)$/) ? ($lname1 = $1) : ($lname1 = $list1);
($list2 =~ /\/?([^\/]+)$/) ? ($lname2 = $1) : ($lname2 = $list2);
my $basename = join '_', ($lname1, $lname2, $match1, $over1);

foreach (1..2) {	# ensure printable values
	$all{$_} = 0;
	$cumatch{$_} = 0;
	$skipped{$_} = 0;
	$unmatched{$_} = 0;
}

open IN1, $list1 || die "Cannot open $list1!\n";
print "Processing $lname1...\n";
while (<IN1>) {
	$all{1}++;
	my ($group, $start, $end, $else);
	$line = $_;
	$line =~ s/[\n\r]//g;
	if ($matching) { 
		my ($group, $start, $end, $sub, $else) = split /\t/, $line, 5 || die "Subgroup matching requires 5 columns of data minimum!\n";
		push @{ $pools{$group}{$sub}{1} }, "$start\t$end\t$else";
	} else {
		my ($group, $start, $end, $else) = split /\t/, $line, 4;
		push @{ $pools{$group}{1} }, "$start\t$end\t$else";
	}
}
close IN1;

open IN2, $list2 || die "Cannot open $list2!\n";
print "Processing $lname2...\n";
while (<IN2>) {
	$all{2}++;
	my ($group, $start, $end, $else);
	$line = $_;
	$line =~ s/[\n\r]//g;
	if ($matching) {
		my ($group, $start, $end, $sub, $else) = split /\t/, $line, 5 || die "Subgroup matching requires 5 columns of data minimum!\n";
		push @{ $pools{$group}{$sub}{2} }, "$start\t$end\t$else";
	} else {
		my ($group, $start, $end, $else) = split /\t/, $line, 4;
		push @{ $pools{$group}{2} }, "$start\t$end\t$else";
	}
}
close IN2;

if ($matching) {
	print "Matching with subgroups...\n";
	foreach $group (sort keys %pools) {
		foreach $sub (sort keys %{ $pools{$group} }) {
			my @list1 = sort {$a <=> $b} @{ $pools{$group}{$sub}{1} };
			my @list2 = sort {$a <=> $b} @{ $pools{$group}{$sub}{2} };
			if (scalar @list1 > 0 && scalar @list2 > 0) {
				$scalars{1} += scalar @list1;
				$scalars{2} += scalar @list2;
				($smart) ? (&smart_process( \@list1, \@list2, $group, $sub )) : (&process( \@list1, \@list2, $group, $sub ));
			} elsif (scalar @list2 > 0) {
				&report("Skipping group $group $sub: ", scalar @list1, " (list 1) vs ", scalar @list2, " (list 2)\n");	# reporter
				push @incomparable2, "$group\t$sub\t$_\n" foreach @list2;
				$skipped{2} += scalar @list2;
			} else {
				&report("Skipping group $group $sub: ", scalar @list1, " (list 1) vs ", scalar @list2, " (list 2)\n");	# reporter
				push @incomparable1, "$group\t$sub\t$_\n" foreach @list1;
				$skipped{1} += scalar @list1;
			}
			reset(X);
		}
	}
} else {
	print "Matching without subgroups...\n";
	foreach $group (sort keys %pools) {
		my @list1 = sort {$a <=> $b} @{ $pools{$group}{1} };
		my @list2 = sort {$a <=> $b} @{ $pools{$group}{2} };
		if (scalar @list1 > 0 && scalar @list2 > 0) {
			$scalars{1} += scalar @list1;
			$scalars{2} += scalar @list2;
			($smart) ? (&smart_process( \@list1, \@list2, $group )) : (&process( \@list1, \@list2, $group ));
		} elsif (scalar @list2 > 0) {
			&report("Skipping group $group: ", scalar @list1, " (list 1) vs ", scalar @list2, " (list 2)\n");	# reporter
			push @incomparable2, "$group\t$_\n" foreach @list2;
			$skipped{2} += scalar @list2;
		} else {
			&report("Skipping group $group: ", scalar @list1, " (list 1) vs ", scalar @list2, " (list 2)\n");	# reporter
			push @incomparable1, "$group\t$_\n" foreach @list1;
			$skipped{1} += scalar @list1;
		}
		reset(X);
	}
}

open OUT1, "> matched_$basename.txt" or die "Cannot open matched_$basename.txt!";
print OUT1 "Chr\tStart1\tEnd1\tVol1\t\tStart2\tEnd2\tVol2\t\tShift\tExpansion\n";
print OUT1 @output;
close OUT1;

if (@unmatched1 || @unmatched2) {
	open OUT2, "> unmatched_$basename.txt";
	print OUT2 "Unmatched from $lname1:\n";
	print OUT2 @unmatched1;
	print OUT2 "\n\nUnmatched from $lname2:\n";
	print OUT2 @unmatched2;
	close OUT2;
}

if (@incomparable1 || @incomparable2) {
	open OUT3, "> incomparable_$basename.txt";
	print OUT3 "Incomparables from $lname1:\n";
	print OUT3 @incomparable1;
	print OUT3 "\n\nIncomparables from $lname2:\n";
	print OUT3 @incomparable2;
	close OUT3;
}

if (@errors) {
	open OUT4, "> errors_$basename.txt";
	print OUT4 @errors;
	close OUT4;
}

&report("\nThe Numbers:\t\tList 1\tList 2\n");
&report("Filenames:\t\t$lname1\t$lname2\n");
&report("Elements compared:\t$all{1}\t$all{2}\n");
&report("Elements skipped:\t$skipped{1}\t$skipped{2}\n");
&report("Elements matched:\t$cumatch{1}\t$cumatch{2}  (", scalar @output, " matches)\n");
&report("Elements unmatched:\t$unmatched{1}\t$unmatched{2}\n");
&report("\nMatch Types:\n");
&report("$_:\t$mtypes{$_}\n") foreach (sort keys %mtypes);

open LOG, "> match_stats_$basename.txt";
print LOG @log;
close LOG;

print "\n\n";
exit;


########################################     SUBROUTINES     ########################################
sub smart_process {
#####  NEEDS TO BE TESTED AND PROBABLY REWRITTEN SOMEWHAT
	my ($listref1, $listref2, $group, $sub);
	(scalar @_ == 4) ? ( ($listref1, $listref2, $group, $sub) = @_ ) : ( ($listref1, $listref2, $group) = @_ );
	my $psub = " $sub" if (defined $sub);
	&report("Matching group $group$psub: ", scalar @$listref1, " against ", scalar @$listref2, "\n");
	my $prev_i = 'HERE WE GO';
	$X_first_marker{$prev_i} = 0;
	my $max1 = scalar @$listref1 - 1;
	my $max2 = scalar @$listref2 - 1;
	foreach my $line (@$listref2) {
		my ($start, $end, $else) = split /\t/, $line, 3;
		push @X_starts, $start;
		push @X_ends, $end;
	}
	foreach my $i (0..$max1) {
		my ($start1, $end1, $else1) = split /\t/, $$listref1[$i], 3;
		my $marker = $X_first_marker{$prev_i};
		foreach my $j ($marker..$max2) {
			my $start2 = $X_starts[$j];
			my $end2 = $X_ends[$j];
			next if ($end2 < $start1);
			last if ($start2 > $end1);
			$X_first_marker{$i} = $j unless (exists $X_first_marker{$i});
			&matchmaker($start1, $end1, $start2, $end2, $lead, $element, $$listref2[$j]);
		}
		$prev_i = $i;
	}
}

sub process {
	my ($listref1, $listref2, $group, $sub, %matched);
	(scalar @_ == 4) ? ( ($listref1, $listref2, $group, $sub) = @_ ) : ( ($listref1, $listref2, $group) = @_ );
	my $psub = " $sub" if (defined $sub);
	&report("Matching group $group$psub: ", scalar @$listref1, " (list 1) against ", scalar @$listref2, " (list 2)\n");	# reporter
	($matching) ? ($lead = "$group\t$sub") : ($lead = $group);
	foreach my $line (@$listref2) {
		my ($start, $end, $else) = split /\t/, $line, 3;
		push @X_starts, $start;
		push @X_ends, $end;
	}
	$i = 0;
	foreach $element (@$listref1) {
#		($matching) ? (print "$group: $sub: $element\t:\t", scalar @X_starts, "\n") : (print "$group: $element\t:\t", scalar @X_starts, "\n");	# reporter
		my ($start1, $end1, $else1) = split /\t/, $element, 3;
		foreach $j (0..$#X_starts) {
#			print "\t\t$j\n";	# reporter
			my $start2 = $X_starts[$j];
			my $end2 = $X_ends[$j];
			next if ($end2 < $start1);			# not there yet
			last if ($start2 > $end1);			# gone too far
#			print "\t$$listref2[$j]\n";	# reporter
			$matched{1}{$i}++;
			$matched{2}{$j}++;
			&matchmaker($start1, $end1, $start2, $end2, $lead, $element, $$listref2[$j]);
		}
		$i++;
	}
	$cumatch{$_} += scalar (keys %{ $matched{$_} }) foreach (1..2);
	if (scalar @X_temp > 0) {
		foreach (@X_temp) {
#			push @output, $_ unless ($_ =~ /\|/);	# DATASET-SPECIFIC FILTERING...  This one kills intergenics, if > 1 entry, if matching against YGIG dataset.
			push @output, $_;
		}
	}
	foreach $i (0..$#$listref1) {
		unless (exists $matched{1}{$i}) {
			push @unmatched1, "$lead\t$$listref1[$i]\n";
			$unmatched{1}++;
		}
	}
	foreach $j (0..$#$listref2) {
		unless (exists $matched{2}{$j}) {
			push @unmatched2, "$lead\t$$listref2[$j]\n";
			$unmatched{2}++;
		}
	}
}

sub report {
	my $msg = join '', @_;
	push @log, $msg;
	print $msg;
}

sub matchmaker {
	($start1, $end1, $start2, $end2, $lead, $element, $lr2j) = @_;
	if ($start1 == $start2 && $end1 == $end2) {
		push @X_temp, "$lead\t$element\tmatched\t$lr2j\t:\t1.0\t1.0\n";
		$mtypes{'Exact'}++;
	} elsif ($start1 <= $start2 && $end1 >= $end2) {
		$expansion = ($end1 - $start1 + 1) / ($end2 - $start2 + 1);
		push @X_temp, "$lead\t$element\tmatched\t$lr2j\t:\t\t$expansion\n";
		$mtypes{'Eclipse 1'}++;
	} elsif ($start1 >= $start2 && $end1 <= $end2) {
		$expansion = ($end1 - $start1 + 1) / ($end2 - $start2 + 1);
		push @X_temp, "$lead\t$element\tmatched\t$lr2j\t:\t\t$expansion\n";
		$mtypes{'Eclipse 2'}++;
# ?		$matched1{$i}++;
# ?		$matched2{$j}++;
	} elsif ($start1 >= $start2 && $start1 <= $end2) {
		$shift = ($end2 - $start1 + 1) / ($end2 - $start2 + 1);
		push @X_temp, "$lead\t$element\tmatched\t$lr2j\t:\t+$shift\t\n";	# don't forget the "+$shift" printed here!  Indicates downstream shift
		$mtypes{'Overlap 1'}++;					# 1 or 2 indicates which element is 3'
	} elsif ($end1 >= $start2 && $end1 <= $end2) {
		$shift = ($end1 - $start2 + 1) / ($end2 - $start2 + 1);
		push @X_temp, "$lead\t$element\tmatched\t$lr2j\t:\t-$shift\t\n";	# don't forget the "-$shift" printed here!  Indicates upstream shift
		$mtypes{'Overlap 2'}++;					# 1 or 2 indicates which element is 3'
	} else {						# problem
		push @errors, "$lead\t$element\t\t$lr2j\n";
	}
}
