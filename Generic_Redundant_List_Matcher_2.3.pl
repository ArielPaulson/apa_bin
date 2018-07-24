#!/usr/bin/perl

### This script is a generalized version of the Solexa interpreter algorithm, which can take two ranked lists of coordinates and return matches, with or without overlaps.
### Element matching is accelerated by position tracking within the arrays, so that distant upstream/downstream elements are not needlessly searched.
### Both lists can have redundant and/or overlapping entries without creating problems.
### Column structure of the lists is based on track file data, i.e. {CHR#\tSTART\tEND\tDATA} or {CHR#\tSTART\tEND\tSTRAND\tDATA}.
### Generalized, these become {GROUP\tSTART\tEND\tDATA} and {GROUP\tSTART\tEND\tSUBGROUP\tDATA}, respectively.
### Use of subgroup matching is optional and can be used if the [strand / subgroup] column is present.
### Overlaps may be reported as such (option 'o'), and return # overlapping bp, or they may be included with matches (option 'm'), or simply ignored (option 'i').
### Nearest Neighbors for unmatched elements can be added to the unmatched list by using ReportNeighbors (flag 6) = "y".
### 
### COMMAND-LINE SUMMARY: > perl Generic_Redundant_List_Matcher.pl [List1 (first file)] [List2 (second file)] [y/n (subgroup matching)] [y/n (restrict to exacts)] [y/n (Neighbor reporting)]
### E.g. 'perl Generic_Redundant_List_Matcher.pl peaks.txt yeast_genes.txt n y o y'
### All five options must be used (and used correctly) for command-line execution, else script will die.  If no options given, script will prompt for input.

########################################     SET PARAMETERS     ########################################

if (scalar @ARGV > 0) {
	($#ARGV == 4) ? (($list1, $list2, $match1, $exact1, $neighbor1) = @ARGV) : (die "Inputs not equal to 5!  Cannot proceed.  5 elements required:\nList1 List2 SubgroupMatch?(y/n) ReportExactsOnly?(y/n) ReportNeighbors?(y/n)\n");
	if ($match1 =~ /^y$/i) {
		$matching = 1;
	} elsif ($match1 =~ /^n$/i) {
	} else {
		die "Subgroup matching choice must be [yY] or [nN] only!\n[yY]: match subgroups within groups\n[nN]: ignore subgroups during matching\n";
	}
	if ($exact1 =~ /^y$/i) {
		$exactonly = 1;
	} elsif ($exact1 =~ /^n$/i) {
		if ($neighbor1 =~ /^y$/i) {
			$neighbors = 1;
		} elsif ($neighbor1 =~ /^n$/i) {
		} else {
			die "Nearest-neighbor reporting choice must be [yY] or [nN] only!\n";
		}
	} else {
		die "Report only exact matches choice must be [yY] or [nN] only!\n";
	}
} else {
	print "\nLists must be in 'GROUP\\tSTART\\tEND\\tSUBGROUP\\tELSE\\n' format and MUST HAVE HEADERS.\n\nSubgroup matching is optional; field may be omitted ('.' = ignore subgrouping)\n";
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
	print "\nReport only exact matches? (y/n)\n";
	{
		chomp($exacts1 = <STDIN>);
		if ($exacts1 =~ /^y$/i) {
			$exactonly = 1;
		} elsif ($exacts1 =~ /^n$/i) {

			print "\nReport nearest neighbors for unmatched entries? (y/n)\n";
			{
				chomp($neighbor1 = <STDIN>);
				if ($neighbor1 =~ /^y$/i) {
					$neighbors = 1;
				} elsif ($neighbor1 =~ /^n$/i) {
				} else {
					print "Enter [yY] or [nN] only!\n";
					redo;
				}
			}

		} else {
			print "Enter [yY] or [nN] only!\n";
			redo;
		}
	}


}

########################################     PREPARE FOR PROCESSING     ########################################

my (%pools, @output, %mtypes, $lname1, $lname2, %headers, %hsplit);
$mtypes{$_} = 0 foreach ('Exact', 'Eclipse 1', 'Eclipse 2', 'Overlap 1', 'Overlap 2');	# ensure printable values
($list1 =~ /\/?([^\/]+)$/) ? ($lname1 = $1) : ($lname1 = $list1);
($list2 =~ /\/?([^\/]+)$/) ? ($lname2 = $1) : ($lname2 = $list2);
my $basename = join '_', ($lname1, $lname2, $match1, $exact1);

foreach (1..2) {	# ensure printable values
	$all{$_} = 0;
	$cumatch{$_} = 0;
	$skipped{$_} = 0;
	$unmatched{$_} = 0;
}

open IN1, $list1 || die "Cannot open $list1!\n";
print "\nProcessing $lname1...\n";
while (<IN1>) {
     next if ($_ =~ /^#/ || $_ =~ /^track/);
	$all{1}++;
	my ($group, $start, $end, $else);
	$line = $_;
	$line =~ s/[\n\r]//g;
	if ($all{1} == 1) {	# header
		$headers{1} = $line;
	} elsif ($matching) { 
		my ($group, $start, $end, $sub, $else) = split /\t/, $line, 5 || die "Subgroup matching requires 5 columns of data minimum!\n";
		push @{ $pools{$group}{$sub}{1} }, "$start\t$end\t$else";
		$hsplit{1} = 5;
	} else {
		my ($group, $start, $end, $else) = split /\t/, $line, 4;
		push @{ $pools{$group}{1} }, "$start\t$end\t$else";
		$hsplit{1} = 4;
	}
}
close IN1;

open IN2, $list2 || die "Cannot open $list2!\n";
print "Processing $lname2...\n";
while (<IN2>) {
     next if ($_ =~ /^#/ || $_ =~ /^track/);
	$all{2}++;
	my ($group, $start, $end, $else);
	$line = $_;
	$line =~ s/[\n\r]//g;
	if ($all{2} == 1) {	# header
		$headers{2} = $line;
	} elsif ($matching) {
		my ($group, $start, $end, $sub, $else) = split /\t/, $line, 5 || die "Subgroup matching requires 5 columns of data minimum!\n";
		push @{ $pools{$group}{$sub}{2} }, "$start\t$end\t$else";
		$hsplit{2} = 5;
	} else {
		my ($group, $start, $end, $else) = split /\t/, $line, 4;
		push @{ $pools{$group}{2} }, "$start\t$end\t$else";
		$hsplit{2} = 4;
	}
}
close IN2;

my @header1 = split /\t/, $headers{1}, $hsplit{1};
my @header2 = split /\t/, $headers{2}, $hsplit{2};
my $intro;
($hsplit{1} == 5) ? ($intro = "GROUP\tSUBGROUP") : ($intro = "GROUP");
my $allheader = join "\t", ($intro, "START1\tEND1", @header1[$hsplit{1}-1], "MATCHED\tSTART2\tEND2", @header2[$hsplit{2}-1], ":\tOVERLAP\tMATCHTYPE");
(my $allheaderE = $allheader) =~ s/\tMATCHED\tSTART2/\tERROR\tSTART2/;

if ($matching) {
	print "Matching with subgroups...\n";
	foreach $group (sort keys %pools) {
		foreach $sub (sort keys %{ $pools{$group} }) {
			my @list1 = sort {$a <=> $b} @{ $pools{$group}{$sub}{1} };
			my @list2 = sort {$a <=> $b} @{ $pools{$group}{$sub}{2} };
			if (scalar @list1 > 0 && scalar @list2 > 0) {
				$scalars{1} += scalar @list1;
				$scalars{2} += scalar @list2;
				&process( \@list1, \@list2, $group, $sub );
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
			&process( \@list1, \@list2, $group );
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

open OUT1, "> matched_$basename.txt" or die "Cannot open matched_$basename.txt";
print OUT1 "$allheader\n",@output;
close OUT1;

if (@unmatched1 || @unmatched2) {
	open OUT2, "> unmatched_$basename.txt";
	print OUT2 "Unmatched from $lname1:\n";
	print OUT2 "$headers{1}\n";
	print OUT2 @unmatched1;
	print OUT2 "\n\nUnmatched from $lname2:\n";
	print OUT2 "$headers{2}\n";
	print OUT2 @unmatched2;
	close OUT2;
}

if (@incomparable1 || @incomparable2) {
	open OUT3, "> incomparable_$basename.txt";
	print OUT3 "Incomparables from $lname1:\n";
	print OUT3 "$headers{1}\n";
	print OUT3 @incomparable1;
	print OUT3 "\n\nIncomparables from $lname2:\n";
	print OUT3 "$headers{2}\n";
	print OUT3 @incomparable2;
	close OUT3;
}

if (@errors) {
	open OUT4, "> errors_$basename.txt";
	print OUT4 "$allheaderE\n",@errors;
	close OUT4;
}

if ($exactonly && @mismatched) {
	open OUT5, "> mismatched_$basename.txt";
	print OUT5 "$allheader\n",@mismatched;
	close OUT5;
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

sub process {
	my ($listref1, $listref2, $group, $sub, %matched, %nn);
	(scalar @_ == 4) ? ( ($listref1, $listref2, $group, $sub) = @_ ) : ( ($listref1, $listref2, $group) = @_ );
	my $psub = " $sub" if (defined $sub);
	&report("Matching group $group$psub: ", scalar @$listref1, " (list 1) against ", scalar @$listref2, " (list 2)\n");	# reporter
	($matching) ? ($lead = "$group\t$sub") : ($lead = $group);
	foreach my $line (@$listref2) {
		my ($start, $end, $else) = split /\t/, $line, 3;
		push @X_starts, $start;
		push @X_ends, $end;
	}
	my %nnspacer;
	$nnspacer{1} = $nnspacer{2} = "\t";	# for "$distN" entry
	$nnspacer{1} .= "\t" while $$listref1[0] =~ /\t/g;	# for remaining entries
	$nnspacer{2} .= "\t" while $$listref2[0] =~ /\t/g;

	## Match List 1 elements to List 2 elements
	$i = $first = 0;
	foreach $element (@$listref1) {
#		($matching) ? (print "$group: $sub: $element\t:\t", scalar @X_starts, "\n") : (print "$group: $element\t:\t", scalar @X_starts, "\n");	# reporter
		my ($start1, $end1, $else1) = split /\t/, $element, 3;
		my $hit;
		foreach $j ($first..$#X_starts) {
#			print "\t\t$j\n";			# reporter
			my $start2 = $X_starts[$j];
			my $end2 = $X_ends[$j];
			if ($end2 < $start1) {			# not there yet
				unless ($matched{2}{$j}) {	# unmatched $j; get neighbors
					my $dist3 = $end1 - $start2;	# 3' neighbor guaranteed
					$nn{ $$listref2[$j] }{3} = "$dist3\t$$listref1[$i]";
					if ($i == 0) {		# no 5' neighbor
						$nn{ $$listref2[$j] }{5} = $nnspacer{1};
					} else {
						my $i2 = $i;
						{
							$i2--;	# move backwards until nearest 5' neighbor is found, or list beginning is reached
							my ($start1B, $end1B, $else1B) = split /\t/, $$listref1[$i2], 3;
							if ($end1B < $start2) {
								my $dist5 = $start2 - $end1B;
								$nn{ $$listref2[$j] }{5} = "$dist5\t$$listref1[$i2]";
							} else {
								($i2 > 0) ? (redo) : ($nn{ $$listref2[$j] }{5} = $nnspacer{1});	# continue 5', or quit
							}
						}
					}
				}
				if ($j == $#X_starts && !$hit) {	# if we are here, then list 2 ran out before list1 (5' neighbor only)
					my $dist5 = $start1 - $end2;
					$nn{ $$listref1[$i] }{5} = "$dist5\t$$listref2[$j]";
					$nn{ $$listref1[$i] }{3} = $nnspacer{2};
				}
				next;
			} elsif ($start2 > $end1) {		# gone too far
				unless ($hit) {			# unmatched $i; get neighbors
					my $dist3 = $start2 - $end1;	# 3' neighbor guaranteed
					$nn{ $$listref1[$i] }{3} = "$dist3\t$$listref2[$j]";
					if ($j == 0) {		# no 5' neighbor
						$nn{ $$listref1[$i] }{5} = $nnspacer{2};
					} else {
						my $j2 = $j;
						{
							$j2--;	# move backwards until nearest 5' neighbor is found, or list beginning is reached
							my ($start2B, $end2B, $else2B) = split /\t/, $$listref2[$j2], 3;
							if ($end2B < $start1) {
								my $dist5 = $start1 - $end2B;
								$nn{ $$listref1[$i] }{5} = "$dist5\t$$listref2[$j2]";
							} else {
								($j2 > 0) ? (redo) : ($nn{ $$listref1[$i] }{5} = $nnspacer{2});	# continue 5', or quit
							}
						}
					}
				}
				last;
			} else {				# match of some kind
				$first = $j unless $hit;	# next time, start at the first match position
				$hit++;
#				print "\t$$listref2[$j]\n";	# reporter
				$matched{1}{$i}++;
				$matched{2}{$j}++;
				&matchmaker($start1, $end1, $start2, $end2, $lead, $element, $$listref2[$j]);
			}
		}
		$i++;
	}

	## Neighborize all trailing unmatched List 2 elements
	if ($neighbors) {
		$i--;	# was 1 too far
		foreach $j ($first..$#X_starts) {
			unless ($matched{2}{$j}) {
				my ($start1, $end1, $else1) = split /\t/, $$listref1[$i], 3;
				my ($start2, $end2, $else2) = split /\t/, $$listref2[$j], 3;
				my $dist5 = $start2 - $end1;
				$nn{ $$listref2[$j] }{5} = "$dist5\t$$listref1[$i]";	# if unmatched $j, get 5' neighbors (since there are no 3's)
			}
		}
	}

	## Tally stuff and store outputs
	$cumatch{$_} += scalar (keys %{ $matched{$_} }) foreach (1..2);
	if (scalar @X_temp > 0) {
		foreach (@X_temp) {
#			push @output, "$_\n" unless ($_ =~ /\|/);	# DATASET-SPECIFIC FILTERING...  This one kills intergenics, if > 1 entry, if matching against YGIG dataset.
			if ($exactonly) {
				($_ =~ /Exact$/) ? (push @output, "$_\n") : (push @mismatched, "$_\n");
			} else {
				push @output, "$_\n";
			}
		}
	}
	foreach $i (0..$#$listref1) {
		unless (exists $matched{1}{$i}) {
			my $string = "$lead\t$$listref1[$i]";
			($neighbors) ? (push @unmatched1, "$string\t\|\t$nn{ $$listref1[$i] }{5}\t\|\t$nn{ $$listref1[$i] }{3}\n") : (push @unmatched1, "$string\n");
			$unmatched{1}++;
		}
	}
	foreach $j (0..$#$listref2) {
		unless (exists $matched{2}{$j}) {
			my $string = "$lead\t$$listref2[$j]";
			($neighbors) ? (push @unmatched2, "$string\t\|\t$nn{ $$listref2[$j] }{5}\t\|\t$nn{ $$listref2[$j] }{3}\n") : (push @unmatched2, "$string\n");
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
	my ($dist, $flag);
	if ($start1 == $start2 && $end1 == $end2) {
		($dist, $flag) = ($end1 - $start1 + 1, 'Exact');
	} elsif ($start1 >= $start2 && $end1 <= $end2) {
		($dist, $flag) = ($end1 - $start1 + 1, 'Eclipse 1');
	} elsif ($start1 <= $start2 && $end1 >= $end2) {
		($dist, $flag) = ($end2 - $start2 + 1, 'Eclipse 2');
	} elsif ($start1 >= $start2 && $start1 <= $end2) {
		($dist, $flag) = ($end2 - $start1 + 1, 'Overlap 1');	# 1 or 2 indicates the 3' element
	} elsif ($end1 >= $start2 && $end1 <= $end2) {
		($dist, $flag) = ($end1 - $start2 + 1, 'Overlap 2');	# 1 or 2 indicates the 3' element
	} else {					# problem
		push @errors, "$lead\t$element\tERROR\t$lr2j\n";
	}
	push @X_temp, "$lead\t$element\tmatched\t$lr2j\t:\t$dist\t$flag";
	$mtypes{$flag}++;
}