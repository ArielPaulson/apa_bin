#!/usr/bin/env perl
use strict;


## Takes a multisample SPIA output location and makes a navigation HTML
## SPIA dir is output location of apa_tools.R::run.SPIA(), or at least, has same structure

## FUTURE: NDE mouseover popup with DE genes

my ($dir, $org, $build) = @ARGV;

die "$0: input location '$dir' does not exist!\n" unless -d $dir;
die "$0: KEGG organism '$org' not specified!\n" unless $org;
die "$0: KEGG build '$build' not specified!\n" unless $build;
my $keggdir = "/n/data1/KEGG/$org/$build";
die "$0: expected KEGG build '$keggdir' does not exist!\n" unless -d $keggdir;
chomp($dir = `readlink -f $dir`);

## Some HTML shortcuts
my $tr_o  = "  <tr>\n";     # table row open
my $tr_c  = "  </tr>\n";    # table row close
my $td_o  = "    <td>\n";   # table data open
my $td_c  = "    </td>\n";  # table data close
my $tc_o  = "      ";       # table data content "open"  (prefix)
my $tc_c  = "<br>\n";       # table data content "close" (suffix)
my $tr_og = "  <tr class=\"g\">\n";     # table row open, tr grey class
my $td_ob = "    <td class=\"b\">\n";   # table data open, td bold class
my $td_or = "    <td class=\"r\">\n";   # table data open, td right class

## Find files
my %data;
foreach my $contrast (glob "$dir/*") {
    next unless -d $contrast;
    $contrast =~ s/$dir\///;
    next if $contrast eq 'data' || $contrast eq 'pathview_data';  # ignore PathView data dir
    foreach my $file (glob "$dir/$contrast/*") {
        if ($file =~ /SPIAPerturbationPlots/) {
            ($data{$contrast}{SPP} = $file) =~ s/$dir\///;
        } elsif ($file =~ /2way_evidence/) {
            ($data{$contrast}{EVI} = $file) =~ s/$dir\///;
        } elsif ($file =~ /SPIA_results/) {
            ($data{$contrast}{RES} = $file) =~ s/$dir\///;
            ($data{$contrast}{RES2} = $data{$contrast}{RES}) =~ s/\.txt$/.html/;
        } elsif ($file =~ /sig-genes/) {
            (my $lfile = $file) =~ s/$dir\///;
            push @{ $data{$contrast}{SIG} }, $lfile;
        } elsif ($file =~ /all-genes/) {
            (my $lfile = $file) =~ s/$dir\///;
            push @{ $data{$contrast}{ALL} }, $lfile;
        } else {
            print "$0: Don't know what to do with file: '$file'!\n";
        }
    }
}

## Start master HTML document
my $html_init = "<html>
<head>
<style \"text/css\">
  td { padding: 5px 5px 5px 5px; }
  td.b { text-align:center; vertical-align:center; font-weight:bold; }
  td.r { text-align:right; }
  tr.g { background-color:#DDD; }
</style>
<title>SPIA Output Navigation</title>
</head>

<body>
<h2>SPIA Output Navigation Page</h2>
This HTML page should live here: <a href=\"http://bioinfo$dir\">$dir</a><br>
KEGG build $org/$build used for analysis: <a href=\"http://bioinfo$keggdir\">$keggdir</a><br>
<br>
<table width=\"100%\" border=\"1\">
"   . $tr_o 
    . &make_td('CONTRAST','B',0) 
    . &make_td('QC','B',0) 
    . &make_td('PATHWAY PLOTS<br>DE GENES','B',0) 
    . &make_td('PATHWAY PLOTS<br>ALL GENES','B',0) 
    . $tr_c;   # add table header row

open HTML, '>', "$dir/SPIA_index.html";
print HTML $html_init;


## Name ID pSize NDE PctDE pNDE tA pPERT pG pGFdr pGFWER Significant Status KEGGLink
my @sci = (5,7,8,9,10);
my @pct = (4);
my @flt = (6);

## Process contrasts into HTML rows
foreach my $contrast (sort keys %data) {
    
    ## First, read in SPIA table / get data, and convert results table to HTML
    my %pathlabels;

    if (exists $data{$contrast}{RES}) {
        
        ## Initialize results-table HTML
        open OUT, '>', $data{$contrast}{RES2};
        my $line = "<html>
<head>
<style \"text/css\">
  td { padding: 5px 5px 5px 5px; }
  td.b { text-align:center; vertical-align:center; font-weight:bold; }
  td.r { text-align:right; }
  tr.g { background-color:#DDD; }
</style>
<title>SPIA Results Table</title>
</head>

<body>
<table border=\"1\">
";
        print OUT $line;
        
        ## Read results / fill in results-table HTML
        open IN, '<', $data{$contrast}{RES};
        my $haspct;
        while (my $line = <IN>) {
            chomp($line);
            my @fields = split /\t/, $line;
            
            if ($. == 1) {
                
                ## header row
                $haspct = $fields[4] eq 'PctDE';
                unless ($haspct) {
                    @fields = (@fields[0..3], 'PctDE', @fields[4..$#fields]);
                }
                print OUT $tr_o;
                print OUT &make_td($fields[$_],'B',0) foreach (0..$#fields);
                print OUT $tr_c;
                
            } else {
                
                ## data rows
                $pathlabels{$fields[1]} = "$fields[1] $fields[0]";
                unless ($haspct) {
                    @fields = (@fields[0..3], 100*$fields[3]/$fields[2], @fields[4..$#fields]);
                }
                my @fields2 = @fields;
                $fields2[$_] = sprintf("%0.2e",$fields[$_]) foreach @sci;
                $fields2[$_] = sprintf("%0.2f",$fields[$_]) foreach @pct;
                $fields2[$_] = sprintf("%0.6f",$fields[$_]) foreach @flt;
                $fields2[-1] = "<a href=\"$fields[-1]\">KEGG Pathway Map</a>";
                $fields2[1] = "<a href=\"http://www.genome.jp/dbget-bin/www_bget?pathway+$fields[1]\">$fields[1]</a>";
                print OUT $. % 2 == 0 ? $tr_og : $tr_o;
                print OUT &make_td($fields2[$_],'',0) foreach (0,1);
                print OUT &make_td($fields2[$_],'R',0) foreach (2..10);
                print OUT &make_td($fields2[$_],'',0) foreach (11..$#fields2);
                print OUT $tr_c;
            }
        }
        close IN;
        
        ## Finalize results-table HTML
        print OUT "</table>\n<br>\n</body>\n</html>\n";
        close OUT;
        
        
        ## Second, add master-HTML row for contrast
        
        my $qc_td = $td_o
            . &make_href($data{$contrast}{SPP},"Perturbation Plots (PDF)",3)  # sub out 'dir' to get local-link version
            . &make_href($data{$contrast}{EVI},"2-Way Evidence (PNG)",3)
            . &make_href($data{$contrast}{RES2},"Results Table (HTML)",3)
            . &make_href($data{$contrast}{RES},"Results Table (TXT)",1)
            . $td_c;
        
        my (@sig_td, @all_td);
        if (exists $data{$contrast}{SIG}) {
            foreach my $sig (@{ $data{$contrast}{SIG} }) {
                my ($pathid) = ($sig =~ /($org\d{5})\D/);
                push @sig_td, &make_href($sig,$pathlabels{$pathid},0);
            }
        }
        if (exists $data{$contrast}{ALL}) {
            foreach my $all (@{ $data{$contrast}{ALL} }) {
                my ($pathid) = ($all =~ /($org\d{5})\D/);
                push @all_td, &make_href($all,$pathlabels{$pathid},0);
            }
        }
        my $sig_td = @sig_td ? join("<br>\n", @sig_td) : "${tc_o}N/A";
        my $all_td = @all_td ? join("<br>\n", @all_td) : "${tc_o}N/A";
        
        my $html_row = $tr_o
            . &make_td($contrast,'B',0)
            . $qc_td
            . &make_td($sig_td,'',1)
            . &make_td($all_td,'',1)
            . $tr_c;
        print HTML $html_row;
    }
    
}

## Finalize master HTML
print HTML "</table>\n<br>\n</body>\n</html>\n";
close HTML;





exit;



sub make_td {
    my ($content, $class, $no_tco) = @_;
    my $tdo = ($class eq 'B' ? $td_ob : $class eq 'R' ? $td_or : $td_o);
    return $tdo . ($no_tco ? '' : $tc_o) . "$content\n$td_c";    # newline is for source (\n), not page rendering (<br>)
}

sub make_href {
    my ($link, $label, $newline) = @_;
    ## $newline: 0 for none, 1 for source only, 2 for page render only, 3 for source AND page rendering
    return "$tc_o<a href=\"$link\">$label</a>" . ($newline==0 ? "" : $newline==1 ? "\n" : $newline==2 ? "<br>" : $tc_c);
}


