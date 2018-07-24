#!/usr/bin/env perl
# OLD: #!/n/local/stage/perlbrew/perlbrew-0.43/perls/perl-5.16.0/bin/perl
# OLDER: #!/n/site/inst/Linux-x86_64/sys/bin/perl
use File::Path qw/ make_path remove_tree /;
use Storable qw/ nstore retrieve dclone /;
use Scalar::Util qw(looks_like_number);
use Capture::Tiny ':all';
use Data::Dumper;
use File::Spec;
use strict;
use Roman;
chomp(my $tilde = `readlink -f ~`);





## Author: Ariel Paulson (apa@stowers.org)
## A Perl module of common custom functions





### SCRIPT SETUP/QC AND I/O FUNCTIONS:
## validate          : validates files; creates and validates temp directories, with clobber/no-clobber handling
## create_writedir   :  -- temporarily deprecated, I think --
## test_writability  : ensure location or prefix is writable
## is_compressed     : returns array: ( 1|0 if properly compressed , compression error message , linux read command, linux write command )
## open2             : open() wrapper with handling for compressed files, stdio
## generic_io        : opens basic input and output handles, tests validity of args, given $input, $output, and $stdio
## backticks_bash    : backticks-substitute wired to use bash; allows commands with single quotes, unlike "bash -c '$cmd'"
## execute           : system() wrapper that sends command to screen, optional logfile, and runs command
## logreport         : send a message to screen & logfile simultaneously
## print_table       : print array-of-arrays object to a filehandle using a delimiter, or sprintf-ed
## nfsAnalyze        : temp directory not getting removed?  Run this on it to see what files are locking and why (prints to STDERR)
## validate_genoanno : testing validity of $geno/$anno in /n/data1/genomes/indexes, also STAR indexes
## validate_bam      : checks if a bam file is corrupted

### GENERIC FUNCTIONS:
## timestamp         : three different types of platform-independent timestamps
## runtime           : takes 2 timestamps and returns elapsed time, in various scales
## monthconv         : convert month names to numbers and vice versa
## benchmark         : platform-independent system resource stamps
## numjust           : pad a digit with leading zeroes
## unique            : return unique values from an array, in order of appearance
## transpose         : transposes a matrix, stored in an array of arrays
## pairs2subgraphs   : given a 2-level hash describing a connectivity matrix (e.g. $hash{$key1}{$key2}, preferably reciprocal so that $hash{$key2}{$key1} also always exists), return all subgraphs.
## validate_memory   : compares requested RAM with server RAM; if too much either downsizes or dies
## canonicalize      : "readlink -f" without the link-following
## extract_names_from_paths  : port of the apa_tools.R function.  Finds the smallest defining substrings from a list of paths.

### MATH FUNCTIONS:
## mean              : mean of an array
## median            : median of an array
## stdev             : standard deviation of an array
## strip_NA          : remove NA entries from array

### BIOLOGICAL FUNCTIONS:
## revcomp           : reverse complement for DNA/RNA, with full degeneracy/masking support
## revcomp_regex     : reverse complement a DNA/RNA regular expression, with full degeneracy/masking support
## word_exp          : get the expected frequency of a DNA/RNA/AA word, given background character frequencies
## coding_potential  : calculate coding potential for a DNA sequence, using one of several methods    ########## NOT WORKING PROPERLY ##########
## translate         : translate a DNA sequence
## blockify          : convert a sequence to a fasta block
## chrsort           : sort a list of chromosomes the way I prefer.  Special handling for roman chrnums (yeast) and drosophila chrs
## get_memedata      : get a standard hash of data from a MEME run (DOES NOT WORK FOR 4.8.1)
## UTR_weld          : take a hash of CDS entries + hash of UTR entries and return a hash of exons

### UNDER CONSTRUCTION: may have some functionality at this time, or may not
## readFile          : read a tabular text file into an array-of-arrays or hash-of-arrays
## writeFile         : write an array, array-of-arrays, or hash-of-arrays back to a tabular text file
## readFasta         : read a fasta file into a hash
## writeFasta        : write a fasta file from a hash
## get_fimodata      : get a standard hash of data from a FIMO run (DOES NOT WORK FOR 4.8.1)
## get_tomdata       : get a standard hash of data from a Tomtom run (DOES NOT WORK FOR 4.8.1)
## length_regex      : calculate length of a regular expression, also report if length is fixed or minimum
## vcf_orderheader   : ensures VCF header is "properly" sorted; optional genome arg will ensure contig lines are correct and ordered.
## read_org_metadata : reads organism metadata file into a hash (for genome build prep)
## read_VARIABLES    : reads genome/transcriptome prep _VARIABLES/_TVARIABLES file into a a hash





## A few key locations for GenomeIndexes stuff
my $indexes = '/n/data1/genomes/indexes';
my $apabin = $ENV{SIMR_APA_BIN};
my $buildroot = $ENV{SIMR_BUILDROOT};





sub validate {
    
    ## validates files
    ## creates and validates temp directories, with clobber/no-clobber handling.
    ## returns rooted path to input location
    ## INPUT is an array with 4 elements: 
    ## 1. an object type              ('file'|'dir')
    ## 2. an object path              (e.g. "path/to/file")
    ## 3. an optional object label    (e.g. "input file")
    ## 4. an optional binary argument (1|0) which is 'required' for files, and 'clobber' for dirs.
    
    my ($type, $path, $label, $arg) = @_;
    my ($fail, $error, $fullpath);
    my $blurb = $label ? "$label '$path'" : "'$path'";
    
    if ($type eq 'file') {
        
        my $need = $arg;
        
        if ($need) {
            $fail .= "$0: $blurb does not exist!\n" unless -e $path;  # required file
        } else { 
            $fail .= "$0: $blurb does not exist!\n" if defined $path && ! -e $path;  # not required, but specified and non-existent
        }
        
    } elsif ($type eq 'dir') {
        
        my $clobber = $arg;
        
        if ($path && $path ne '.') {  ##### NEVER CLOBBER '.' !!!!!
            if ($clobber) {
                remove_tree($path, {error => \$error});
                sleep 1;
                $fail .= "$0: Cannot remove existing $blurb: $error\n" if -d $path;
            } else {
                $fail .= "$0: $blurb already exists!\n" if -d $path;
            }
            make_path($path, {error => \$error});
            sleep 1;
            if (!$error && ! -d $path) {
                sleep 5;
                make_path($path, {error => \$error});  # try once more, if unsuccessful
                sleep 1;
            }
            $fail .= "$0: Failed to create $blurb: $error!\n" unless -d $path;
        } else {
            $path = '.';
        }
        
    }
    
    if ($fail) {
        die $fail;
    } else {
        # ROOT THE PATH BEFORE RETURNING
        #$fullpath = File::Spec->rel2abs($path) if $path;  # this module is trash; cannot return true full path, only "relative" full path
        chomp($fullpath = `readlink -f $path`);            # 'readlink' actually works (unless on Mac, but all my code is intended to be run serverside, i.e. CentOS)
        return $fullpath;
    }
}





sub validate_build {
    
    ## Given $geno and $anno (and perhaps a sub-build type), ensures the dataset exists in /n/data1/genomes/indexes
    
    my ($geno, $anno, $sub) = @_;

    my ($prefix, $valid);
    my $genodir = "$indexes/$geno";
    if ($sub) {
        $prefix = $anno ? "$genodir/$anno/$sub/$geno.$anno.$sub" : "$genodir/$sub/$geno.$sub";
        ## sub-builds always have fastas
        ## annotation sub-builds: 'ribosomes', 'funcRNA'
        ## genome sub-builds: 'repeats'
        chomp($valid = `ls $genodir/$anno/$sub/$geno.$anno.$sub.fa 2>/dev/null | wc -l`);
    } else {
        $prefix = "$genodir/$anno/$geno.$anno";
        ## main builds: check for genedata.txt
        chomp($valid = `ls $genodir/$anno/$geno.$anno.genedata.txt 2>/dev/null | wc -l`);
    }
    if ($valid) {
        return $prefix;
    } elsif ($sub && $anno) {
        die "$0: '$sub' sub-build for genome '$geno' annotation '$anno' was not found!\n";
    } elsif ($sub) {
        die "$0: '$sub' sub-build for genome '$geno' was not found!\n";
    } else {
        die "$0: annotation build '$anno' for genome '$geno' was not found!\n";
    }
}





sub create_writedir {
    
    ## A special-case wrapper to &validate, when script has MEME-style output args:
    ##   -o <outdir_don't_clobber>  OR  -oc <outdir_clobber>
    ## Resolves and creates an output directory, or dies
    ## EITHER $outdir_nocl (-o) OR $outdir_cl (-oc) are the output paths, but not both:
    ##   IF BOTH $outdir_nocl AND $outdir_cl: dies.
    ##   IF $outdir_nocl AND NOT $outdir_cl: any existing location ($outdir_nocl) WILL NOT be clobbered.
    ##   IF $outdir_cl AND NOT $outdir_nocl: any existing location ($outdir_cl) WILL be clobbered.
    ## $label is for a friendly die message, if dying is warranted
    ## $default is an OPTIONAL default value for $outdir
    ## returns fully rooted path to output location.
    ## example: $outdir = &create_writedir($o, $oc, 'output location', '.');
    
    my ($outdir_nocl, $outdir_cl, $label, $default) = @_;
    my ($outdir, $clobber);
    
    if ($outdir_nocl && $outdir_cl) {
        die "$0: Cannot specify both -o and -oc!\n";
    } elsif ($outdir_cl) {
        $outdir = $outdir_cl;
        $clobber = 1;
    } elsif ($outdir_nocl) {
        $outdir = $outdir_nocl;
        $clobber = 0;
    } elsif ($default) {
        $outdir = $default;
        $clobber = 0;  # NEVER CLOBBER A DEFAULT LOCATION: MAY BE '.'
    } else {
        die "$0: Failed to specify output directory, and no default!\n";
    }
    
    $label = 'output location' unless $label;
    my $outdir2 = &validate('dir', $outdir, $label, $clobber);
    return $outdir2;
    
}


sub test_writability {
    
    ## Input a path or pathed filename prefix
    ## Can either return 0|1 for writability, or just die if unwritable
    
    my ($location, $die) = @_;
    my $testfile = -d $location ? "$location/test" : "$location.test";
    my $testpath = $location;
    $testpath =~ s/[^\/]+$// unless -d $testpath;
    
    system "touch $testfile";
    sleep 1;
    my $pass = -e $testfile ? 1 : 0;
    system "rm -f $testfile";
    
    die "$0: Output location '$testpath' is not writable!\n" if !$pass && $die;
    return $pass;
}


sub is_compressed {
    
    ## INPUT: 1. filename, 2. OPTIONAL 1|0 indicating to use bgzip for write command
    ## OUTPUT: array: ( is-compressed flag , error message , read command, write command )
    ## recognizes gzip/bgzip/bzip2/xz/zip compression, also plain text and tar
    ## NOTE: read and write commands ALREADY CONTAIN PIPES, so just use like:
    ##       open(IN, $read_command)
    ##       open(OUT, $write_command)
    
    my ($file, $bgzip_flag) = @_;
    my $is_comp;    # 1|0 is file compressed (or will it be compressed)
    my $comp_err;   # error message, if any
    my @iocmd;      # (read command, write command)
    my $ext_type;   # compression type according to file extension
    my $file_type;  # actual compression type (according to 'file' command)
    (my $fext = "\L$file") =~ s/^.*\.//;    # file extension, ensure lowercase
    (my $unfile = $file) =~ s/\.$fext$//i;  # $file without the compression suffix (unused if no compression suffix)
    
    ## Extension -> IO read type lookup
    my @ordtypes = qw/ text bgzip gzip bzip2 zip xz /;  # for test-running purposes: ensure bgzip is tested BEFORE gzip
    my %iotypes = ('xz'=>'xz', 'gz'=>'gzip', 'tgz'=>'gzip', 'bgz'=>'bgzip', 'bz2'=>'bzip2', 'zip'=>'zip');  # '.gz' could also be bgzip, but this is corrected for elsewhere
    $ext_type = $iotypes{$fext};
    
    ## Data for recognized compression types.
    ## Also includes the uncompressed type, WHICH ASSUMES FORMAT IS PLAIN TEXT.
    ## Elements = 1. expected 'file' output regexp, 2. integrity test command, 3. read command, 4. write command.
    my $testx = "$file 2>&1 | grep -P \"\S\"";
    my %typedata = (
        'text'  => { 'REGEX'=>' text$',                         'TEST'=>"",                               'READ'=>"cat        $file |",  'WRITE'=>"| cat      > $file" },
        'xz'    => { 'REGEX'=>'^xz compressed',                 'TEST'=>"bash -c 'xz    -t $testx'",  'READ'=>"unxz -q -c $file |",  'WRITE'=>"| xz -q -f > $file" },
        'gzip'  => { 'REGEX'=>'^gzip compressed',               'TEST'=>"bash -c 'gzip  -t $testx'",  'READ'=>"gunzip  -c $file |",  'WRITE'=>"| gzip  -f > $file" },
        'bgzip' => { 'REGEX'=>'^gzip compressed, extra field',  'TEST'=>"bash -c 'gzip  -t $testx'",  'READ'=>"gunzip  -c $file |",  'WRITE'=>"| bgzip -f > $file" },  # no 'bgunzip', neither is there 'bgzip -t'
        'bzip2' => { 'REGEX'=>'^bzip2 compressed',              'TEST'=>"bash -c 'bzip2 -t $testx'",  'READ'=>"bunzip2 -c $file |",  'WRITE'=>"| bzip2 -f > $file" },
        'zip'   => { 'REGEX'=>'^Zip archive',                   'TEST'=>"bash -c 'zip   -T $testx'",  'READ'=>"unzip   -p $file |",  'WRITE'=>"| cat > $unfile && rm -f $file && zip $file $unfile && rm -f $unfile" }
    );
    ## 'zip' write command is crazy because zip is designed to archive existing files, NOT write from a pipe.
    ##   It CAN write from a pipe, e.g. 'echo "blah" | zip test.zip -', but if you unzip 'test.zip', the filename is '-' not 'test'.
    ##   Thus, command is: 1. write to non-zip temp, 2. destroy existing zip (otherwise it appends), 3. zip temp, 4. destroy temp.
    ##   Yes, it actually works.
    
    ## As written above, commands are heavily justified to make them clearer.
    ## Here, collapse all that whitespace before actually running any of them!  
    ## Not that leaving it in would break anything, but if errors are thrown, it could make reading the message more confusing...
    foreach my $type (keys %typedata) {
        $typedata{$type}{$_} =~ s/ +/ /g foreach keys %{ $typedata{$type} };
    }
    
    ## Determine $file_type
    if (-e $file) {
        
        ## See what 'file' says about the format of $file
        chomp(my $filemsg = `file $file | cut -f2 -d: | sed 's/^ //'`);
        
        foreach my $type (@ordtypes) {    
            ## MUST be ordered!  Test bgzip BEFORE gzip!!
            if ($filemsg =~ /$typedata{$type}{REGEX}/) {
                $file_type = $type;
                @iocmd = ($typedata{$type}{READ}, $typedata{$type}{WRITE});
                last;
            }
        }
        
        ## Compare $file_type and $ext_type
        if ($file_type eq 'text') {
            
            ## some kind of text file (UTF-8, ASCII, etc)
            ## no issues, unless it looks like it SHOULD be compressed...
            $comp_err = "WARNING: filename looks compressed, but data is not!" if $ext_type;
            
        } elsif ($file_type) {
            
            ## then data is compressed
            $is_comp = 1;
            ## test archive integrity
            chomp($comp_err = `$typedata{$file_type}{TEST}`);
            $comp_err = undef if $file_type eq 'zip' && $comp_err eq "test of $file OK";
            
            if ($comp_err) {
                
                ## corrupted archive
                @iocmd = ();  # do not pretend we can read corrupted files
                
            } else {
                
                ## intact archive
                if ($ext_type eq 'gzip' && $file_type eq 'bgzip') {
                    ## OK; this is an unavoidable issue caused by bgzip's non-unique file extension
                } elsif ($ext_type eq $file_type) {
                    ## OK; this is the expected scenario
                } else {
                    ## mismatched intentions?
                    $comp_err = "WARNING: actual compression format '$file_type' is not compatible with file extension '$fext'!\n";
                }
                
            }
            
        } elsif ($filemsg eq 'data') {
            
            ## unknown file format: no handling
            ## no error either, unless it looks like it SHOULD be compressed...
            $comp_err = "Compression format not recognized by 'file'!  Archive is probably corrupted." if $ext_type;
            
        } else {
            
            ## some other kind of file; no handling or error
            
        }
        
    } else {
        
        ## FILE DOES NOT EXIST YET
        ## no errors; have to believe $ext_type since no $file_type
        if ($ext_type) {
            ## file will be compressed if $ext_type is valued
            $is_comp = 1;
        } else {
            ## file will NOT be compressed
            $ext_type = 'text';
        }
        @iocmd = ($typedata{$ext_type}{READ}, $typedata{$ext_type}{WRITE});
        $iocmd[1] = $typedata{bgzip}{WRITE} if $ext_type eq 'gzip' && $bgzip_flag;   # switch write command, if requested
        
    }
    
    ## Return: 1|0 if valid compressed, error if compression problem, linux commands to read/write file
    return [$is_comp, $comp_err, @iocmd];
    
}


sub open2 {
    
    ## open() wrapper with handling for compressed files, stdio
    ## DOES NOT SUPPORT bidirectional file access, e.g. '+<', 'r+'
    ## OUTPUT: a filehandle
    ## INPUT: array with 3-4 elements:
    ## 1. mode: one of R, W, A (read, write, append), or a &open descriptor (like '-|', '>', etc)
    ## 2. file: a filename, or a command if opening a pipe, or undef.  Most compressed files are recognized and handled.
    ## 3. label: an optional file label.
    ## 4. bgzip: OPTIONAL, 1|0 value.  Read/Write compressed using 'bgzip' instead of 'gz'
    ## NOTE: 'append' mode not compatible with compressed files
    ## NOTE: 'bgzip' is ignored if no file is given
    ## handles gzip, bgzip, bzip2, and zip compression formats
    
    my ($mode, $file, $label, $bgzip) = @_;
    
    $label = $label ? "$label '$file'" : "'$file'";  # upgrade $label
    my %rwamodes = ('r'=>'<', 'w'=>'>', 'a'=>'>>');
    $mode = $rwamodes{"\L$mode"} if $mode =~ /^[rwa]$/i;
    
    my ($exists, $iscomp, $comperr, $readcmd, $writecmd);
    if ($file) {
        $exists = -e $file;
        ($iscomp, $comperr, $readcmd, $writecmd) = @{ &is_compressed($file, $bgzip) };
    }
    
    my ($IN, $OUT);
    if ($mode eq '-' || $mode eq '<-') {
        ## STDIN
        if ($file) {
            die "$0: do not specify a file when requesting STDIN access!\n";
        } else {
            open($IN, $mode) or die "$0: cannot open STDIN: $!\n";
        }
    } elsif ($mode eq '>-') {
        ## STDOUT
        if ($file) {
            die "$0: do not specify a file when requesting STDOUT access!\n";
        } else {
            open($OUT, $mode) or die "$0: cannot open STDOUT: $!\n";
        }
    } elsif ($mode eq '-|') {
        ## pipe in  ($file should be a command)
        if ($file) {
            open($IN, $mode, $file) or die "$0: cannot open pipe from $label: $!\n";
        } else {
            die "$0: must specify a command when requesting pipe access!\n";
        }
    } elsif ($mode eq '|-') {
        ## pipe out  ($file should be a command)
        if ($file) {
            open($OUT, $mode, $file) or die "$0: cannot open pipe to $label: $!\n";
        } else {
            die "$0: must specify a command when requesting pipe access!\n";
        }
    } elsif ($mode eq '<') {
        ## read (file)
        if ($file) {
            if ($iscomp) {
                open($IN, $readcmd) or die "$0: problems executing '$readcmd': $!\n";
            } else {
                open($IN, $mode, $file) or die "$0: cannot read from file $label: $!\n";
            }
        } else {
            die "$0: must specify a file when requesting a '$mode' handle!\n";
        }
    } elsif ($mode eq '>' || $mode eq '>>') {
        ## write or append (file)
        my $blurb = $mode eq '>' ? 'write' : 'append';
        if ($file) {
            if ($iscomp) {
                die "$0: cannot append to compressed files!  Please use another method.\n" if $mode eq '>>';
                open($OUT, $writecmd) or die "$0: problems executing '$writecmd': $!\n";
            } else {
                open($OUT, $mode, $file) or die "$0: cannot $blurb to file $label: $!\n";
            }
        } else {
            die "$0: must specify a file when requesting a '$mode' handle!\n";
        }
    } else {
        die "$0: 'open2' access mode value '$mode' not supported!\n";
    }
    
    if ($IN) {
        return $IN;
    } elsif ($OUT) {
        return $OUT;
    }
}


sub generic_io {
    
    ## designed for simplistic situations with very direct, single-input -> single-output workflows, e.g. pipes
    ## opens (and returns) input and output handles
    ## resolves any conflict between $input and $stdio
    ## chooses output mode based on $output
    ## ensures handles are openable
    ## 
    ## usage: "my ($IN, $OUT) = &generic_io($input, $output, $stdio);", where args are from GetOptions, @ARGV, etc
    ##
    ## NOTE: be sure that print-to-screen messages in your script go to STDERR, since $output may be STDOUT and thus captured elsewhere.
    
    my ($INPUT, $OUTPUT, $STDIO) = @_;
    
    my ($IN, $OUT);
    if ($INPUT) {
        die "$0: do not specify both '-' and '-i'!\n" if $STDIO;
        $IN = &open2('R', $INPUT, 'input VCF');
    } elsif ($STDIO) {
        $IN = *STDIN;
    } else {
        die "$0: no input specified!\n";
    }
    if ($OUTPUT) {
        $OUT = &open2('W', $OUTPUT, 'output VCF');
    } else {
        $OUT = *STDOUT;
    }
    
    return ($IN, $OUT);
    
}


sub backticks_bash {
    
    ## thanks to http://www.perlmonks.org/?node_id=691255
    
    ## need to integrate this? $stderr = tee_stderr \&code;
    
    my ($cmd, $login) = @_;
    my $pipe;
    if ($login) {
        open($pipe, '-|', '/bin/bash', '--login -c', $cmd) or return;
    } else {
        open($pipe, '-|', '/bin/bash', '-c', $cmd) or return;
    }
    local $/ = wantarray ? $/ : undef;
    <$pipe>
        
}


sub execute {
    
    ## system() wrapper that sends command to screen, optional logfile, and runs command
    ## INPUT: array with 4 elements:
    ## 1. command to run
    ## 2. visualize command on std-what? (1=out,2=err)
    ## 3. optional: path to logfile, or open filehandle
    ## 4. optional: capture messages from std-what? (1=out,2=err,3=both)  # BE SURE THIS DOES NOT CONFLICT WITH COMMAND OUTPUT
    ## 5. optional: 1|0; execute with 'bash -c'
    
    my ($cmd, $std, $fh, $cap, $bash) = @_;
    
    chomp $cmd;  # just in case
    
    if ($std == 1) {
        print STDOUT "$cmd\n";
    } elsif ($std == 2) {
        print STDERR "$cmd\n";
    }
    $cmd =~ s/^\n?RUNNING: //;  # goofy way to allow "headerization" of commands; later versions will be more formal
    
    if ($fh) {
        if (defined fileno($fh)) {
            ## filehandle
            print $fh "$cmd\n";
        } else {
            ## file
            open my $LOG, '>>', $fh or die "$0: Failed to append to log file '$fh': $!\n";
            print $LOG "$cmd\n";
            close $LOG;
        }
    }
    
    my $cmd2_sh = sub { system($cmd) };
    my $cmd2_bash = sub { system("bash -c '$cmd'") };
    my $cmd2 = $bash ? $cmd2_bash : $cmd2_sh;
    
    my $msg;
    if ($bash && $cmd =~ /'/) {
        $msg = &backticks_bash($cmd);
    } else {
        if ($cap == 1) {
            $msg = $std ? tee_stdout \&$cmd2 : capture_stdout \&$cmd2;
        } elsif ($cap == 2) {
            $msg = $std ? tee_stderr \&$cmd2 : capture_stderr \&$cmd2;
        } elsif ($cap == 3) {
            $msg = $std ? tee_merged \&$cmd2 : capture_merged \&$cmd2;
        } else {
            &$cmd2;
        }
    }

    if ($msg =~ /\S/) {
        if ($fh) {
            if (defined fileno($fh)) {
                ## filehandle
                print $fh "$msg\n";
            } else {
                ## file
                open my $LOG, '>>', $fh or die "$0: Failed to append to log file '$fh': $!\n";
                print $LOG "$msg\n";
                close $LOG;
            }
        }
    }
    
}





sub logreport {
    
    ## prints a message to STDOUT and to a given logfile 
    
    my ($msg, $out, $fh, $nochomp) = @_;
    
    chomp $msg unless $nochomp;  # so you don't have to remember if it needs newlines or not
    
    if (defined fileno($fh)) {
        ## $fh is a handle
        print $fh "$msg\n";
    } else {
        ## $fh is a logfile
        open my $LOG, '>>', $fh or die "$0: Failed to append to log file '$fh': $!\n";
        print $LOG "$msg\n";
        close $LOG;
    }
    
    if ($out == 1) {
        print STDOUT "$msg\n";
    } elsif ($out == 2) {
        print STDERR "$msg\n";
    }
}





sub print_table {
    
    ## print array-of-arrays object to a filehandle using a delimiter, or sprintf-ed
    ## input: 1=reference to array-of-arrays, 2=1|0 use sprintf, 3=filehandle, 4=delimiter string
    ## writes directly to filehandle
    ## only understands three column datatypes: string, integer and float
    ## strings columns are left-justified, numeric columns are right-justified
    ## float formatting specifies number of digits to the right of decimal by the first float value seen per column
    ##  - thus float values should be sprintf-ed BEFORE being stored in the object, to avoid formatting inconsistencies
    
    my ($ref, $do_sprintf, $fh, $delim) = @_;
    $fh = *STDOUT unless $fh;
    my @AOA = @$ref;
    my @rows = (0..$#AOA);
    my @cols = (0..$#{ $AOA[0] });
    
    if ($do_sprintf) {
        
        $delim = '  ' unless defined $delim;  # default sprintf delimiter: two spaces
        my $fmt;
        
        ## first pass: get column widths and numeric/non status
        ## construct format on the fly
        foreach my $col (@cols) {
            my ($width, $right);
            my $numeric = 1;  # initial assumption
            my $float = 0;    # initial assumption
            foreach my $row (@rows) {
                my $w = length($AOA[$row][$col]);
                $width = $w if $w > $width;
                next if $row == 0;    # don't format-test the header row
                next unless $numeric; # still have breakable assumption(s)
                if ($AOA[$row][$col] =~ /[^0-9.-]/) {
                    $numeric = $float = 0;  # column is not numeric
                } elsif ($AOA[$row][$col] =~ /\./) {
                    $float = 1;       # irrelevant unless $numeric
                    $right = length( (split /\./, $AOA[$row][$col])[-1] );
                }
            }
            $fmt .= '%'.($numeric?'':'-').$width.($numeric?($float?".${right}f":'i'):'s');
            $fmt .= $col == $#cols ? "\n" : $delim;
        }
        (my $hfmt = $fmt) =~ s/(i|\.\d+f)(\s)/s$2/g;  # header format: convert numeric 'i' or '.#f' to string 's'; preserve following space OR newline
        $hfmt =~ s/\%([^-])/%-$1/g;   # ensure headers all left-justify
        
        ## second pass: print
        print $fh sprintf($hfmt, @{ $AOA[0] });
        foreach my $row (1..$#rows) {
            print $fh sprintf($fmt, @{ $AOA[$row] });
        }
        
    } else {
        
        $delim = "\t" unless defined $delim;  # default non-sprintf delimiter: tab
        
        foreach my $row (@rows) {
            print $fh join($delim, @{ $AOA[$row] }), "\n";
        }
        
    }
    
}





sub nfsAnalyze {
    
    ## reports on undeletable files (i.e. files that, when deleted, result in undeletable .nfs files).
    ## for troubleshooting temp directories that cannot be removed by the script.
    
    my $dir = shift;
    my @files = glob "$dir/*";
    my $maxw;
    foreach (@files) {
        $maxw = length($_) if length($_) > $maxw;
    }
    my $fmt = 'FILE: %-'.$maxw.'s | %s'."\n";
    my %already;
    
    print STDERR "\n\n.NFS FILE ANALYSIS\nTHIS PID: $$\nCONTENTS OF '$dir':\n";
    foreach my $file (@files) {
        my $fuser = `fuser $file 2>/dev/null`;
        $fuser =~ s/\s//g;
        my $msg = $fuser ? "IN USE: $fuser" : 'not in use';
        system "rm -f $file";
        my @nfs = glob "$dir/.nfs*";
        foreach my $nfs (@nfs) {
            next if $already{$nfs};  # some previous file's .nfs file
            $already{$nfs} = 1;      # this is a new .nfs file
            $msg .= " | NOT DELETABLE: $nfs";
        }
        print STDERR sprintf($fmt, $file, $msg);
    }
    print "\n\n";
    
}





sub validate_genoanno {
    
    my ($GENO, $ANNO, $STARLEN) = @_;
    
    my $Gdir  = "$indexes/$GENO";
    my $Gpref = "$indexes/$GENO/$GENO";
    my $Adir  = "$indexes/$GENO/$ANNO";
    my $Apref = "$indexes/$GENO/$ANNO/$GENO.$ANNO";
    
    die "$0: Expected dataset path '$Gdir' does not exist!\n" unless -d $Gdir;
    die "$0: Expected dataset path '$Adir' does not exist!\n" if $ANNO && ! -d $Adir;
    
    my @out = ($Gdir, $Gpref, $Adir, $Apref);
    
    if ($STARLEN) {
		## if $STARLEN (read length), then must be an $ANNO STAR idx, because $geno idx has no indicated read length
		my $Sdir = "$indexes/$GENO/$ANNO/STAR_${STARLEN}bp";
		die "$0: Expected STAR index '$Sdir' does not exist!\n" if $STARLEN && ! -d $Sdir;
		push @out, $Sdir;
    }
    
    return @out;
}





sub validate_bam {
    
    ## Returns non-null if the bam file is corrupted
    
    my $BAM = shift;  # /path/to/file.bam
    my $TMP = "validate_bam.$$.tmp";
    system "samtools view $BAM 2> $TMP | head > /dev/null";
    chomp(my $ERR = `cat $TMP`);
    system "rm -f $TMP";
    return $ERR;
}





sub timestamp {
    
    ## returns date, time, or date+time timestamps
    
    my $timetype = shift;    # FULL = date + time | DATE = date | TIME = time
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
        #    if ($seconds > $timescales{HOUR}) {
        $timediff{HOUR} = $seconds/$timescales{HOUR};
        $seconds %= $timescales{HOUR};
        #    }
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
    
    my ($pid, $platform, $file, $msg, $script) = @_;    # $msg, $script, $file can be null
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
        $stdfreq *= $freqs{$letter} foreach (1..$nondeg{W}{$letter});    # weighted product of frequencies for non-degenerates * number of occurrances
    }
    foreach my $letter (keys %{ $nonstd{W} }) {
        $nonfreq *= $freqs{$letter} foreach (1..$nonstd{W}{$letter});    # weighted product of frequencies for non-standards * number of occurrances
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
    
    my %CODONS = %{ &translate({'sequence'=>$NTSEQ, 'mode'=>$tmode, 'ascodons'=>1}) };   # $NTSEQ, $tmode, $stats, $keepcase, $orfs, $stopchar, $ascodons
    
    my (%COM, %STOP);
    foreach my $frame (sort keys %CODONS) {
        my (%NUM, %PSF, %AG, %CGA);
        my ($NUMsum, %PSFsum, %AGsum, %CGAsum, @PSFfail);
        #    print "$frame: @{ $CODONS{$frame} }\n";
        #    print "$N, $C, ", (scalar @{ $CODONS{$frame} }), "\n";
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
            #        print "$frame: ($AG{A}*$AG{G}) / ($CGA{C}*$CGA{G}*$CGA{A}+$STOP{$frame}+$W) = $COM{$frame}\n";
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
                #        print "$frame: $codon: $NUM{$codon} / $DEN{$codon}\n";
                #        print "$frame: $codon: $NUM{$codon} / $DEN{$codon} = ", ($NUM{$codon} / $DEN{$codon}), "\n";
                $COM{$frame} += $NUM{$codon} / $DEN{$codon};
            }
            #        print "$frame: COM = $COM{$frame}\n";
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
    my $SEQ = $data[0] =~ /^SCALAR\(0x[a-f0-9]+\)$/ ? ${$data[0]} : $data[0];
    my $WIDTH = $data[1] ? $data[1] : 50;  # default width
    
    my (@lines, $start);
    my $blocks = length($SEQ) / $WIDTH;
    $blocks++ if (length($SEQ) % $WIDTH != 0);
    foreach (1..$blocks) {
        push @lines, substr($SEQ, $start, $WIDTH);
        $start += $WIDTH;
    }
    #    my $seqblock = join("\n", @lines)."\n";
    my $seqblock = join("\n", @lines);
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
    ## Designed with UCSC/Ensembl chromosome naming conventions in mind; also works with some others.
    ## Output order is basically: 1. numerics (incl. Romans), 2. alphas, 3. randoms, 4. scaffolds, 5. unplaced, 6. mito, 7. else
    
    ## FIXME ASAP:
    ## Up-front post-suffix detection like for '_random', 'Het', etc.
    ## Thence to suffix tranches, unless interleaving.
    ## Also ensure interleavable tranches sort in the same chrom order as the canonicals!!
    
    # $dataref = chr names array ref
    # $interleave: 0 = sort randoms into their own tranch, 1 = sort randoms next to their canonical counterparts
    
    my ($dataref, $interleave, $scaf_limit, $verbose, $troubleshoot) = @_;
    my $scaf_limit = 5 unless defined $scaf_limit;  # if any non-chr prefix exceeds this, it will be classified as a scaffold
    
    my %troubleshooting;
    my %chr_prefixes = map {($_=>1)} qw/ chr group /;            # if they have any prefixes
    my %mito_names = map {($_=>1)} qw/ M MT chrM chrMT MtDNA /;  # or otherwise matches /mito/i
    my %plastid_names = map {($_=>1)} qw/ PT /;                  # or otherwise matches /(chloro|plast)/i
    my %plasmid_names = map {($_=>1)} qw/ 2-micron 2micron /;    
    my %unplaced_names = map {($_=>1)} qw/ U Un chrU chrUn /;    # or otherwise matches /^un/i
    #my $Bclass = '\d._-';
    
    ## Set some hardcoded sorting patterns for Drosophila
    my $dros_grep1 = '^(2L|2R|3L|3R)';        # no longer including 'Het'
    my $dros_grep2 = '(^|\D)(2L|2R|3L|3R)$';  # "
    
    ## Name-test function
    sub name_type_test {
        my ($name, $hashref, $regex) = @_;
        my $matched;
        if ($regex && $name =~ /$regex/i) {
            $matched = 1;
        } else {    
            foreach my $key (keys %$hashref) {
                if ($name =~ /^$key$/i) {
                    $matched = 1;
                    last;
                }
            }
        }
        return $matched;
    }
    
    ## Separate sequences into tranches based on names
    my $inputs = scalar @$dataref;
    #print "INPUTS: @$dataref\n" if $troubleshoot;
    
    ## First pass: run obvious sequence identity tests first
    ## Separate "early" tranches
    print "Stage 1\n" if $verbose;
    my (%tranches, @remainder);
    foreach my $name (@$dataref) {
        if ( &name_type_test($name, \%mito_names, 'mito') ) {
            ## Mitochondrion
            push @{ $tranches{mito} }, $name;
            $troubleshooting{$name} .= "Early=mito; ";
        } elsif ( &name_type_test($name, \%plastid_names, '(chloro|plast)') ) {
            ## Plastid / chloroplast
            push @{ $tranches{plastid} }, $name;
            $troubleshooting{$name} .= "Early=plastid; ";
        } elsif ( &name_type_test($name, \%plasmid_names) ) {
            ## Plasmids
            push @{ $tranches{plasmid} }, $name;
            $troubleshooting{$name} .= "Early=plasmid; ";
        } elsif ( &name_type_test($name, \%unplaced_names, '^(chr)?(u)') ) {
            ## Unknown/unplaced sequences
            push @{ $tranches{unk} }, $name;
            $troubleshooting{$name} .= "Early=unk; ";
        } else {
            push @remainder, $name;
        }
    }
    @remainder = sort @remainder;  # critical for later ops
    #print "$_\t".scalar(@{ $tranches{$_} })."\n" foreach sort keys %tranches;
    
    ## Second pass: find obvious prefixes
    ## Prefixes are discovered by taking the longest non-numeric strings from the left.
    ## e.g. 'chr2' -> 'chr' , 'scaffold_1000' -> 'scaffold_' , 'DS484752.1' -> 'DS'
    print "Stage 2\n" if $verbose;
    my %prefixes;
    foreach my $name (@remainder) {
        if ($name =~ /^(chr)(.*)$/) {
            $prefixes{$1} = 0;
        } elsif ($name =~ /^([^\d]+)(\d[^A-Za-z]+)$/) {
            $prefixes{$1} = 0;
        }
    }
    
    ## Third pass: classify sequences by prefix
    ## MATCH LONGER PREFIXES FIRST
    ## Also identify potential Drosophila chromosomes
    ## Also identify potential Roman-numeral chromosomes:
    ##   Chr suffix (or entire name) must be a roman numeral, BUT NOT 'M', THAT IS MITO (and who uses Romans up to 1000 anyway???)
    ##   However 'X' can be counted, since 'chrX' could actually be chr10.
    ##   Known Roman users: Saccharomyces, Caenorhabditis, stickleback, pombe
    print "Stage 3\n" if $verbose;
    my (%byPrefix, %prefix_lookup, @noPrefix, %droschrs, %romanchrs);
    my @ordpref = sort { length($b) <=> length($a) } keys %prefixes;
    foreach my $name (@remainder) {
        my $has_prefix = 0;
        foreach my $prefix (@ordpref) {
            if ($name =~ /^$prefix(.*)/) {
                my $suffix = $1;
                $byPrefix{$prefix}{$suffix} = $name;  # store by suffix, for sorting
                $prefix_lookup{$name} = $prefix;
                $prefixes{$prefix}++;
                $has_prefix = 1;
                $troubleshooting{$name} .= "Prefix=$prefix; Suffix=$suffix; ";
                if ($name =~ /$dros_grep1/i || $name =~ /$dros_grep2/i) {
                    $droschrs{$prefix}{$suffix} = $name;
                    $troubleshooting{$name} .= "maybe dros ($suffix); ";
                }
                if ($suffix ne 'M' && isroman($suffix)) {
                    $romanchrs{$prefix}{$suffix} = $name;
                    $troubleshooting{$name} .= "maybe roman ($suffix); ";
                }
                last;   # first (longest) match is what we take
            }
        }
        unless ($has_prefix) {
            my $prefix = '';
            $byPrefix{$prefix}{$name} = $name;   # store by suffix, for sorting
            if ($name =~ /$dros_grep1/i || $name =~ /$dros_grep2/i) {
                $droschrs{$prefix}{$name} = $name;
                $troubleshooting{$name} .= "maybe dros ($name); ";
            }
            if ($name ne 'M' && isroman($name)) {
                $romanchrs{$prefix}{$name} = $name;
                $troubleshooting{$name} .= "maybe roman ($name); ";
            }
        }
    }
    
    #print Dumper(\%byPrefix),"\n" if $troubleshoot;
    
    ## Fourth pass:
    ## Cull minor prefixes and reassign contents to the null prefix
    print "Stage 4\n" if $verbose;
    ## First: what is the smallest acceptable N chroms for a prefix to be considered valid?
    my $min_prefix = $inputs > 100 ? 5 : 0;
    foreach my $prefix (keys %byPrefix) {
        if (scalar(keys %{ $byPrefix{$prefix} }) < $min_prefix) {
            print "CULL: $prefix\n" if $verbose;
            ## Too small to be a "real" sequence prefix
            $byPrefix{''}{$_} = $byPrefix{$prefix}{$_} foreach keys %{ $byPrefix{$prefix} };
            $romanchrs{''}{$_} = $romanchrs{$prefix}{$_} foreach keys %{ $romanchrs{$prefix} };
            $droschrs{''}{$_} = $droschrs{$prefix}{$_} foreach keys %{ $droschrs{$prefix} };
            foreach my $suffix (keys %{ $byPrefix{$prefix} }) {
                my $name = $byPrefix{$prefix}{$suffix};
                $troubleshooting{$name} .= "Prefix culled; Prefix=''; ";
            }
            delete $byPrefix{$prefix};
            delete $romanchrs{$prefix};
            delete $droschrs{$prefix};
        }
    }
    
    ## Fifth pass:
    ## Decide if Romans were really Roman, on a prefix-by-prefix basis
    ## Decide if Drosophila chrs are really pertinent
    print "Stage 5\n" if $verbose;
    my (%has_roman, %has_dros);
    foreach my $prefix (keys %romanchrs) {
        if (scalar(keys %{ $romanchrs{$prefix} }) >= 3) {
            ## at least 3, and you are considered Roman
            $has_roman{$prefix} = 1;
            foreach my $suffix (keys %{ $byPrefix{$prefix} }) {
                my $name = $byPrefix{$prefix}{$suffix};
                $troubleshooting{$name} .= "roman confirmed; " if $romanchrs{$prefix}{$suffix};
            }
        }
    }
    foreach my $prefix (keys %droschrs) {
        if (scalar(keys %{ $droschrs{$prefix} }) >= 3) {
            ## at least 3, and you are considered Drosophila
            $has_dros{$prefix} = 1;
            foreach my $suffix (keys %{ $byPrefix{$prefix} }) {
                my $name = $byPrefix{$prefix}{$suffix};
                $troubleshooting{$name} .= "dros confirmed; " if $droschrs{$prefix}{$suffix};
            }
        }
    }
    
    ## Sixth pass: 
    ## Find any *_random/*Het/etc where "*" is an independent chromosome name.
    ## Thus, say "chr1" will be the "interleaver", and "chr1_random" will be the "interleavee".
    print "Stage 6\n" if $verbose;
    my (%interleavees, %inter_lookup, %bySuffix);
    if ($interleave) {
        foreach my $prefix (keys %byPrefix) {
            my @names = sort values %{ $byPrefix{$prefix} };  # MUST BE SORTED!!  Must ensure ordering like "chr1" < "chr1_random".
            foreach my $i (0..$#names-1) {
                foreach my $j ($i+1..$#names) {
                    if ($names[$j] =~ /^$names[$i]([.:_-])?([^A-Za-z0-9]|Het|alt|random)(.*)/) { 
                        $inter_lookup{I2J}{$names[$i]}{$names[$j]} = "$2$3";  # "$2$3" is everything after the interleaver's name, e.g. "_random" or "Het" etc.
                        $inter_lookup{J2I}{$names[$j]}{$names[$i]} = "$2$3";
                        $interleavees{$names[$j]} = 1;  # these $names[$j] will sort subordinate to some $names[$i]
                        #$bySuffix{$1}{$names[$j]} = $names[$i];
                        $troubleshooting{$names[$j]} .= "can interleave with $names[$i]; ";
                    }
                }
            }
        }
    }
    
    ## Seventh pass:
    ## Separate chroms into final tranches
    ## Make sure Drosophila main chroms go into {dros} which will print first, and not {mixed} which will print last!
    print "Stage 7\n" if $verbose;
    my (%byTranch, %interleaved);
    my %byTranch = ();
    my %interleaved = ();
      foreach my $prefix (keys %byPrefix) {
        foreach my $suffix (keys %{ $byPrefix{$prefix} }) {
            my $name = $byPrefix{$prefix}{$suffix};
            next if $interleavees{$name};   ## these will sort subordinate to something else
            
            ## Assign tranches
            my $tranch;
            if (!$interleave && exists $inter_lookup{J2I}{$name}) {
                print "AUTO-MIX: $name\n";
                ## Interleavable with more-canonical-sounding chrom name, but !$interleave:
                ## Automatic 'mixed' assignment
                $tranch = 'mixed';
            } elsif ($has_dros{$prefix} && $name =~ /Het/) {
                print "DROS-MIX: $name\n";
                ## Drosophila 'Het' chr -- current handling for these is very wonky IMO and indicative of a larger suffix-detection issue which is not being addressed...
                $tranch = 'mixed';
            } elsif ($has_roman{$prefix} && $romanchrs{$prefix}{$suffix}) {
                ## Roman tranch
                ## potential issue: tranch spuriously defined as Roman, and this $suffix was not Roman
                $tranch = 'roman';
                $suffix = arabic($suffix);  # REPLACE $suffix with sortable arabic-numeral equivalent
            } elsif ($has_dros{$prefix} && $droschrs{$prefix}{$suffix}) {
                ## Drosophila tranch
                $tranch = 'dros';
            } elsif (looks_like_number($suffix)) {
                ## Numeric tranch
                $tranch = 'numeric';
            } elsif ($suffix =~ /^[A-Za-z]+$/) {
                ## Alphabetical tranch
                $tranch = 'alpha';
            } elsif ($name =~ /_(het|hap|alt|random|patch)\b/) {
                
                ## FIXME: FIND SOMETHING ELSE TO DO WITH THESE
                
                ## Mixed tranch
                $tranch = 'mixed';
            } else {
                ## Mixed tranch
                $tranch = 'mixed';
            }
            
            ## Store chrom by prefix, tranch, and sortable suffix
            $byTranch{$prefix}{$tranch}{$suffix} = $name;
            $troubleshooting{$name} .= "tranch=$tranch; ";
            
            ## If $name is an interleaver, then add interleavees
            if ($interleave && exists $inter_lookup{I2J}{$name}) {
                my $interleave_i = 0;
                foreach my $jname (sort keys %{ $inter_lookup{I2J}{$name} }) {
                    $byTranch{$prefix}{$tranch}{ $suffix.'.'.$interleave_i++ } = $jname;  # .i suffix ensures these sort after main chr
                    $troubleshooting{$jname} .= "interleaved with $name; ";
                }
            }
        }
    }
    
    
    ## Tranch sort function
    sub sort_tranches {
        my $ref = shift;
        my @sorted;
        foreach my $tranch (qw/ dros numeric roman alpha mixed /) {
            if (exists $$ref{$tranch}) {
                if ($tranch eq 'numeric' || $tranch eq 'roman') {
                    push @sorted, map { $$ref{$tranch}{$_} } sort { $a <=> $b } keys %{ $$ref{$tranch} };
                } else {
                    push @sorted, map { $$ref{$tranch}{$_} } sort { $a cmp $b } keys %{ $$ref{$tranch} };
                }
            }
        }
        return @sorted;
    }
    
    ## Eigth pass: final sorting and output
    print "Stage 8\n" if $verbose;
    my (@final, $haschr);
    ## Do unprefixed names actually seem to contain real chrom names, or only {mixed} scaffold junk?
    ## To reliably detect unprefixed chromosomes, they must have been sorted into a *stringent* tranch, like numeric OR alpha OR roman.
    my $null_prefix_meaningful = exists $byTranch{''} && (exists $byTranch{''}{numeric} || exists $byTranch{''}{alpha} || exists $byTranch{''}{roman});
    
    ## Sort (obvious) chromosomes first
    if (exists $byTranch{chr}) {
        ## Most-obvious chromosome prefix (typical UCSC)
        $haschr = 1;
        push @final, &sort_tranches($byTranch{chr});
        delete $byTranch{chr};
    } elsif (exists $byTranch{group}) {
        ## Rare chromosome prefix
        $haschr = 1;
        push @final, &sort_tranches($byTranch{group});
        delete $byTranch{group};
    } elsif ($null_prefix_meaningful) {
        ## Ensembl chromosome prefix (no prefix)
        $haschr = 1;
        push @final, &sort_tranches($byTranch{''});
        delete $byTranch{''};
    }
    
    ## Sort (obvious) scaffolds second
    if (exists $byTranch{scaffold}) {
        ## Most-obvious scaffold prefix
        $haschr = 1;
        push @final, &sort_tranches('scaffold');
        delete $byTranch{scaffold};
    }
    ## Sort remaining prefixes by the number of sequences in each, decreasing
    foreach my $prefix (sort { scalar(keys %{ $byTranch{$b} }) <=> scalar(keys %{ $byTranch{$a} }) } keys %byTranch) {
        push @final, &sort_tranches($prefix);
        delete $byTranch{$prefix};
    }
    
    ## Add "early" tranches, in desired order
    foreach my $early (qw/ unk mito plastid plasmid /) {
        push @final, @{ $tranches{$early} } if exists $tranches{$early};
    }
    
    ## Test that input == output
    my $outputs = scalar(@final);
    if ($inputs != $outputs) {
        my $msg = "$0 &chrsort: $outputs outputs != $inputs inputs\n";
        if ($troubleshoot) {
            print $msg;
        } else {
            die $msg;
        }
    }
    
    #print Dumper(\%byPrefix),"\n" if $troubleshoot;
    
    ## All done!
    print "Post-production\n" if $verbose;
    if ($troubleshoot) {
        my (@final2, %already);
        foreach my $name (@final) {
            push @final2, "$name : $troubleshooting{$name}";
            $already{$name} = 1;
        }
        foreach my $name (@$dataref) {
            next if $already{$name};
            push @final2, "$name : LOST : $troubleshooting{$name}";
        }
        return \@final2;
    } else {
        return \@final;
    }
    
}





sub readFile {
    
    ## reads a text file into a standard format (array-of-arrays or hash-of-arrays) with various options
    
    # $filename = filename | $delim is the file delimiter ("\t", ",", etc) | $chomp = [01] remove terminal [\r\n]+ 
    # $header = integer: capture $header first lines as header | $skip = integer: skip $skip first lines
    # $nfields = integer: return row broken into $nfields fields (as "split /$delim/, $_, $nfields")
    # $keycol = integer: use data in column $keycol (USE 0-BASED!!!) as hash key for storing row (changes output from @-of-@ to %-of-@)
    my ($filename, $delim, $chomp, $header, $skip, $nfields, $keycol) = @_;
    $chomp = 1 if !defined $chomp;
    my (@header, @data, %data);
    open my $IN, '<', $filename or die "Can't open '$filename': $!\n";
    while (<$IN>) {
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
    close $IN;
    if ($keycol) {
        return(\%data, \@header);
    } else {
        return(\@data, \@header);
    }
}





sub writeFile {
    
    ## writes a text file from a contents block (must be single string, so join "\n" any arrays)
    
    # $filename = filename | $contents = string to write | $mode = '>' or '>>'
    my ($filename, $contents, $mode) = @_;
    $mode = '>' unless $mode;
    open my $OUT, $mode, $filename or die "Can't write to '$filename': $!\n";
    print $OUT $contents;
    close $OUT;
}





sub readFasta {
    
    ## reads a fasta into a hash (keys=headers, values=sequence)
    
    # $fh = filehandle | $chomp = 0|1 remove terminal [\r\n]+ (affects sequence block only)
    my ($fh, $nochomp, $toupper) = @_;
    my ($header, @headers, %data);
    while (<$fh>) {
        next if $_ =~ /^#/;
        $_ =~ s/[\n\r]+$// unless $nochomp;
        if ($_ =~ /^>(.*)/) {
            $header = $1;
            push @headers, $header;  # record of header order
        } else {
            $data{$header} .= $toupper ? "\U$_" : $_;
        }
    }
    return(\%data, \@headers);
}





sub writeFasta {
    
    ## writes a sequence hash (keys=headers, values=sequence) to a fasta file
    
    # $filename = filename | $contents = string to write | $mode = '>' or '>>' | $width = line width
    my ($filename, $dataref, $mode, $width) = @_;
    $mode = '>' unless $mode;
    $width = 50 unless $width;
    open my $OUT, $mode, $filename or die "Can't write to '$filename': $!\n";
    print $OUT ">$_\n", ${ &blockify($$dataref{$_}, $width) } foreach keys %$dataref;
    close $OUT;
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
                $start = index($consensus, '[', 0);        # keep restarting from 0 because matches get eliminated
                last if $start == -1;
                my $end = index($consensus, ']', $start);
                substr($consensus, $start, $end-$start+1, $dchars[$idx]);
                $idx++;
            }
            $mdata{$reflag}{CONSENSUS} = $consensus;
            $reflag = undef;
        } elsif ($_ =~ /^\s+Motif (\d+) regular expression/) {
            $reflag = $1;
        } elsif ($_ =~ /\s+Motif (\d+) sites sorted by position p-value/) {
            $mbdflag = $1;    # prep for capture
        } elsif ($mbdflag && $_ !~ /^-/) {
            next if $_ =~ /^Sequence name\s/;
            $mbdflag2 = 1;    # begin capture
            my ($seq, $strand, $start, $pval, $etc) = split /\s+/, $_, 5;
            $mdata{$mbdflag}{POS}{$seq}++;
            $mdata{$mbdflag}{PVAL}{$pval}++;
        } elsif ($mbdflag2 && $_ =~ /^-/) {
            $mbdflag = $mbdflag2 = undef;    # end capture
        } elsif ($_ =~ /^MOTIF\s+(\d+)\s+width =\s+(\d+)\s+sites =\s+(\d+)\s+llr =\s+(\d+)\s+E-value =\s+(\S+)/) {
            $motif = $1;
            $mdata{$motif}{WIDTH} = $2;
            $mdata{$motif}{NSITES} = $3;
            $mdata{$motif}{LOGLR} = $4;
            $mdata{$motif}{EVAL} = $5;
        } elsif ($_ =~ /^\s+Motif (\d+) position-specific probability matrix/) {   # begin matrix capture
            ($pssmflag, $dashed, $line) = ($1, 0, 0);
            $mdata{$pssmflag}{PSSM} = {'A'=>[], 'C'=>[], 'G'=>[], 'T'=>[]};  # initialize
        } elsif ($dashed == 2) {    # second dashed line: terminate PSSM capture
            $dashed = $pssmflag = undef;
        } elsif ($pssmflag) {        # capturing lines; $pssmflag = motif ID
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
        $temp{$motif}{obs_sites} += $mdata{$motif}{POS}{$_} foreach keys %{ $mdata{$motif}{POS} };
        if ($mdata{$motif}{NSITES} != $temp{$motif}{obs_sites}) {
            print "get_memedata error: motif $motif: stated sequence count '$mdata{$motif}{NSITES}' != observed sequence count '$temp{$motif}{obs_sites}'!\n";
        }
    }
    
    my ($inflag, $enflag, $motif);
    open IN, $memehtml;
    while (<IN>) {
        $_ =~ s/[\n\r]//g;
        if ($_ =~ />Motif (\d+)\s*</i) {
            $motif = $1;
        } elsif ($_ =~ />Information Content\s*</) {
            $inflag = 1;
        } elsif ($_ =~ />Relative Entropy\s*</) {
            $enflag = 1;
        } elsif ($_ =~ /([\d.]+)_\(bits\)/) {
            if ($inflag) {
                $mdata{$motif}{INFO} = $1;
                $inflag = 0;
            } elsif ($enflag) {
                $mdata{$motif}{ENTRO} = $1;
                $enflag = 0;
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
        #    $fail = 1;
        #    print "UTR_weld error: no UTR coords supplied!  Nothing to do.\n";
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
    my @clean = $NArm ? @$ref : @{ &strip_NA($ref, $NArm) };
    
    if ($clean[0] eq 'NA') {  # had NAs, and didn't specify removal
        return 'NA';
    } else {
        my $avg;
        $avg += $_ foreach @clean;
        $avg /= scalar @clean;
        return $avg;
    }
}





sub median {
    
    ## returns the median of an array
    
    my ($ref, $NArm) = @_;
    my @clean = $NArm ? @$ref :@{ &strip_NA($ref, $NArm) };
    
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
    
    ## translates NT -> AA sequence
    ## input: hash ref (defined below)
    ## output: hash ref; hash keys = translated frames; subkeys = sequence, ORF stats, etc
    ## key-value pairs: 
    ##  sequence => NT sequence; only mandatory parameter.
    ##  stopchar => stop-codon character; default '*'.  If 'names', uses { # $ % } for { ochre amber opal }.
    ##  mode => 1|3|6|bestof3|bestof6 (default 1); frames to translate in.  Modes 'bestofN' do N-frame translation but only return the frame(s) with longest ORFs.
    ##  frames => list with any of { +0 +1 +2 -0 -1 -2 }, indicating specific frame(s) to be translated in; overrides 'mode'.
    ##  from_frame => frame that this sequence begins in; will translate only in that frame.  Overrides 'frames' and 'mode'.
    ##  orfs => 1|0 (default 0); also return longest ORF sequence(s) per frame?
    ##  stats => 1|0 (default 0); include stats with returned seq?  Returns frame-level stats by default; if 'orfs=>1' also returns ORF stats.
    ##  keepcase => 1|0 (default 0); leave AA in lowercase if codons had any lowercase NT-sequence letters?
    ##  autotrim => 1|0 (default 0); Automatically trim any trailing 1 or 2 bases off the end of each frame?
    ##  ascodons => 1|0 (default 0); return array of codons instead of amino acid sequence?  Incompatible with and overrides 'stats=>1'.
    
    my %PARAMS = %{ $_[0] };
    my %TRANS;  # return object
    
    $PARAMS{stats} = 0 if $PARAMS{ascodons};   ##########  NO STATS IN THIS MODE  ##########
    @{ $PARAMS{frames} } = ($PARAMS{from_frame}) if defined $PARAMS{from_frame};
    $PARAMS{mode} = 'FRAMES' if defined $PARAMS{frames};
    $PARAMS{mode} = 1 unless $PARAMS{mode};

    $PARAMS{stopcodons} = { 'TAA'=>'#', 'TAG'=>'$', 'TGA'=>'%' };  ## ochre, amber, opal/umber
    
    $PARAMS{aminoacids} = {
        'TTT' => 'F',    'TTC' => 'F',    'TTA' => 'L',    'TTG' => 'L',
        'TCT' => 'S',    'TCC' => 'S',    'TCA' => 'S',    'TCG' => 'S',
        'TAT' => 'Y',    'TAC' => 'Y',    'TAA' => '*',    'TAG' => '*',
        'TGT' => 'C',    'TGC' => 'C',    'TGA' => '*',    'TGG' => 'W',
        'CTT' => 'L',    'CTC' => 'L',    'CTA' => 'L',    'CTG' => 'L',
        'CCT' => 'P',    'CCC' => 'P',    'CCA' => 'P',    'CCG' => 'P',
        'CAT' => 'H',    'CAC' => 'H',    'CAA' => 'Q',    'CAG' => 'Q',
        'CGT' => 'R',    'CGC' => 'R',    'CGA' => 'R',    'CGG' => 'R',
        'ATT' => 'I',    'ATC' => 'I',    'ATA' => 'I',    'ATG' => 'M',
        'ACT' => 'T',    'ACC' => 'T',    'ACA' => 'T',    'ACG' => 'T',
        'AAT' => 'N',    'AAC' => 'N',    'AAA' => 'K',    'AAG' => 'K',
        'AGT' => 'S',    'AGC' => 'S',    'AGA' => 'R',    'AGG' => 'R',
        'GTT' => 'V',    'GTC' => 'V',    'GTA' => 'V',    'GTG' => 'V',
        'GCT' => 'A',    'GCC' => 'A',    'GCA' => 'A',    'GCG' => 'A',
        'GAT' => 'D',    'GAC' => 'D',    'GAA' => 'E',    'GAG' => 'E',
        'GGT' => 'G',    'GGC' => 'G',    'GGA' => 'G',    'GGG' => 'G' };
    
    $PARAMS{sequence} =~ s/^\s+//;  # clip leading whitespace
    $PARAMS{sequence} =~ s/\s+$//;  # clip trailing whitespace
    $PARAMS{sequence} =~ s/(.*)/\U$1/ unless $PARAMS{keepcase};     # all caps
    $PARAMS{sequence} =~ s/U/T/g;   # RNA->DNA
    warn "Sequence contains non-ACGTN characters! (\"$1\")\n" if $PARAMS{sequence} =~ /([^ACGTNacgtn])/;   # nonstandard letters exist
    
    if ($PARAMS{stopchar} eq 'names') {
        $PARAMS{aminoacids}{$_} = $PARAMS{stopcodons}{$_} foreach keys %{ $PARAMS{stopcodons} };
    } elsif ($PARAMS{stopchar} && $PARAMS{stopchar} ne '*') {
        $PARAMS{aminoacids}{$_} = $PARAMS{stopchar} foreach keys %{ $PARAMS{stopcodons} };
    }
    
    if ($PARAMS{stopchar} && $PARAMS{stopchar} ne 'names') {
        ## custom-defined stop char names: validate!  cannot overlap with AA letters
        my %AAnonStop = my %AAonlyStop = %{ $PARAMS{aminoacids} };
        delete $AAnonStop{$_} foreach keys %{ $PARAMS{stopcodons} };  ## amino acids only
        delete $AAonlyStop{$_} foreach keys %AAnonStop;  ## final stop codon names only
        my %AAchar = map {($_=>1)} values %AAnonStop;
        my %crossover;
        foreach (keys %{ $PARAMS{stopcodons} }) {
            $crossover{$_} = 1 if exists $AAchar{ $AAonlyStop{$_} };  ## we should NOT have any crossover between stop codon names and amino acid residues
        }
        if ($PARAMS{keepcase}) {
            my %AAcharLC = map { "\L$_" } keys %AAchar;
            foreach (keys %{ $PARAMS{stopcodons} }) {
                $crossover{$_} = 1 if exists $AAcharLC{ $AAonlyStop{$_} };  ## (might use lowercase AA characters) + (custom stops could be lowercase) = (ensure no overlap)
            }
        }
        if (%crossover) {  # bad choice of stop codon names
            my $fail = join ', ', sort keys %crossover;
            warn "With the given settings, the following stop characters overlap with legal AA residue characters: '$fail'.\n Switching stop characters to { # $ % }...\n";
            $PARAMS{aminoacids}{$_} = $PARAMS{stopcodons}{$_} foreach keys %{ $PARAMS{stopcodons} };
        }
    }
    my %stopchars = map { $PARAMS{aminoacids}{$_} } keys %{ $PARAMS{stopcodons} };
    $TRANS{stopcodons} = $PARAMS{stopcodons};
    $PARAMS{stopclass} = '[' . (join '', keys %stopchars) . ']';   # regex
    
    ## Determine $modeN given 'mode', 'frames', or 'from_frame' arguments; vet frames if 'frames' || 'from_frame'
    my (@framelist, $modeN, @useframes, %useframes, $useframes);
    my @knownframes = ('+0','+1','+2','-0','-1','-2');
    my %knownframes = map {($_=>1)} @knownframes;
    my $knownframes = join " ", @knownframes;
    my $modeN = ($PARAMS{mode} =~ /^bestof([36])$/) ? $1 : $PARAMS{mode};
    if ($modeN eq 'FRAMES') {
        @useframes = @{ $PARAMS{frames} };
        %useframes = map {($_=>1)} @useframes;
        $useframes = join " ", @useframes;
        my $failframe = 0;
        foreach my $frame (@useframes) {
            $failframe = 1 unless exists $knownframes{$frame};
        }
        die "\nIllegal frame names specified!\n Given: $useframes\n Allowed: $knownframes\n\n" if $failframe;
        if (exists $useframes{'-0'} || exists $useframes{'-1'} || exists $useframes{'-2'}) {
            $modeN = 6;
        } elsif (exists $useframes{'+1'} || exists $useframes{'+2'}) {
            $modeN = 3;
        } else {
            $modeN = 1;
        }
    }
    
    my @bases = split //, $PARAMS{sequence};
    my $L = scalar @bases;
    
    # Start of translatable portion of + frames; 0-based
    $PARAMS{pCodonStart}{$_} = $_ foreach (0..2);
    # End of translatable portion of - frames; 0-based
    $PARAMS{nCodonEnd}{$_} = $L - $_ - 1 foreach (0..2);
    
    if ($PARAMS{autotrim}) {
        # End of translatable portion of + frames; 0-based
        $PARAMS{pCodonEnd}{$_} = $_ + 3*int( ($L-$_)/3 ) - 1 foreach (0..2);
        # Start of translatable portion of - frames; 0-based
        $PARAMS{nCodonStart}{$_} = $L - $_ - 3*int( ($L-$_)/3 ) foreach (0..2);
    } else {
        # End of translatable portion of + frames; 0-based
        $PARAMS{pCodonEnd}{$_} = $#bases foreach (0..2);
        # Start of translatable portion of - frames; 0-based
        $PARAMS{nCodonStart}{$_} = 0 foreach (0..2);
    }
    
    if ($PARAMS{from_frame}) {
        
        if ($PARAMS{from_frame} eq '+0') {
            $PARAMS{frameseq}{'+0'} = join '', @bases[0..$PARAMS{pCodonEnd}{0}];
        } elsif ($PARAMS{from_frame} eq '+1') {
            $PARAMS{frameseq}{'+1'} = join '', @bases[1..$PARAMS{pCodonEnd}{1}];
        } elsif ($PARAMS{from_frame} eq '+2') {
            $PARAMS{frameseq}{'+2'} = join '', @bases[2..$PARAMS{pCodonEnd}{2}];
        } elsif ($PARAMS{from_frame} eq '-0') {
            my $temp = join '', @bases[$PARAMS{nCodonStart}{0}..$PARAMS{nCodonEnd}{0}];
            $PARAMS{frameseq}{'-0'} = ${ &revcomp(\$temp) };
        } elsif ($PARAMS{from_frame} eq '-1') {
            my $temp = join '', @bases[$PARAMS{nCodonStart}{1}..$PARAMS{nCodonEnd}{1}];
            $PARAMS{frameseq}{'-1'} = ${ &revcomp(\$temp) };
        } elsif ($PARAMS{from_frame} eq '-2') {
            my $temp = join '', @bases[$PARAMS{nCodonStart}{2}..$PARAMS{nCodonEnd}{2}];
            $PARAMS{frameseq}{'-2'} = ${ &revcomp(\$temp) };
        }
        
    } else {
        
        $PARAMS{frameseq}{'+0'} = join '', @bases[0..$PARAMS{pCodonEnd}{0}];
        
        if ($modeN == 1) {
            
            @framelist = ('+0');
            
        } else {
            
            $PARAMS{frameseq}{'+1'} = join '', @bases[1..$PARAMS{pCodonEnd}{1}];
            $PARAMS{frameseq}{'+2'} = join '', @bases[2..$PARAMS{pCodonEnd}{2}];
            
            if ($modeN == 3) {
                @framelist = ('+0','+1','+2');
            } elsif ($modeN == 6) {
                
                @framelist = ('+0','+1','+2','-0','-1','-2');
                my $temp0 = join '', @bases[$PARAMS{nCodonStart}{0}..$PARAMS{nCodonEnd}{0}];
                my $temp1 = join '', @bases[$PARAMS{nCodonStart}{1}..$PARAMS{nCodonEnd}{1}];
                my $temp2 = join '', @bases[$PARAMS{nCodonStart}{2}..$PARAMS{nCodonEnd}{2}];
                $PARAMS{frameseq}{'-0'} = ${ &revcomp(\$temp0) };
                $PARAMS{frameseq}{'-1'} = ${ &revcomp(\$temp1) };
                $PARAMS{frameseq}{'-2'} = ${ &revcomp(\$temp2) };
            }
            
        }
    }
    
    if ($PARAMS{mode} =~ /^bestof([36])$/) {
        
        my (%TEMP, %ORF);
        my $frame = '+0';
        my $stats = $PARAMS{stats};  # keep original copy
        $PARAMS{stats} = 1;   # stats must be on 
        
        @framelist = ('+0','+1','+2');
        push @framelist, ('-0','-1','-2') if $modeN == 6;
        $TEMP{$_} = &subtranslate($_, \%PARAMS) foreach @framelist;
        
        $ORF{ $TEMP{$_}{MAX_ORF_AA} }{$_} = 1 foreach @framelist;
        my $max_orf = (sort {$b <=> $a} keys %ORF)[0];
        foreach my $frame (keys %{ $ORF{$max_orf} }) {
            if (!$stats) {  # stats were originally off
                delete $TEMP{$frame}{$_} foreach qw/ LENGTH_NT LENGTH_AA N_STOPS N_NON_ACGT N_CODON_FAILS MAX_ORF_AA MAX_ORF_NT MAX_ORF_COUNT /;
            }
            $TRANS{$frame} = dclone($TEMP{$frame});
        }
        
    } else {
        
        @framelist = @{ $PARAMS{frames} } if $PARAMS{mode} eq 'FRAMES';
        $TRANS{$_} = &subtranslate($_, \%PARAMS) foreach @framelist;
        
    }
    
    return \%TRANS;
}





sub subtranslate {
    
    my ($frame, $PARAMS) = @_;    
    my %FRAME;
    
    my (@translation, @codon_table, $stops, $untranslatable, $unknown);
    my $seqlen = length($$PARAMS{frameseq}{$frame});
    $unknown += length($_) foreach split /[ACGT]+/i, $$PARAMS{frameseq}{$frame};
    my $startmax = int($seqlen/3);
    
    foreach my $start (0..$startmax) {      # one for each codon
        $start *= 3;
        if ($start == $seqlen) {
            # done
            $FRAME{MESSAGE} = 'Sequence was correct length for full translation';
        } elsif (abs($seqlen-$start) <= 2) {    # reached end of translatable sequence
            my $off = abs($seqlen-$start) == 2 ? 1 : 2;
            my $msg = "Sequence was $off base";
            $msg .= "s" if ($off > 1);
            $msg .= " too short for full translation";
            $FRAME{MESSAGE} = $msg;
        } else {
            my ($codon, $residue);
            if ($$PARAMS{keepcase}) {
                my $Rcodon = substr($$PARAMS{frameseq}{$frame}, $start, 3);
                my $Ucodon = $codon = "\U$Rcodon";
                $residue = $$PARAMS{aminoacids}{$Ucodon};
                $stops++ if $$PARAMS{stopcodons}{$Ucodon};
                $residue = "\L$residue" if $Rcodon ne $Ucodon;  # raw had lowercase letters
            } else {
                $codon = substr($$PARAMS{frameseq}{$frame}, $start, 3);
                $residue = $$PARAMS{aminoacids}{$codon};
                $stops++ if $$PARAMS{stopcodons}{$codon};
            }
            if ($residue) {  # don't keep codons for untranslatable material
                push @translation, $residue;
                push @codon_table, "$codon\t$residue";
            } else {        # codon failed to translate: probably due to Ns or use of extended IUPAC alphabet
                #        push @errors, "header $HEADER | residue $residue | codon $codon | start $start | frame $frame | len $seqlen | seq $$PARAMS{frameseq}{$i}\n";
                #        print "FAIL: '$codon' = '$residue'\n";
                push @translation, 'X';
                push @codon_table, "$codon\tX";
                $untranslatable++;
            }
        }
    }
    
    if ($$PARAMS{ascodons}) {
        
        $FRAME{CODON_TABLE} = \@codon_table;
        
    } else {
        
        my $AASEQ = join '', @translation;
        $FRAME{SEQUENCE_NT} = $$PARAMS{frameseq}{$frame};
        $FRAME{SEQUENCE_AA} = $AASEQ;
        my ($orfmax, $maxorfcount, %orflens, @orftest);
        my $stopclass = $$PARAMS{stopclass};
        my @orftest = split /$stopclass+/, $AASEQ;    # split on stops
        my (%orflens, $maxorfcount);
        $orflens{ length($orftest[$_]) }{$orftest[$_]} = $_ foreach (0..$#orftest);
        my $orfmax = (sort {$b <=> $a} keys %orflens)[0];
        my $maxorfrank = $#orftest;

        if ($$PARAMS{stats}) {
            $FRAME{LENGTH_NT} = $seqlen;
            $FRAME{LENGTH_AA} = $startmax;
            $FRAME{N_STOPS} = $stops || 0;
            $FRAME{N_NON_ACGT} = $unknown || 0;
            $FRAME{N_CODON_FAILS} = $untranslatable || 0;
            $FRAME{MAX_ORF_AA} = $orfmax;
            $FRAME{MAX_ORF_NT} = $orfmax*3;
            $FRAME{MAX_ORF_COUNT} = scalar keys %{ $orflens{$orfmax} };
        }

        if ($$PARAMS{orfs}) {
            while (my ($orfseq, $orfrank) = each %{ $orflens{$orfmax} }) {    # allows for > 1 longest ORF
                $maxorfcount++;
                ##         $FRAME{ORFS}{$maxorfcount}{ORF_SEQUENCE_NT} = ;   # eventually
                $FRAME{ORFS}{$maxorfcount}{ORF_SEQUENCE_AA} = $orfseq;
                if ($$PARAMS{stats}) {
                    $FRAME{ORFS}{$maxorfcount}{ORF_LENGTH_AA} = $FRAME{MAX_ORF_AA};
                    $FRAME{ORFS}{$maxorfcount}{ORF_LENGTH_NT} = $FRAME{MAX_ORF_NT};
                    $FRAME{ORFS}{$maxorfcount}{ORF_IS_5P_TERMINAL} = $orfseq eq $orftest[0] ? 1 : 0;     # ORF is 5'-terminal?
                    $FRAME{ORFS}{$maxorfcount}{ORF_IS_3P_TERMINAL} = $orfseq eq $orftest[-1] ? 1 : 0;    # ORF is 3'-terminal?
                    $FRAME{ORFS}{$maxorfcount}{ORF_STARTS_WITH_M} = $orfseq =~ /^M/i ? 1 : 0;            # ORF starts with an M?
                    ## this last is a little complex to determine after we've split on stop codons
                    $FRAME{ORFS}{$maxorfcount}{ORF_ENDS_WITH_STOP} = ($orfrank < $maxorfrank) ? 1 : ($AASEQ =~ /$stopclass$/i) ? 1 : 0;  # ORF ends with a stop?
                }
            }
        }
        
    }
    
    return \%FRAME;
}





sub transpose {
    
    ## takes a reference to a 2D matrix, which must be stored as an array-of-arrays
    ## returns a reference to the transposed matrix (still A-O-A)
    
    my @matrix = @{ $_[0] };
    my @trans;
    
    my $rows = $#matrix;
    my $cols = $#{ $matrix[0] };
    
    foreach my $i (0..$cols) {
        my @newrow;
        foreach my $j (0..$rows) {
            push @newrow, $matrix[$j][$i];
        }
        push @trans, \@newrow;
    }
    
    my $trows = $#trans;
    my $tcols = $#{ $trans[0] };
    
    if ($cols == $trows && $rows == $tcols) {
        return \@trans;
    } else {
        print "Matrix transposition failed! ($rows,$cols) => ($trows,$tcols)\n";
        return [];  # empty array ref
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



#sub collapse_memory {
#    ## FUTURE FUNCTION
#    ## takes a memory bytes value and converts to string like '10G'
#    ## Issue: don't want '10.5G', want '1050M', logic is tedious, not writing it now...
#    
#    my $NUM = shift;
#    my $FINAL;
#
#    die "$0: memory value must be >= 1!\n" if $NUM < 1;
#    if ($NUM > 1024^4) {
#        $FINAL = sprintf("$NUM/1024^4
#    } elsif ($NUM > 1024^3) {
#    } elsif ($NUM > 1024^2) {
#    } elsif ($NUM > 1024^1) {
#    } elsif ($NUM > 1024^0) {
#    } else {
#    }
#}





sub validate_memory {
    
    my ($ASK, $FIX) = @_;
    my $FINAL;
    
    my ($NUM, $ORD) = ($ASK =~ /(\d+)([TtGgMmKkBb]?)/);
    
    ## Test memory order-of-magnitude character, if present
    my %ORDS = ('T',1024^4, 'G',1024^3, 'M',1024^2, 'K',1024^1, 'B',1024^0);
    if ($ORD && !exists $ORDS{"\U$ORD"}) {
        die "$0: Memory string '$ASK' order-of-magnitude character '$ORD' is invalid!  Must match [TtGgMmKkBb].\n";
    }
    
    ## Test numeric component
    if (looks_like_number($NUM)) {
        die "$0: Numeric component of '$ASK' must be > 1!\n" if $NUM < 1;
        my $FULL = $ORD ? $NUM * $ORDS{"\U$ORD"} : $NUM;
        chomp(my $SNAME = `hostname`);
        chomp(my ($SRAM, $SORD) = split / /, `head -1 /proc/meminfo | awk '{ print \$2" "\$3 }'`);
        $SORD =~ s/B$//;
        my $SFULL = $SRAM * $ORDS{"\U$SORD"};
        if ($SFULL < $FULL) {
            my $msg = "Memory string '$ASK' exceeds server '$SNAME' available memory of $SFULL bytes!";
            if ($FIX eq 'die') {
                die "$0: $msg\n";
            } elsif ($FIX eq 'force') {
                $FINAL = $ASK;
                print STDERR "WARNING: $msg\nContinuing anyway...\n";
            } elsif ($FIX eq 'downsize') {
                $FINAL = int(0.75*$SFULL);
                print STDERR "WARNING: $msg\nDownsizing to $FINAL\n";
            } else {
                die "$0: fix method '$FIX' must be either 'force', 'downsize', or 'die'!\n";
            }
        }
    } else {
        die "$0: Memory string '$ASK' does not contain a recognizable numeric component!\n";
    }
    
    return $FINAL;
}





sub vcf_orderheader {
    
    ## Orders the header of a VCF file; input is a reference to an array of strings, 
    ##   which holds ONLY header data, each line being one array element.
    ## Can take a genome argument $GENO, in which case:
    ##   Adds/replaces contig lines in the header; also adds reference line if it doesn't already exist.
    ##   The file /n/data1/genomes/indexes/$GENO/$GENO.chrom.sizes must exist.
    
    my ($VCF, $GENO) = @_;
    
    my (@contigs, $greference);
    if ($GENO) {
        my $GPREF = "$indexes/$GENO/$GENO";
        my $chrsz = "$GPREF.chrom.sizes";
        my $CS = &open2('R', $chrsz, 'required chrom.sizes file');
        while (<$CS>) {
            chomp;
            next unless $_;
            my ($chr, $len) = split /\t/, $_;
            push @contigs, "##contig=<ID=$chr,length=$len>";
        }
        close $CS;
        $greference = "##reference=file://$GPREF.fa";
    }
    
    if ($VCF =~ /^ARRAY\(/) {
        my ($first, @known, @reference, @unknown, $last);
        foreach (@{$VCF}) {
            s/[\n\r]+$//;
            next unless $_;
            if (/^##fileformat/) {
                $first = $_;
            } elsif (/^#CHROM\s/) {
                $last = $_;
            } elsif (/^##(INFO|FILTER|FORMAT|ALT|SAMPLE|PEDIGREE)=/) {
                push @known, $_;
            } elsif (/^##reference=/) {
                push @reference, $_;
            } elsif (/^##contig=/) {
                push @contigs, $_ unless $GENO;
                ## if $GENO, do not retain contig lines; will replace any extant contig lines with guaranteed-correct ones in correct order.
            } else {
                push @unknown, $_;
            }
        }
        my @FINAL;
        foreach ($first, (sort @unknown), (sort @known), @contigs, $greference, @reference, $last) {  # @contigs: inject new, or reposition existing
            push @FINAL, "$_\n" if $_;
        }
        return \@FINAL;
    } else {
        die "$0: &vcf_orderheader was not given an array ref!\n";
    }
    
}





sub read_org_metadata {
    
    ## Reads the _ORGANISM_METADATA file
    ## Returns a hash with keys = taxon IDs, subkeys = columns names and values = cell values
    ## OPTIONAL: given a txon ID, returns ONLY the matching subhash (if any matching)
    
    my $TAXON = shift;
    my $ORGMETA = "$buildroot/code/data/_ORGANISM_METADATA";
    my %HASH;
    
    open my $IN, '<', $ORGMETA or die "$0: could not open required file '$ORGMETA': $!\n";
    my ($i, @colnames);
    while (<$IN>) {
        chomp;
        next if /^#/;
        my @data = split /\t/, $_;
        if (++$i == 1) {
            @colnames = @data;
        } else {
            $HASH{ $data[0] }{ $colnames[$_] } = $data[$_] foreach (0..$#data);  # $data[0] is the taxon ID
        }
    }
    close $IN;
    
    if ($TAXON) {
        if (exists $HASH{$TAXON}) {
            return \%{ $HASH{$TAXON} };
        } else {
            print "No records found for taxon ID '$TAXON'!\n";
        }
    } else {
        return \%HASH;
    }
}





sub read_VARIABLES {
    
    ## Reads a genome/transcriptome prep _VARIABLES/_TVARIABLES file
    ## Returns a hash with key=>{value, label, comment, line num} structure
    
    my ($GENO, $ANNO) = @_;
    my (%HASH, $VARS, $gatype);
    if ($ANNO) {
        $VARS = "$buildroot/preps/$GENO/$ANNO/_TVARIABLES";
        $gatype = 'transcriptome';
    } else {
        $VARS = "$buildroot/preps/$GENO/_VARIABLES";
        $gatype = 'genome';
    }
    
    open my $IN, '<', $VARS or die "$0: cannot read $gatype variables file '$VARS': $!\n";
    my $i;
    while (<$IN>) {
        chomp;
        next if /^#/;
        my ($keyval, $desc) = split /#/, $_, 2;
        my ($key, $val) = split /=/, $keyval, 2;
        my ($label, $comment) = ($desc =~ /\s+\[(.*?)\]\s+(.*)/);
        $comment =~ s/\s*$//;  # strip trailing whitespace
        $val =~ s/\s*$//;  # strip trailing whitespace
        $val =~ s/"//g;    # strip quotes
        $HASH{$key} = {'VALUE'=>$val, 'LABEL'=>$label, 'COMMENT'=>$comment, 'LINE'=>++$i};
    }
    close $IN;
    
    ## ADD NEW KEY: home-base URL (if possible)
    if (exists $HASH{home_db} && $HASH{home_db}{VALUE}) {
        chomp(my $url = `grep $HASH{home_db}{VALUE}_url $buildroot/code/data/_DATABASES_AND_URLS`);
        $url =~ s/.*=//;
        $url =~ s/ .*//;
        $HASH{home_url}{VALUE} = $url;
    }
    
    ## expanded (readme-friendly) ribosome subunits field
    if (exists $HASH{ribounits} && $HASH{ribounits}{VALUE}) {
        (my $ribo_in = $HASH{ribounits}{VALUE}) =~ s/"//g;
        my %nameconv = ('L','Large', 'S','Small', 'O','', 'ML','Mito Large', 'MS','Mito Small');
        my @ribo_out;
        foreach my $group (split /;/, $ribo_in) {
            my ($name, $value) = split /:/, $group;
            foreach my $size (split /,/, $value) {
                my $string = "${size}S";
                $string .= " ($nameconv{$name})" if $nameconv{$name};
                push @ribo_out, $string;
            }
        }
        $HASH{RIBOUNITS_FULL} = {'VALUE'=>join(', ',@ribo_out), 'LABEL'=>'Ribosome Subunits', 'COMMENT'=>'Pretty-formatted ribosome subunit weights', 'LINE'=>$HASH{ribounits}{LINE}};
    }
    
    return \%HASH;
}

   
sub canonicalize {
    
    ## canonicalize a path
    ## basically 'readlink' that doesn't follow links
    
    my $path = shift;
    
    ## Correct start of path
    my $tmp;
    if ($path =~ /^\//) {
        ## Then ok
        $tmp = $path;
    } elsif ($path =~ /^~/) {
        (my $path2 = $path) =~ s/^~//;
        $tmp = "$tilde/$path2";
    } else {
        chomp(my $wdir = `pwd`);
        $tmp = "$wdir/$path";
    }
    
    ## Expand any '..' segments
    if ($tmp =~ /\.\./) {
        my ($path1) = ($tmp =~ /^(.*\.\.)\/[^.]/);
        (my $path2 = $tmp) =~ s/$path1\///;
        chomp(my $path1f = `readlink -f "$path1"`);
        return "$path1f/$path2";
    } else {
        return $tmp;
    }
}

   
sub extract_names_from_paths {
    
    my @files = @{ $_[0] };
    
    sub crop_longest_extension {
        my ($x, $dir) = @_;
        die "$0: crop_longest_extension 'dir' must be \"left\" or \"right\"!\n" if ($dir ne 'left' && $dir ne 'right');
        my @y = map { [split(/\./,$_)] } @$x;
        if ($dir eq 'left') {
            $_ = [reverse(@$_)] foreach @y;
        }
        my $continue = 1;
        while ($continue) {
            my %firsts = map {($_=>1)} (map { $_->[0] } @y);
            if (scalar(keys %firsts)==1) {
                shift @$_ foreach @y;  # if first words all same, delete
            } else {
                $continue = 0;
            }
        }
        if ($dir eq 'left') {
            $_ = [reverse(@$_)] foreach @y;
        }
        $_ = join('.', @$_) foreach @y;
        return \@y;
    }
    
    sub soft_root_path {
        my $file = shift;
        if ($file =~ /^~/) {
            $file = s/^~/$tilde/;
        } elsif ($file =~ /\.\./) {
            (my $postdots = $file) =~ s!.*../!!;
            if ($postdots eq $file) {
                ## then do nothing
            } else {
                (my $path = $file) =~ s/$postdots//;
                chomp($path = `readlink -f $path`);
                $file = "$path/$postdots";
            }
        }
        return $file;
    }
    
    my (@paths, @paths2, $ml);
    foreach (@files) {
        die "$0: some input files to 'extract_names_from_paths' do not exist!\n" unless -e $_;
        $_ = &soft_root_path($_);
        my @parts = reverse(split(/\//, $_));
        push @paths, \@parts;
        $ml = scalar(@parts) if scalar(@parts) > $ml;
    }
    my $P = scalar @paths;
    foreach (@paths) {
        my @blank = map {""} (1..$ml);
        @blank[0..$#{$_}] = @$_;
        push @paths2, \@blank;
    }
    my @var = map {''} (1..$ml);

    my (@w1, @w2, $v1, $v2);
    foreach my $i (1..$ml) {
        my %pending = map {($_=>1)} (map { $_->[$i-1] } @paths2);
        if (scalar(keys %pending)==1) {
            ## no variation at this level in the file paths
            $var[$i-1] = 0;
        } elsif (scalar(keys %pending)==$P) {
            ## this level can fully distinguish all files
            $var[$i-1] = 2;
            $v2++;
            push @w2, $i-1;
        } else {
            ## this level can partially distinguish files; retain just in case
            $var[$i-1] = 1;
            $v1++;
            push @w1, $i-1;
        }
    }
    
    my @final;
    if ($v2) {
        
        ## this is all we need
        @final = map { $_->[$w2[0]] } @paths2;  # just take single right-most, if > 1 (vectors are still reversed)
        if ($w2[0]==0) {
            @final = @{ &crop_longest_extension(\@final,"left") };   ## w2 set included the filename itself; crop extensions
            @final = @{ &crop_longest_extension(\@final,"right") };
        }
        
    } elsif ($v1) {
        
        ## use all of these; will test if this is sufficient
        foreach my $ref (@paths2) {
            push @final, join("/", reverse(@$ref[@w1]));
        }
        my %test = map {($_=>1)} @final;
        
        if (scalar(keys %test)==$P) {
            ## sufficient to distinguish all files
            if ($w1[-1]==$P) {
                @final = @{ &crop_longest_extension(\@final,"left") };   ## w1 set included the filename itself; crop extensions
                @final = @{ &crop_longest_extension(\@final,"right") };
            }
        }
        
    }
    return \@final;
    
}

1;
