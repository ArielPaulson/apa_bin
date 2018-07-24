#!/bin/bash

# DO NOT set -e	in function libraries -- turns command-line environment deadly!



#############
## ALIASES ##
#############

alias ls="LC_ALL=C LC_COLLATE=C ls --color=auto"
alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
#alias mkdir='mkdir -p'
#alias vi='vim'

alias rmbak='rm -f *~'
alias la='ls -la'
alias latr='ls -latr'
alias R='R --no-save --no-restore'
alias calc="libreoffice --calc"
alias myrun="mysql -h mysql-dev -u apa -p"

alias showTaxa=/home/apa/local/bin/showTaxa
alias GO_Tools=/home/apa/local/bin/GO_Tools
alias scratchspace=/home/apa/local/bin/scratchspace
alias editsrc=/home/apa/local/bin/editsrc


###############
## FUNCTIONS ##
###############

function join_by { 
        
    ## array-joining function
    ## call: join_by $delim ${array[@]}
    ## from https://stackoverflow.com/questions/1527049/join-elements-of-an-array
    
    local IFS="$1"
    shift
    echo "$*"
}


function root_symlink { 
    
    ## convert a symlink from a relative path to an absolute path
    ## call: root_symlink $link
    
    ln -sf $(readlink -f $1) $1
}


function restring_symlink { 
    
    ## substitute strings in a symlink
    ## call: restring_symlink $file $oldstr $newstr
    
    file=$1
    export OLDSTR=$2
    export NEWSTR=$3
    
    oldsym=$(readlink -f $file)
    if [ -z "$oldsym" ]; then oldsym=$(readlink $file); fi
    newsym=$(echo $oldsym | perl -pe 's/$ENV{OLDSTR}/$ENV{NEWSTR}/')
    ln -sf $newsym $file
}


function canonicalize { 
    
    ## canonicalize a path
    ## basically 'readlink' that doesn't follow links
    
    path=$1
    
    ## Correct start of path
    if [[ "$path" =~ ^"/" ]]; then
        ## Then ok
        tmp="$path"
    elif [[ "$path" =~ ^"~" ]]; then
        path2=${path#~/}
        tilde=$(readlink -f ~)
        tmp="$tilde/$path2"
    else
        wdir=$(pwd)
        tmp="$wdir/$path"
    fi
    
    ## Expand any '..' segments
    if [[ "$tmp" =~ ".." ]]; then
        path1=$(echo "$tmp" | perl -pe 's!(.*[^.]/|^)([./]+)/([^.].*)!$1$2!')
        path2=$(echo "$tmp" | perl -pe 's!(.*[^.]/|^)([./]+)/([^.].*)!$3!')
        path1f=$(readlink -f "$path1")
        echo "$path1f/$path2"
    else
        echo "$path2"
    fi
}


function match_STARidx {
    
    ## given a read length, genome and annotation labels, set the variable "staridx" to the best available STAR index path (if any).
    ## returns nothing if no suitable STAR index found.
    ## for instance, call: match_STARidx 65 mm10 Ens_80 && echo $staridx
    
    readlen=$1
    geno=$2
    anno=$3
    
    idx=/n/data1/genomes/indexes/$geno
    if [ -z $anno ]; then starpref=$idx/STAR_; else starpref=$idx/$anno/STAR_; fi
    staresc=${starpref//\//\\\/}
    starlen=($(ls -d $starpref* | sed "s/$staresc//" | sed 's/bp//' | sort -nr))  # SORT DECREASING
    Nstarlen=${#starlen[@]}
    if [ $Nstarlen -eq 1 ]; then matchblurb="ONLY MATCH"; else matchblurb="BEST MATCH"; fi
    
    uselen=0
    staridx=''
    for bp in ${starlen[@]}
    do
        if [ $readlen -le $bp ]; then
            uselen=$bp   # smallest STAR length that is >= read length
        fi
    done
    
    echo "Found STAR indexes:"
    for bp in ${starlen[@]}
    do
        if [ $bp -eq $uselen ]; then
	          echo "$starpref${bp}bp   * $matchblurb for ${readlen}bp*"
	          staridx=$starpref${uselen}bp
	      else
	          echo "$starpref${bp}bp"
        fi
    done
    
    if [ $uselen -eq 0 ]; then
        echo -e "Failed to match read length '$readlen' to any STAR index lengths!"   # : ( ${starlen[@]} )!"
    fi
    echo ""
}


function check_stage () {
    
    ## Stage-completion test function.
    ## Designed for GATK pipeline, but can be used for anything.
    ## Requires $STAGE (last-run stage number) and $skip_fwd (1|0) to be already specified in the environment.
    ## $skip_fwd=1 means all prior stage output files exist; chain is not broken.  If chain is broken, $skip_fwd=0; must rerun all stages from this point.
    ## Sets value of $run_stage; may reset value of $skip_fwd.
    
    this_stage=$1     # an integer
    final_output=$2   # a filename
    
    if [ -e $final_output ] && [ $skip_fwd -eq 1 ]; then file_skip=1; else file_skip=0; fi    # does final output exist AND skip_fwd still true?
    if [ $STAGE -gt $this_stage ]; then stage_skip=1; else stage_skip=0; fi                   # is last-completed stage greater than pending stage?
    
    if [ $stage_skip -ne 1 ] || [ $file_skip -ne 1 ]; then run_stage=1; else run_stage=0; fi  # run pending stage or not?
    if [ $run_stage -eq 1 ]; then skip_fwd=0; fi                                              # if stage must be re-run, we can no longer skip_fwd
}


