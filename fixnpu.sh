#!/bin/bash

exit;

npu=$1
np=${npu%.underflow}

zero=$(grep "@ 0 peaks" $npu | wc -l)
if [ $zero -eq 1 ]; then
    peaks=$(cat $np | wc -l)
    if [ $peaks -eq 100000 ]; then
        rm -f $npu
    else
        p=$(grep -oP "-p \S+" $npu | cut -f2 -d' ')
        pct1=$(echo "scale=3; 100*$peaks/100000" | bc)
        pct2=$(sprintf "%0.2f" $pct1)
        echo "Underflow -p $p @ $peaks peaks ($pct2% of 100000 quota)." > $npu
    fi
fi

