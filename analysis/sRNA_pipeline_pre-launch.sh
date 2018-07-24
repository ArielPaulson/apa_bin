#!/bin/bash

echo -e '#!/bin/bash\n' > launch.sh
for fq in ../data/fastq/*.gz; do samp=${fq##*/}; samp=${samp%.fastq.gz}; echo "./run.sh $samp" >> launch.sh; done
echo -e "\nls ../data/output/*/*/*.sh | awk '{ print $1\" &\" }' > all_RGLQ.sh" >> launch.sh
echo "./all_RGLQ.sh" >> launch.sh

