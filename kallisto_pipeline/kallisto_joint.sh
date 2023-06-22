#!/bin/bash

INDEX=~/knowles_lab/index/kallisto/mus_musculus/with_precursor/index.idx
OUTDIR=../kallisto_joint

pushd .

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/SRR

kallisto bus -i $INDEX -x SPLIT-SEQ -t 10 -o ${OUTDIR} *.fastq 

cd ${OUTDIR}

bustools sort -o sorted.bus output.bus

bustools text -p sorted.bus | gzip > sorted.txt.gz

bustools whitelist -f 100 -o whitelist100.txt sorted.bus

bustools correct -o corrected100.bus -w whitelist100.txt sorted.bus 
 
bustools sort -o corrected100_sorted.bus corrected100.bus # is this necessary? Binary files corrected100.bus and corrected100_sorted.bus differ according to diff, and this is fast anyway. 

bustools text -p corrected100_sorted.bus | gzip > corrected100_sorted.txt.gz

popd
