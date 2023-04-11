#!/bin/bash

module load bustools/0.39.3

#INDEX=~/knowles_lab/index/kallisto/mus_musculus/transcriptome.idx

pushd .

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/kallisto_output

#which bustools
bustools merge -o merged2 SRR67500*_kb_out # this worked with bustools 0.39

# bustools merge -o merged2 -t SRR6750043_kb_out/transcripts.txt -e SRR6750043_kb_out/matrix.ec SRR67500*_kb_out/sorted.bus  
# this failed with bustools 0.42: 
# bus records read:    14113825
# bus records written: 0
# bus records lost:    13464519

popd
