#!/bin/bash

module load bustools/0.39.3

#INDEX=~/knowles_lab/index/kallisto/mus_musculus/transcriptome.idx

pushd .

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/kallisto_output/merged

#which bustools
# bustools merge -o merged2 SRR67500*_kb_out # this worked with bustools 0.39
# hmm guess i don't know if the output of this is sorted? 

# this failed with bustools 0.42: 
# bustools merge -o merged2 -t SRR6750043_kb_out/transcripts.txt -e SRR6750043_kb_out/matrix.ec SRR67500*_kb_out/sorted.bus  
# bus records read:    14113825
# bus records written: 0
# bus records lost:    13464519

bustools count sorted.bus -o count.txt -g ~/knowles_lab/index/kallisto/mus_musculus/transcripts_to_genes.txt -e matrix.ec -t transcripts.txt -m 

bustools count correct100.bus -o genecounts -g ~/knowles_lab/index/kallisto/mus_musculus/transcripts_to_genes.txt -e matrix.ec -t transcripts.txt -m --genecounts

popd
