#!/bin/bash

module load bustools/0.39.3

#INDEX=~/knowles_lab/index/kallisto/mus_musculus/transcriptome.idx

pushd .

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/kallisto_output/merged


bustools count sorted.bus -o count.txt -g ~/knowles_lab/index/kallisto/mus_musculus/transcripts_to_genes.txt -e matrix.ec -t transcripts.txt -m 

bustools count correct100.bus -o genecounts -g ~/knowles_lab/index/kallisto/mus_musculus/transcripts_to_genes.txt -e matrix.ec -t transcripts.txt -m --genecounts

popd
