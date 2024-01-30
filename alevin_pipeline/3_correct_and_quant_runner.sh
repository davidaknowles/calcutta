#!/bin/bash

SRADIR=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/SRA

for file in $SRADIR/*_1.fastq.gz; do
    # Extract the part matching the *
    filename=$(basename "$file")
    SRR="${filename%%_1.fastq.gz}"
    echo $SRR
    sbatch 3_correct_and_quant_param.sh $SRR
done