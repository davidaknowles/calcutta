#!/bin/bash

tsv_file=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/sublib_SRR_map.tsv

# Loop over the lines in the TSV file
while IFS=$'\t' read -r sublib srr1 srr2; do
    echo Sublib: $sublib SRR1: $srr1 SRR2: $srr2
    sbatch 3_correct_and_quant_param.sh $sublib
done < "$tsv_file"