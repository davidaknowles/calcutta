#!/bin/bash

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate pyroe

# module load bedtools/2.29.0 # needs 2.30! 

pushd .

cd ~/knowles_lab/index/kallisto/mus_musculus/with_precursor
 
pyroe make-spliceu --dedup-seqs GRCm39.primary_assembly.genome.fa gencode.vM32.basic.annotation.gtf ~/knowles_lab/index/salmon/mus_spliceu/

salmon index -t spliceu.fa -i spliceu -p 16

popd