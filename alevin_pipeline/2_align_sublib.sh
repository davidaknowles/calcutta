#!/bin/bash

#SBATCH --job-name=align_microglialess
#SBATCH --cpus-per-task=18
#SBATCH --mem=164G
#SBATCH --time=4-00:00:00

# command arg 1: name of sublib
# arg 2: first SRR
# arg 3: second SRR

pushd .

cd /gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/SRA

#INDEX=/gpfs/commons/groups/knowles_lab/index/salmon/mm10_cdna
#INDEX=/gpfs/commons/groups/knowles_lab/index/salmon/Mus_musculus.GRCm38.cdna
INDEX=/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/index_dir

# need V2 for this data: actually used Parse kit? 
salmon alevin -i $INDEX -l A -1 $2_1.fastq.gz $3_1.fastq.gz -2 $2_2.fastq.gz $3_2.fastq.gz -p 20 --splitseqV2 -o ../salmon_per_sublib/$1 --sketch

popd