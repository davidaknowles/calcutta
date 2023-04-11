#!/bin/bash

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/SRR

#INDEX=/gpfs/commons/groups/knowles_lab/index/salmon/mm10_cdna
INDEX=/gpfs/commons/groups/knowles_lab/index/salmon/Mus_musculus.GRCm38.cdna

pushd .

salmon alevin -i $INDEX -l A -1 *_1.fastq -2 *_corrected_2.fastq -p 20 --splitseqV1 -o ../salmon_GRCm38 --sketch

popd