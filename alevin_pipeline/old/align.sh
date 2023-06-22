#!/bin/bash

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/SRR

INDEX=/gpfs/commons/groups/knowles_lab/index/salmon/mm10_cdna

pushd .

while read SRR; do
    echo ${SRR}
    salmon alevin -i $INDEX -l A -1 ${SRR}_1.fastq -2 ${SRR}_corrected_2.fastq -p 32 --splitseqV1 -o ../salmon_output/${SRR}_salmon --sketch
done < SRR_Acc_List.txt

popd