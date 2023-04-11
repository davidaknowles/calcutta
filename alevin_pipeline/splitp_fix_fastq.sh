#!/bin/bash

source "$HOME/.cargo/env"
cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/SRR

pushd .

while read SRR; do
    echo ${SRR}
    splitp -r ${SRR}_2.fastq -b ~/calcutta/alevin_splitseq/oligo_hex_bc_mapping.txt -s 87 -e 94 -o > ${SRR}_corrected_2.fastq
done < SRR_Acc_List.txt

popd