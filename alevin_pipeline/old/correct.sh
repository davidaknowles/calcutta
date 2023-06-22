#!/bin/bash

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate alevinfry

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823

INDEX=/gpfs/commons/groups/knowles_lab/index/salmon/mm10_cdna

pushd .

while read SRR; do
    echo ${SRR}
    
    SALMONDIR=salmon_output/${SRR}_salmon
    PERMITDIR=${SALMONDIR}/out_permit_knee

    alevin-fry generate-permit-list -d both -i ${SALMONDIR} --output-dir $PERMITDIR -k
    alevin-fry collate -r ${SALMONDIR} -t 16 -i $PERMITDIR
    alevin-fry quant -m splici_index_reference/transcriptome_splici_fl61_t2g_3col.tsv -i ./SRR6750057_out_permit_knee -o ./SRR6750057_counts -t 16 -r cr-like-em --use-mtx # Don't have t2g mapping :(

done < SRR/SRR_Acc_List.txt

popd