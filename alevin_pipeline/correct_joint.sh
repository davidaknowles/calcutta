#!/bin/bash

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate alevinfry

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823

pushd .
    
SALMONDIR=salmon_GRCm38
PERMITDIR=${SALMONDIR}/out_permit_knee

alevin-fry generate-permit-list -d both -i ${SALMONDIR} --output-dir $PERMITDIR -k
alevin-fry collate -r ${SALMONDIR} -t 16 -i $PERMITDIR
# alevin-fry quant -m splici_index_reference/transcriptome_splici_fl61_t2g_3col.tsv -i ./SRR6750057_out_permit_knee -o ./SRR6750057_counts -t 16 -r cr-like-em --use-mtx. # Don't have t2g mapping :(

popd

# Command to sort by CB > UMI > read. 
# zcat map.txt.gz | grep "DIR:true" | sort -t$'\t' -k4,4 -k5,5 -k1,1 -k2,2 | gzip >  sorted.txt.gz