#!/bin/bash

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate alevinfry

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823

pushd .
    
SALMONDIR=salmon_spliceu
PERMITDIR=${SALMONDIR}/out_permit_knee
OUTDIR=${SALMONDIR}/quant_permit_knee
INDEX=/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/

#alevin-fry quant -m ${INDEX}spliceu_t2g_3col.tsv -i $PERMITDIR -o $OUTDIR -t 16 -r cr-like-em --use-mtx --dump-eqclasses
alevin-fry quant -m ${INDEX}t2t.tsv -i $PERMITDIR -o $OUTDIR -t 16 -r cr-like-em --use-mtx --dump-eqclasses

popd

# Command to sort by CB > UMI > read. 
# zcat map.txt.gz | grep "DIR:true" | sort -t$'\t' -k4,4 -k5,5 -k1,1 -k2,2 | gzip >  sorted.txt.gz