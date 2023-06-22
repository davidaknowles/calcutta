#!/bin/bash

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate alevinfry

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823

pushd .
    
SALMONDIR=salmon_spliceu
POSSIBLE_BARCODES=~/calcutta/calcutta/splitseqv1_meta/splitseqv1_bc_all_combo.txt
OUTDIR=${SALMONDIR}/quant_t2t
INDEX=/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/
USE_KNOW_BARCODES=true

if [ "$USE_KNOW_BARCODES" = true ] ; then
    PERMITDIR=${SALMONDIR}/out_permit_known
    alevin-fry generate-permit-list -d both -i ${SALMONDIR} --output-dir $PERMITDIR --unfiltered-pl $POSSIBLE_BARCODES --min-reads 100
else
    PERMITDIR=${SALMONDIR}/out_permit_knee
    alevin-fry generate-permit-list -d both -i ${SALMONDIR} --output-dir $PERMITDIR -k  # without known list
fi

alevin-fry collate -r ${SALMONDIR} -t 16 -i $PERMITDIR

# use spliceu_t2g.tsv instead of t2t.tsv for gene (rather than isoform) level quantification
# OR spliceu_t2g_3col.tsv for USA mode
alevin-fry quant -m ${INDEX}t2t.tsv -i $PERMITDIR -o $OUTDIR -t 16 -r cr-like-em --use-mtx --dump-eqclasses

popd

# Command to sort by CB > UMI > read. 
# zcat map.txt.gz | grep "DIR:true" | sort -t$'\t' -k4,4 -k5,5 -k1,1 -k2,2 | gzip >  sorted.txt.gz