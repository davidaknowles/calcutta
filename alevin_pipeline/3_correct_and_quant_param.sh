#!/bin/bash

#SBATCH --job-name=align_microglialess
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00

SALMONDIR=salmon_per_sublib/$1

POSSIBLE_BARCODES_V1=~/calcutta/calcutta/splitseqv1_meta/splitseqv1_bc_all_combo.txt
POSSIBLE_BARCODES_V2=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/Parse_meta/Parse_expanded_barcodes.txt
POSSIBLE_BARCODES=$POSSIBLE_BARCODES_V2

INDEX=/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/

USE_KNOWN_BARCODES=true

MAPPING=spliceu_t2g

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate alevinfry

pushd .

#cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823
cd /gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice

if [ "$USE_KNOWN_BARCODES" = true ] ; then
    PERMITDIR=${SALMONDIR}/out_permit_known
    alevin-fry generate-permit-list -d both -i ${SALMONDIR} --output-dir $PERMITDIR --unfiltered-pl $POSSIBLE_BARCODES --min-reads 100
else
    PERMITDIR=${SALMONDIR}/out_permit_knee
    alevin-fry generate-permit-list -d both -i ${SALMONDIR} --output-dir $PERMITDIR -k  # without known list
fi

alevin-fry collate -r ${SALMONDIR} -t 20 -i $PERMITDIR

OUTDIR=${PERMITDIR}/quant_$MAPPING

# use spliceu_t2g.tsv instead of t2t.tsv for gene (rather than isoform) level quantification
# OR spliceu_t2g_3col.tsv for USA mode
alevin-fry quant -m ${INDEX}${MAPPING}.tsv -i $PERMITDIR -o $OUTDIR -t 20 -r cr-like-em --use-mtx --dump-eqclasses

popd

# Command to sort by CB > UMI > read. 
# zcat map.txt.gz | grep "DIR:true" | sort -t$'\t' -k4,4 -k5,5 -k1,1 -k2,2 | gzip >  sorted.txt.gz