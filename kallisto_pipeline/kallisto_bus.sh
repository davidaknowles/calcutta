#!/bin/bash

INDEX=~/knowles_lab/index/kallisto/mus_musculus/transcriptome.idx
OUTDIR=../kallisto_splitp/

pushd .

cd /gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/SRR

while read SRR; do
    #kallisto bus -i $INDEX -x 1,10,18,1,48,56,1,78,86:1,0,10:0,0,0 -o ${OUTDIR}${SRR}_kb_out ${SRR}_1.fastq ${SRR}_corrected_2.fastq
    kallisto bus -i $INDEX -x SPLIT-SEQ -o ${OUTDIR}${SRR}_kb_out ${SRR}_1.fastq ${SRR}_corrected_2.fastq
    bustools sort -o ${OUTDIR}${SRR}_kb_out/sorted.bus ${OUTDIR}${SRR}_kb_out/output.bus
done < SRR_Acc_List.txt

popd

else if (opt.technology == "SPLIT-SEQ") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(0,0,0));
        busopt.umi.push_back(BUSOptionSubstr(1,0,10));
        busopt.bc.push_back(BUSOptionSubstr(1,10,18));
        busopt.bc.push_back(BUSOptionSubstr(1,48,56));
        busopt.bc.push_back(BUSOptionSubstr(1,78,86));
        -s 87 -e 94
        
