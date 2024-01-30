#!/bin/bash

module load sratoolkit/3.0.0

set -e

pushd .

cd /gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/SRA_with_tech

while read SRR; do
    echo $SRR
    fasterq-dump $SRR --include-technical
    break
done < SRR_Acc_List.txt

popd