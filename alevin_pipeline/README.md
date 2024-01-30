## Installing stuff

`pyroe` makes an index that includes introns (splici) or precursors (spliceu). 
libs are too old on the cluster to use conda install
pip install fails under python 3.11 because numba requires <3.11, but seems to work with python 3.10
docs: https://pyroe.readthedocs.io/en/latest/building_splici_index.html

`alevin fry` can be installed using conda. 

I compiled `salmon 1.10.0` from source into `/gpfs/commons/home/daknowles/knowles_lab/software/bin/`

## Isoform aware alevin fry pipeline

1_index.sh builds an index that include both mature mRNA and a precursor "transcript" for each gene, which alevin calls a "spliceu" index (as opposed to a splici index which has each intron represented separately). 
2_align.sh (pseudo)aligns reads creating a "RAD" file.
3_correct_and_quant.sh tries to correct cell barcodes either to a white list (this doesn't work well for some reason) or de novo, and then performs quantification which gives isoform or gene level abundance (expected counts), plus a equivalence class matrix that could be used for downstream analysis. 
4_quant.sh is just the quantification step of (3) standalone (in case you want to re-run just that.)

Confusing thing about average number of UMIs per cell: went down from ~1600 with quant(knee) to 287 with quant_known (min_reads==100). 
How many of the quant(knee) BCs look legit? 98%. 

splitp_fix_fastq.sh is a script to make the matching polyA and random hexamer barcodes match (from https://combine-lab.github.io/alevin-fry-tutorials/2022/split-seq/). I think it is probably preferable to NOT do this and just combine the counts for the half cells at the end. 

### Specifically for the "microglialess" dataset: 
The "srr" versions run each SRR fastq pair separately. "runner.sh" sets the sbatch jobs going. This was necessary to figure out the mapping between SRRs and sublibs via the preprocessed data. 
The "sublib" and "sublib_runner" versions run each sublib (consisting of two SRR IDs for some reason) separately. 

## Details

Main issue was rad files are big. 
Good news is the collated file is sorted by BC so can process with a lot less memory. 
Even better news is `alevin fry quant` will make a sparse cells by EC matrix via the --dump-eqclasses option. 

Initially used "mm10_cdna" index from
http://refgenomes.databio.org/v3/genomes/splash/fa159612d40b1bedea9a279eb24999b3d27145f9dd70dcca
This doesn't include precursor mRNA. 

## Notes on USA mode

Unspliced/Spliced/Ambiguous mode. 
This mode is automatic if the file specified by the -m option to `alevin-fry quant` has 3 columns: transcript, gene, S or U. 
e.g. spliceu_t2g_3col.tsv, ~140k entries.  
For U(nspliced) the transcript will be e.g. ENSMUSG00000095366.3-I, note G(ene) not T(ranscript)!
t2g.tsv (2 column!) has 140k entries trancript -> gene. 
Not sure what the point of g2g.tsv is. 

quants_mat_cols.txt ends up having 170k targets: genes (no suffix), unspliced (-U suffix) and ambiguous (-A suffix)
Presumably A is S or U (e.g. reads coming from exons)
These are the targets listed in the ECs. 
These are correct/collapsed UMI counts, NOT TPMs. 

So to get transcript level EC I can just make a dummy `t2g` mapping that is just `t2t`. 
Seems to work fine, see results in quant_t2t. 

# quant_t2t and --dump-eqclasses

- gene_eqclass.txt.gz: 
    - line 1: the number of transcripts (genes normally) (or the number of USA targets)
    - line 2: the number of ECs. 
    - The rest: one line per EC, where the first n-1 elements are transcript (normally gene) IDs, and the last entry is the EC idx (which is not the ordering in the file). 

- geqc_counts.mtx is cells x ECs, with the row labels (cell barcodes being given by quants_mat_rows.txt). The indexing for the columns is the EC idx (the last entry of each line in gene_eqclass.txt.gz). 

## Handling deduplicate transcripts

When pyroe builds the index it has the option to check for (and remove) duplicates in the generated transcriptome fasta. These are stored in duplicate_entries.tsv but are NOT removed from t2g.tsv. Since we include precursor sequences, all single exon genes are duplicated. The solution for now is to make a custom t2g that maps duplicates to a reference transcript. See discussion here: 
https://github.com/COMBINE-lab/pyroe/issues/41

Note that SUPPA will not be aware of this deduping, so it needs to be accounted for when mapping SUPPA events to transcript quantification results. 