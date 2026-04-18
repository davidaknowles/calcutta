
### tealeaf

Other potentially interesting datasets: 
- BrainSTEM: A multi-resolution fetal brain atlas to assess the fidelity of human midbrain cultures https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281535
- Loss of UBE3A impacts both neuronal and non-neuronal cells in human cerebral organoidshttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253230

## Big picture: 

Kallisto gave low read count on the OG SPLIT-seq data, compared to alevin which gives a decent number (1600?) UMIs per cell as long as you do de novo barcode correction (whitelist doesn't work for some reason?)

~~The OG SPLIT-seq data gives weird cell types (Scott also saw this in PBMCs with Parse though?)~~ [I think this was a bug]

Now looking at data from "microglialess" mice, in
/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice
Note this is Parse so "splitSeqV2" by Alevin naming! 

## Understanding Split-seq and Parse barcodes

OG (Science 2019) Split-seq ("V1") and Parse (split-seq "V2") both use 3 8mer barcodes. 

Round 1 is RT: half (48) the wells are polyDT and half are random hexamer primed. The mapping for V1 is in splitseqv1_meta/oligo_hex_bc_mapping.txt. For V2 the individual BCs are in 
/gpfs/commons/groups/knowles_lab/data/sc/splitpool/Parse_meta/Parse_RT_barcodes.txt
It's not documented but the first 48 RT BCs are for polydT, the next 48 for ranhex, and that ordering gives the pairing. 

Rounds 2 and 3 are ligation. These use the same BC in both rounds and across V1 and V2. These are the union of the RT barcodes from V1, also given in 
/gpfs/commons/groups/knowles_lab/data/sc/splitpool/Parse_meta/Parse_ligation_barcodes.txt

The final 24mer BC is [rnd3 BC][rnd2 BC][rnd1 BC]. 

There is also a sneaky 4th BC corresponding to the Illumina sublibrary. alevin doesn't extract this, and surprisingly ~8% of the 3 BCs seem to collide without including the 4th BC. 
Turns out each sublib corresponds to TWO SRR IDs (not sure why they uploaded to SRA like that). Each sublib should be analyzed by alevin separately (since the cells are different) and then concatenated (see _sublib analyses). 

## Permit list

For the OG Split-seq I preferred `knee` (at least cells had a lot more reads), but for microglialess/V2 `known` looks a lot better: end up with way more cells matching with the preprocessed data. 

# Analyses

eda.ipynb: loads kallisto ECs and looks at various properties

## Done: 

- compare cell coverage distribution for (half) cells in the preproc data vs not (assume the latter will be a lot lower): yes

## TODO:

- figure out extracting sublib BC, 8% of data is doublets otherwise!  https://github.com/COMBINE-lab/salmon/issues/874
    - OK I think I have this: by running Alevin fry on each SRR ID individually and cross referencing with the preproc data, each sublib corresponds to two SRR IDs. 
- figure out combining half cells without preproc meta: the mapping is actually simple, the first 48 RT BCs are for polydT, the next 48 for ranhex, and that ordering gives the pairing. Still leaves the question of whether is it good to combine half cells for differential splicing (DS) or we should just use ranhex. 
- Now we have reasonable looking cell-type labels from the preproc we can sum EC counts for cell types -> EM -> SUPPA-esque DS. 

## Programmatic subisoform workflow

- `grouped_em.py` aggregates cells into cell-type pseudobulks or metacells and runs transcript-level EM against an EC-by-transcript compatibility matrix.
- `subisoforms.py` builds a splice-graph-style representation from the GTF by splitting exons into disjoint sub-exonic segments, marking anchor segments, and enumerating transcript subpaths between anchors as "subisoforms". It also materializes the sparse isoform-by-subisoform matrix used to collapse transcript TPMs.
- `differential_usage.py` fits Dirichlet models to subisoform compositions and performs likelihood-ratio tests across cell-types.
- `subisoform_pipeline.py` wires the full workflow together: grouped EM -> GTF-derived subisoform collapse -> Dirichlet differential usage.
- `subisoform_simulation.py` includes both a lightweight Dirichlet-multinomial simulator for the LRT itself and a GTF-driven end-to-end simulator that starts from real transcript models, samples subisoform usage, projects to isoform expression, generates EC counts, and reruns the full pipeline.
- `microglialess.py` provides a real-data loader for the microglialess Parse sublibraries, including barcode matching back to the preprocessed cell-type metadata and cross-sublib EC remapping for selected genes.
- `transcript_utils.py` resolves transcript lengths / aliases against the spliceu reference so EM weights keep working after transcript deduplication.

The intended inputs are:

1. cells x EC counts from `alevin-fry --dump-eqclasses`
2. the matching EC x transcript compatibility matrix
3. transcript IDs / effective-length weights in the same transcript order
4. a GTF describing the mature transcript models

Typical usage is:

1. aggregate cells by cell type or build metacells within each cell type
2. run transcript EM for each aggregated sample
3. collapse transcript TPMs onto subisoforms with the sparse isoform-by-subisoform matrix
4. test each anchor-pair event group for differential subisoform usage with the Dirichlet model

## Tests and notebooks

- Unit tests covering grouped EM, subisoform construction, the Dirichlet wrapper, and simulation-based operating characteristics live in `test_grouped_em.py`, `test_differential_usage.py`, and `test_subisoform_pipeline.py`.
- `subisoform_simulation.ipynb` now uses a real GTF-derived two-isoform event to simulate subisoform usage -> isoform expression -> EC counts -> full-pipeline recovery; larger end-to-end power / type I studies are available programmatically from `subisoform_simulation.py`.
- `subisoform_microglialess.ipynb` runs the full grouped EM -> subisoform collapse -> differential-usage workflow on all microglialess sublibraries for a small panel of real genes.
