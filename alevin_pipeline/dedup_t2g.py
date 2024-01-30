from pathlib import Path
import pandas as pd

indexdir = Path("/gpfs/commons/groups/knowles_lab/index/salmon/mus_spliceu/")

dups = pd.read_csv(indexdir / "duplicate_entries.tsv", sep = "\t") # note RetainedRef has duplicates (DuplicateRef does not)
t2g = pd.read_csv(indexdir / "spliceu_t2g.tsv", sep = "\t", names = ["transcript","gene"])

merged = t2g.merge(dups, left_on = "transcript", right_on = "RetainedRef" , how='inner')

# Filter out the matched entries. Same length as fasta, good
t2g_dedup = t2g[~t2g.index.isin(merged.index)]

t2g_dedup.to_csv(indexdir / "spliceu_t2g_dedup.tsv", index = False, sep = "\t", header=False)

t2t_dedup = t2g_dedup.copy()
t2t_dedup.gene = t2t_dedup.transcript
t2t_dedup.to_csv(indexdir / "spliceu_t2t_dedup.tsv", index = False, sep = "\t", header=False)

dup2t = t2g.merge(dups, how = 'left', left_on = "transcript", right_on = "DuplicateRef")
del dup2t["gene"]
del dup2t["DuplicateRef"]
dup2t.rename(columns={"RetainedRef": 'gene'}, inplace = True)
dup2t["gene"][ dup2t["gene"].isna() ] = dup2t["transcript"][ dup2t["gene"].isna() ]

dup2t.to_csv(indexdir / "t2t_dedup_v2.tsv", index = False, sep = "\t", header=False)

