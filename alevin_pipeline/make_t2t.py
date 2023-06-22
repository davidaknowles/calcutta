import pandas as pd
from pathlib import Path

INDEX = Path("/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/")

# t2g maps transcripts to genes. Since I want the ECs to be sets of transcripts NOT genes, just copy the transcript column into the "gene" column
t2g = pd.read_csv(INDEX / "spliceu_t2g.tsv", sep="\t", names=["t","g"])
t2t = t2g.copy()
t2t.g = t2t.t
t2t.to_csv(INDEX / "t2t.tsv", sep="\t", header=False, index=False)