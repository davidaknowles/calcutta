import pandas as pd
import collections
import plotnine as p9
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import calcutta
import scipy.sparse as sp
from importlib import reload
reload(calcutta)

salmon_dir = Path("/gpfs/commons/groups/knowles_lab/data/parse/original_SPLITseq_GSE110823/salmon_joint/out_permit_knee")

# it seems the rad file is sorted by CB so probably could do this memory efficiently (read in cell data, then summarize)
new_ecs,cell_ec,transcripts = calcutta.read_rad(
    salmon_dir / "map.collated.txt.gz", 
    stranded = False)

np.all(np.array(list(new_ecs.values())) == np.arange(len(new_ecs)))

# is it worth recoding? probably easier from the sparse matrix? 

with open(salmon_dir / "transcripts.txt","w") as f: 
    f.writelines([g+"\n" for g  in transcripts])
    
with open(salmon_dir / "my.ec","w") as f: 
    for g in new_ecs.keys(): 
        f.write(",".join([str(h) for h in g]) + "\n")
        
nnz = sum( len(cell_ec[BC]) for BC in cell_ec )

indices = np.zeros((2,nnz), dtype=int) # cell then EC idx
umi_counts = np.zeros(nnz, dtype=int)

nz_idx = 0
for BC_idx,BC in enumerate(cell_ec): 
    nhere = len(cell_ec[BC])
    indices[0,nz_idx:nz_idx+nhere] = BC_idx
    indices[1,nz_idx:nz_idx+nhere] = list(cell_ec[BC].keys()) # EC idx
    umi_counts[nz_idx:nz_idx+nhere] = list(cell_ec[BC].values())
    nz_idx += nhere
    
cell_ec_sparse = sp.coo_matrix((umi_counts, indices))

sp.save_npz(salmon_dir / "cell_umi.npz", cell_ec_sparse)