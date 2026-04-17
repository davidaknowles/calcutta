import pandas as pd
import collections
import plotnine as p9
import time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

import gzip

from collections import OrderedDict, namedtuple

REVERSER=str.maketrans("AGCT","TCGA")

def sparse_sum(x, dim):
    return np.squeeze(np.asarray(x.sum(dim)))

def reverse_complement(seq):
    return seq.translate(REVERSER)[::-1]

def smart_open(filename, *argene_set, **kwargene_set):
    return gzip.open(filename, *argene_set, **kwargene_set) if filename.suffix==".gz" else open(filename, *argene_set, **kwargene_set)

def get_fasta(fasta_file, first_field = False):
    with smart_open(fasta_file) as f:
        F = f.read().split(">")
    dic = OrderedDict()
    for x in F:
        x = x.split("\n")
        seq =  "".join("".join(x[1:]).split("\r"))
        seq=seq.upper()
        if len(x) <= 1: continue
        if first_field:
            dic[x[0].split()[0].strip()] = seq
        else:
            dic[x[0].strip()] = seq

    return dic



def get_transcript_gene_map(transcript_fn, t2g_fn):
    transcripts = pd.read_csv(transcript_fn, sep = "\t", names = ["transcripts"]).transcripts

    t2g = pd.read_csv(t2g_fn, sep = "\t", names = ("transcript","gene","common"))

    enst2g_dic = { row.transcript : row.gene for row in t2g.itertuples() } # maps ENST string to gene.
    
    t2g_dic = {tid:enst2g_dic.get(t,t) for tid,t in enumerate(transcripts)} # some transcripts are missing?? 
    
    return transcripts, t2g_dic

def get_transcript_lengths(fasta_file):
    cdna = get_fasta(fasta_file, first_field = True)
    return OrderedDict([(transcript,len(seq)) for transcript,seq in cdna.items()])

def read_ec(fn):
    ec_map = {} # mapping from equivalence class ID to sets of transcript IDs
    with open(fn) as f:
        for line in f:
            l = line.split()
            ec = int(l[0])
            ec_map[ec] = [int(x) for x in l[1].split(',')]
    return ec_map

def transcript_to_ec(ec_map):
    ec_reverse = {}
    for ec,transcript_list in ec_map.items(): 
        for transcript in transcript_list: 
            if not transcript in ec_reverse: 
                ec_reverse[transcript] = []
            ec_reverse[transcript].append(ec)
    return ec_reverse

def cell_ec_to_sparse(cell_ec):

    my_ecs = set().union( *[g.keys() for bc,g in cell_ec.items()] )
    
    sorted_ecs = sorted(my_ecs, key = lambda x: (len(x), min(x) if x else -1)) 
    
    ec_id_map = {ec:i for i,ec in enumerate(sorted_ecs)} # mapping from EC to idx
    BCs = list(cell_ec.keys())
    sorted_BCs = sorted(BCs) # probably don't need to do this? 

    nnz = sum( len(cell_ec[BC]) for BC in sorted_BCs )

    indices = np.zeros((2,nnz), dtype=int) # cell then EC idx
    umi_counts = np.zeros(nnz, dtype=int)

    nz_idx = 0
    for BC_idx,BC in enumerate(sorted_BCs): 
        nhere = len(cell_ec[BC])
        indices[0,nz_idx:nz_idx+nhere] = BC_idx
        ecs,counts = zip(*cell_ec[BC].items())
        indices[1,nz_idx:nz_idx+nhere] = [ ec_id_map[g] for g in ecs ]
        umi_counts[nz_idx:nz_idx+nhere] = counts
        nz_idx += nhere

    return sorted_ecs, ec_id_map, sorted_BCs, indices, umi_counts

def to_coo(x, shape = None):
    """x is a list, where each item corresponds to a different row. The item for a row gives the column IDs for the non-zero entries. """
    nnz = sum([len(g) for g in x])

    indices = np.zeros((2,nnz), dtype=int) # cell then EC idx

    nz_idx = 0
    for row_idx,col_ids in enumerate(x): 
        nhere = len(col_ids)
        indices[0,nz_idx:nz_idx+nhere] = row_idx
        indices[1,nz_idx:nz_idx+nhere] = list(col_ids) # might be a set
        nz_idx += nhere
    
    return sp.coo_matrix((np.ones(nnz), indices), shape = shape)


def read_alevin_ec(fn):
    """Read gene_eqclass.txt.gz from `alevin quant --dump-eqclasses`"""

    ecs = collections.OrderedDict()

    with smart_open(fn) as f: 
        for i,l in enumerate(f):
            if i==0: # first line gives number of features (normally genes but can be transcripts)
                num_genes = int(l.decode().strip())
                continue
            if i==1: # second line gives number of ECs
                num_ec = int(l.decode().strip())
                continue
            l = l.decode().strip().split()
            l = [int(g) for g in l]
            ec_idx = l[-1] # last element of line is EC index
            gene_idx = l[:-1]
            ecs[ec_idx] = gene_idx
    return num_genes, num_ec, ecs


def make_cell_halfcell_matrix(BCs, polydT_hex_pairs, dtype = np.float32):
    cb_to_idx = {cb:i for i,cb in enumerate(BCs)}
    polydT_to_hex = { p:h for p,h in zip(polydT_hex_pairs.polydT, polydT_hex_pairs.hex) }
    nnz = len(BCs) # number of half cells. every halfcell should be included exactly once. 
    indices = np.zeros((2,nnz), dtype=int) # cell then half-cell
    new_BCs = []
    missing_hex_bc = 0
    cell_i = 0 
    nz_idx = 0
    for i,bc in enumerate(BCs): # iterate over half cells
        rt_bc = bc[16:]
        #rt_bc = bc[:8]
        if not rt_bc in polydT_to_hex: continue # this was a hex primer

        corresponding_hex_bc = bc[:16] + polydT_to_hex[rt_bc]
        #corresponding_hex_bc =  polydT_to_hex[rt_bc] + bc[8:]

        new_BCs.append(corresponding_hex_bc)

        indices[0,nz_idx] = cell_i
        indices[1,nz_idx] = i
        nz_idx += 1

        if corresponding_hex_bc in cb_to_idx: 
            indices[0,nz_idx] = cell_i
            indices[1,nz_idx] = cb_to_idx[corresponding_hex_bc] 
            nz_idx += 1
        else: 
            missing_hex_bc += 1
        cell_i += 1

    cell_to_half_map = sp.coo_matrix((np.ones(nz_idx, dtype = dtype), indices[:,:nz_idx]))

    unused_hex = np.where(sparse_sum(cell_to_half_map,0)==0)[0]
    for i in unused_hex:
        new_BCs.append(BCs[i])
        indices[0,nz_idx] = cell_i
        indices[1,nz_idx] = i
        nz_idx += 1
        cell_i += 1
    #assert(nz_idx == nnz)
    cell_to_half_map = sp.coo_matrix((np.ones(nnz, dtype = dtype), indices))
    print("Note:",np.sum(sparse_sum(cell_to_half_map,1)==1),"half cells have no pair")
    # 32k cells lack their pair? Much worse (~140k) with first BC
    
    return np.array(new_BCs), cell_to_half_map


def get_feature_weights(features, transcript_lengths_dic,  fragment_size = 300):
    feature_lengths = np.array([transcript_lengths_dic[g] for g in features])
    eff_lens = np.array([ 
        ((g-fragment_size) if (g>fragment_size) else g) 
        for g in feature_lengths ]) # discontinous :(
    return feature_lengths, 1. / eff_lens


def EM(counts, ec_transcript_mat, w, iterations = 30):
    
    n_transcripts = len(w)
    alpha = np.full(n_transcripts,1.0/n_transcripts) # initialize to uniform
    
    alpha_w = ec_transcript_mat.copy()

    for i in range(iterations):
        alpha_w.data = (alpha * w)[ec_transcript_mat.col] # alpha_w[e,t] = ec_transcript_mat[e,t] alpha_t w_t
        ec_sums = calcutta.sparse_sum(alpha_w,1) # ec_sums[e] = sum_{t \in e} alpha_t w_t
        z = sp.diags(counts / ec_sums) @ alpha_w 
        alpha_new = calcutta.sparse_sum(z,0)
        alpha_new /= alpha_new.sum()
        print(i,np.mean(np.abs(alpha - alpha_new)),end="\r")
        alpha = alpha_new
        
    return alpha

def get_neg_hessian(pseudobulk, ec_transcript_mat, w):
    alpha_w = ec_transcript_mat.copy()
    alpha_w.data = (alpha * w)[ec_transcript_mat.col]
    ec_sums = calcutta.sparse_sum(alpha_w,1)
    diag_w = sp.diags(w)
    neg_hessian_factor = sp.diags( np.sqrt(pseudobulk) / ec_sums ) @ ec_transcript_mat @ sp.diags(w)
    neg_hessian = neg_hessian_factor.T @ neg_hessian_factor
    return neg_hessian_factor, neg_hessian
