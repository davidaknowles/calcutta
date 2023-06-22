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

def process_bus(
    bus_file,
    global_store,
    setup_umi_store_func,
    update_umi_store_func,
    update_global_func
):
    """
    Processes a bus file and extracts info
    
    Parameters
    ----------
    bus_file: the filename of the bus file, converted to txt or txt.gz
    global_store: the datastructure where info will be stored (and then returned)
    setup_umi_store_func: function that returns an empty umi_store, datastructure that will store info about the current UMI. () -> umi_store
    update_umi_store_func: function that will update the umi_store. umi_store, ec -> umi_store, where ec is the index. 
    update_global_func: function that will update the global_store. (global_store,umi_store,current_BC) -> None
    """

    previous_bc=None
    previous_umi=None
    
    start_time = time.time()
    
    with smart_open(bus_file) as f:
        
        umi_store = setup_umi_store_func() # store info about the current UMI. e.g. corresponding set of ECs  

        for line_number,line in enumerate(f): 
            
            if line_number > 0 and line_number % 1000 == 0: 
                elapsed_time = time.time() - start_time
                print(
                    line_number, 
                    elapsed_time / line_number,
                    end = "\r")
            
            barcode,umi,ec,count = line.split() # count is ignored: PCR duplicates?
            ec = int(ec)

            if barcode == previous_bc:
                if umi != previous_umi: # i.e. same cell but new UMI
                    update_global_func(global_store,umi_store,previous_bc) 

                    previous_umi = umi
                    umi_store = setup_umi_store_func()
            else:
                if not previous_bc is None: update_global_func(global_store,umi_store,previous_bc)

                previous_bc = barcode
                previous_umi = umi

                umi_store = setup_umi_store_func()
            
            umi_store = update_umi_store_func(umi_store,ec) 
            
        update_global_func(global_store,umi_store,previous_bc)
    
    return global_store

def count_ECs_per_umi(bus_file):
    
    freq = collections.Counter()
    def update_global_func(global_store,umi_store,bc): # need to define here because can't modify an argument in a lambda
        global_store[len(umi_store)] += 1
    
    process_bus(
        bus_file,
        global_store = freq,
        setup_umi_store_func = lambda: set(), # umi_store is a set of ECs
        update_umi_store_func = lambda ec_set,ec: ec_set.union({ec}), 
        update_global_func = update_global_func,
    )
    
    return freq

def count_transcripts_per_umi(bus_file, ec_map):
    
    freq = collections.Counter()
    def update_global_func(global_store,umi_store,bc): # need to define here because can't modify an argument in a lambda
        #global_store[0 if (umi_store is None) else len(umi_store)] += 1
        global_store[len(umi_store)] += 1
    
    process_bus(
        bus_file,
        global_store = freq,
        setup_umi_store_func = lambda: None, # umi_store is a set of transcripts. Init to None so intersection_update will be valid. 
        update_umi_store_func = lambda transcript_set,ec: set(ec_map[ec]) if (transcript_set is None) else transcript_set.intersection(ec_map[ec]), # ec_map maps EC idx to transcript IDs 
        update_global_func = update_global_func,
    )
    return freq

def count_genes_per_umi(bus_file, ec_map, t2g_dic):

    freq = collections.Counter()
    def update_global_func(global_store,umi_store,bc): # need to define here because can't modify an argument in a lambda
        global_store[len(umi_store)] += 1
    def update_umi_store_func(gene_set,ec):
        genes = [ t2g_dic[tid] for tid in ec_map[ec]] # ec_map maps EC idx to transcript IDs, t2g_dic maps to genes
        return set(genes) if (gene_set is None) else gene_set.intersection(genes)
    
    process_bus(
        bus_file,
        global_store = freq,
        setup_umi_store_func = lambda: None, # umi_store is a set of genes
        update_umi_store_func = update_umi_store_func, 
        update_global_func = update_global_func,
    )
    return freq

def collapse_umis_per_cell(bus_file, ec_map, min_umi_per_cell = 0):
    
    cell_ec = collections.defaultdict(lambda: collections.defaultdict(float)) # cell barcode -> EC -> UMI count
    
    def update_global_func(cell_ec,umi_store,bc): # need to define here because can't modify an argument in a lambda
        ec = frozenset(umi_store) # would be more efficient to be building an EC lookup rather than storing the frozenset each time
        cell_ec[bc][ec] += 1
    
    process_bus(
        bus_file,
        global_store = cell_ec,
        setup_umi_store_func = lambda: None, # umi_store is a set of transcripts. Init to None so intersection_update will be valid. 
        update_umi_store_func = lambda transcript_set,ec: set(ec_map[ec]) if (transcript_set is None) else transcript_set.intersection(ec_map[ec]), 
        update_global_func = update_global_func,
    )
    
    if min_umi_per_cell > 0: 
        cell_ec = {bc:g for bc,g in cell_ec.items() if sum(g.values())>=min_umi_per_cell}
    
    return cell_ec

def collapse_umis_per_cell_lowmem(bus_file, ec_map, min_umi_per_cell = 0):
    """Iteratively builds new EC map"""
    GlobalStore = collections.namedtuple("GlobalStore", "cell_ec ec_map")
    
    global_store = GlobalStore( 
        collections.defaultdict(lambda: collections.defaultdict(float)), # cell_ec: cell barcode -> EC -> UMI count
        collections.OrderedDict()
    )
    
    def update_global_func(global_store,umi_store,bc): # need to define here because can't modify an argument in a lambda
        ec = frozenset(umi_store) # would be more efficient to be building an EC lookup rather than storing the frozenset each time
        if not ec in global_store.ec_map: 
            global_store.ec_map[ec] = len(global_store.ec_map)
        global_store.cell_ec[bc][global_store.ec_map[ec]] += 1
    
    process_bus(
        bus_file,
        global_store = global_store,
        setup_umi_store_func = lambda: None, # umi_store is a set of transcripts. Init to None so intersection_update will be valid. 
        update_umi_store_func = lambda transcript_set,ec: set(ec_map[ec]) if (transcript_set is None) else transcript_set.intersection(ec_map[ec]), 
        update_global_func = update_global_func,
    )
    
    if min_umi_per_cell > 0: 
        global_store.cell_ec = {bc:g for bc,g in global_store.cell_ec.items() if sum(g.values())>=min_umi_per_cell}
    
    return global_store.cell_ec, global_store.ec_map

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

def ecs_to_sparse(ecs):
    """ecs should be a list of sets"""
    nnz=sum( len(g) for g in ecs )

    indices = np.zeros((2,nnz), dtype=int) # cell then EC idx
    nz_idx = 0
    for ec_idx,ec in enumerate(ecs): 
        nhere = len(ec)
        indices[0,nz_idx:nz_idx+nhere] = ec_idx
        indices[1,nz_idx:nz_idx+nhere] = list(ec)
        nz_idx += nhere
    # ec_transcript_mat = torch.sparse_coo_tensor(indices, np.ones(nnz))
    return sp.coo_matrix((np.ones(nnz), indices))


def read_rad(
    rad_file, 
    stranded = True, 
    has_transcript_header = False,
    max_lines = np.inf
):
    
    t2idx = collections.OrderedDict()
    if has_transcript_header: 
        transcripts = []
        start_time = time.time()
        print("Getting transcripts")
        with smart_open(rad_file) as f: 
            for line_number,line in enumerate(f): 
                if line_number > 0 and line_number % 1000 == 0: 
                    elapsed_time = time.time() - start_time
                    print(
                        line_number, 
                        elapsed_time / line_number,
                        end = "\r")
                line = line.decode()
                ss = line.split(":")
                if len(ss)>2: break
                transcripts.append(ss[1].strip())

        t2idx = { t:i for i,t in enumerate(transcripts) } # transcript ENS to numeric idx
        num_transcripts = len(transcripts)
    
    # ideally would do this in one loop (to enable reading from stdin)
    is_sense=collections.Counter()
    prev_read_idx = None
    ec = []
    ecs = {}
    start_time = time.time()
    
    # sorting is vvveery slow so maybe try without
    cell_umi = collections.defaultdict(lambda: collections.defaultdict(lambda: [])) # CB -> UMI -> ECs. 
    # Then at the end we will intersect ECs for each UMI. 

    print("Reading rad file")
    with smart_open(rad_file) as f: 
        for line_number,line in enumerate(f): 
            if line_number > max_lines: break
            if line_number > 0 and line_number % 1000 == 0: 
                elapsed_time = time.time() - start_time
                print(
                    line_number, 
                    elapsed_time / line_number,
                    end = "\r")
            line = line.decode()
            ss = line.split("\t")
            if len(ss)==1: continue # still in transcripts
            transcript = ss[-1].strip()
            if not transcript in t2idx: t2idx[transcript] = len(t2idx)
            trans_idx = t2idx[transcript]
            read_id,hi,nh,cb,umi,sense = [ g.split(":")[1] for g in ss[:-1] ]
            read_id = int(read_id)
            sense = sense == "true" 
            is_sense[sense] += 1
            if stranded and (not sense): continue # i think this is right?! 

            ec.append(trans_idx)

            if read_id != prev_read_idx: 
                ec = frozenset(ec)
                if not ec in ecs: ecs[ec] = len(ecs)
                # should we convert UMI to integer to save space? 
                cell_umi[cb][umi].append(ecs[ec])
                ec = []

            prev_read_idx = read_id
            
    print("Collapsing UMIs")
    ec_map = {i:ec for ec,i in ecs.items() }

    new_ecs = {}
    cell_ec = collections.defaultdict(lambda: collections.Counter())
    for cb in cell_umi: 
        for umi in cell_umi[cb]: 
            umi_ecs = [ set(ec_map[ec_i]) for ec_i in cell_umi[cb][umi] ] 
            ec = frozenset(set.intersection(*umi_ecs))
            if not ec in new_ecs: new_ecs[ec] = len(new_ecs)
            cell_ec[cb][new_ecs[ec]] += 1

    return new_ecs,cell_ec,list(t2idx.keys())

def read_alevin_ec(fn):

    ecs = collections.OrderedDict()

    with smart_open(fn) as f: 
        for i,l in enumerate(f):
            if i==0: 
                num_genes = int(l.decode().strip())
                continue
            if i==1: 
                num_ec = int(l.decode().strip())
                continue
            l = l.decode().strip().split()
            l = [int(g) for g in l]
            ec_idx = l[-1]
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