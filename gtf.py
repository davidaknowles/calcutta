import gzip
import sys
from collections import namedtuple

Interval=namedtuple("Interval", ("chr","strand","start","end","gene"))

def my_open(filename, mode):
    """
    Automatically determine whether to open a file using gzip or not
    
    Assumes gzip if the last two characters in filename are "gz"
    """
    return( gzip.open(filename,mode) if (filename[-2:] == "gz") else open(filename,mode) )

def get_exons(gtf_filename):
    """
    Load exons from GTF file
    
    Parameters
    ----------
    gtf_filename : string
        GTF file to parse
        
    Returns
    -------
    exons : set(Intervals)
    transcripts : dictionary transcript_id -> set(exons)
    genes : dictionary transcript_id -> gene_id
    """
    genes = {}
    transcripts={}
    exons=set()

    with my_open(gtf_filename,"r") as f:
        for l in f:
            l = l.decode()
            if l[0] == "#": continue
            (chrom, _, feature, start, end, _, strand, _, the_rest) = l.strip().split("\t")
            if not feature in ("exon","five_prime_utr","three_prime_utr"): continue
            #if not feature in ("exon"): continue
            try: 
                meta = dict( [ g.split(' ') for g in the_rest.split("; ") ]  )
            except ValueError as inst:
                print("Error in GTF: " + l, file=sys.stderr)
                raise(inst)
            
            meta = { k:v.strip('"') for k,v in meta.items() }

            gene_id = meta["gene_id"]
            exon=Interval(chr=chrom, strand=strand, start=int(start), end=int(end), gene=gene_id)
            if exon.start==exon.end: continue # length 0 UTR

            transcript_id = meta["transcript_id"]
            if not transcript_id in transcripts: transcripts[transcript_id]=set()
            transcripts[transcript_id].add(exon)

            genes[transcript_id] = gene_id

            exons.add(exon)
    return(exons,transcripts,genes)

def get_introns(transcripts, genes):
    """
    Construct introns from transcript models
    
    Parameters
    ----------
    transcripts : dictionary transcript_id -> set(exons)
    genes : dictionary transcript_id -> gene_id
    
    Returns
    -------
    introns: set(Intervals)
    """
    introns=set()
    for transcript_id, transcript in transcripts.items():

        # deal with UTRs (probably overkill), i.e. merge UTRs into coding part of exon
        start_ends = {} # dictionary exon_start -> exon_end
        for exon in transcript: # merge exons (intervals) with same start
            start_ends[exon.start]=max(exon.end, start_ends.get(exon.start,0))
        end_starts={}
        for start,end in start_ends.items(): # merge intervals with same end
            end_starts[end]=min(start, end_starts.get(end,start)) # min(start,start)=start
        start_ends = { v:k for k,v in end_starts.items() }

        starts = sorted(start_ends.keys())
        ends = [ start_ends[g] for g in starts ]
        an_exon=transcript.pop() # removes an exon!
        transcript.add(an_exon)
        chrom=an_exon.chr
        strand=an_exon.strand

        for exon_index in range(len(starts)-1):
            intron=Interval(chr=chrom,strand=strand,start=ends[exon_index],end=starts[exon_index+1], gene=genes[transcript_id])
            assert( (intron.end-intron.start) > 0 )
            introns.add(intron)
    return(introns)

def load_gtf(gtf_filename):
    """
    Extract exons and introns from GTF file. 
    
    Parameters
    ----------
    gtf_filename : string
        GTF file to parse
        
    Returns
    -------
    exons : set(Intervals)
    introns : set(Intervals)
    """

    (exons,transcripts,genes) = get_exons(gtf_filename)

    introns = get_introns(transcripts, genes)    

    return(exons, introns)

# gtf_filename = "/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor/gencode.vM32.basic.annotation.gtf"
gtf_filename = "/gpfs/commons/home/daknowles/knowles_lab/index/hg38/gencode.v38.basic.gtf.gz"

(exons,transcripts,genes) = get_exons(gtf_filename)

transcript_lens = {k:len(v) for k,v in transcripts.items()}

genes['ENSMUST00000194332.2'] # so there are single exon genes represented 


(exons, introns) = load_gtf(gtf_filename)


list(introns)[:5]

import numpy as np
import scipy.sparse as sp

splice_sites = set()
for intron in introns: 
    splice_sites.add( (intron.chr, intron.strand, intron.start) )
    splice_sites.add( (intron.chr, intron.strand, intron.end) )
nss = len(splice_sites) # 460924

nnz = len(introns)

splice_sites_lookup = {ss:i for i,ss in enumerate(splice_sites)}

indices = np.zeros((2,nnz), dtype=int) # cell then EC idx
for idx,intron in enumerate(introns): 
    indices[0,idx] = splice_sites_lookup[(intron.chr, intron.strand, intron.start)]
    indices[1,idx] = splice_sites_lookup[(intron.chr, intron.strand, intron.end)]

intron_connectivity = sp.coo_matrix((np.ones(nnz), indices), shape = (nss,nss))
intron_connectivity = intron_connectivity + intron_connectivity.T # not sure if necessary? 

#from scipy.sparse.csgraph import connected_components, degree
n_components, labels = sp.csgraph.connected_components(csgraph=intron_connectivity, directed=False, return_labels=True)
n_components # 216k

import pandas as pd
#pd.Series(labels).value_counts()

def sparse_sum(x, dim):
    return np.squeeze(np.asarray(x.sum(dim)))

degree = sparse_sum(intron_connectivity,0)

component_degrees = [ np.sum(degree[np.where(labels == label)[0]]) for label in np.unique(labels) ]
component_degrees = np.array(component_degrees) / 2

import matplotlib.pyplot as plt 

#plt.hist(component_degrees,30)

cluster_size_dist = pd.Series(component_degrees).value_counts()

#cluster_size_dist = cluster_size_dist.iloc[1:]
plt.scatter(cluster_size_dist.index, np.log10(cluster_size_dist))
plt.ylabel("log10(number clusters)")
plt.xlabel("number junctions in cluster")

    9834   suppaout_A3_strict.ioe
    8363   suppaout_A5_strict.ioe
   29171   suppaout_AF_strict.ioe
    5648   suppaout_AL_strict.ioe
    1595   suppaout_MX_strict.ioe
    4107   suppaout_RI_strict.ioe
   18569   suppaout_SE_strict.ioe

list(introns)[0]
