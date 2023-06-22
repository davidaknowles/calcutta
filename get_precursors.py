import gzip
import sys
from collections import namedtuple
from pathlib import Path
from importlib import reload
import collections
import calcutta
reload(calcutta)

Gene=namedtuple("Interval", ("chrom","strand","start","end"))


def get_genes(gtf_filename):

    genes = {}
    src_counter = collections.Counter()
    
    with calcutta.smart_open(gtf_filename) as f:
        for l in f:
            l = l.decode()
            if l[0] == "#": continue
            (chrom, src, feature, start, end, _, strand, _, the_rest) = l.strip().split("\t")
            
            if not feature == "gene": continue
            src_counter[src] +=1 
            try: 
                meta = dict( [ g.split(' ') for g in the_rest.split("; ") ]  )
            except ValueError as inst:
                print("Error in GTF: " + l, file=sys.stderr)
                raise(inst)
            
            meta = { k:v.strip('"') for k,v in meta.items() }

            gene_id = meta["gene_id"]
            genes[gene_id] = Gene(chrom, strand, int(start), int(end))

    print(src_counter)
    return(genes)

indexdir = Path("/gpfs/commons/groups/knowles_lab/index/kallisto/mus_musculus/with_precursor/")

genome_fasta = calcutta.get_fasta(indexdir / "GRCm39.primary_assembly.genome.fa", first_field = True)

genes = get_genes(indexdir / "gencode.vM32.basic.annotation.gtf.gz")

# check for redundant genes? 

with open(indexdir / "precursors.fa", "w") as f: 
    for gene,tup in genes.items(): 
        seq = genome_fasta[tup.chrom][tup.start:tup.end]
        if tup.strand == "-": seq = calcutta.reverse_complement(seq)
        f.write(">" + gene + "\n")
        f.write(seq + "\n")
        break