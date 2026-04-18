"""Utilities for extracting exon and intron intervals from GTF files."""

from __future__ import annotations

import gzip
import sys
from collections import namedtuple
from pathlib import Path

Interval = namedtuple("Interval", ("chr", "strand", "start", "end", "gene"))


def my_open(filename: str | Path, mode: str):
    """Open plain-text or gzipped files using the same call signature."""
    path = str(filename)
    return gzip.open(path, mode) if path.endswith("gz") else open(path, mode)


def get_exons(gtf_filename: str | Path):
    """Parse exon-like features from a GTF and organize them by transcript."""
    genes = {}
    transcripts = {}
    exons = set()

    with my_open(gtf_filename, "r") as handle:
        for line in handle:
            line = line.decode() if isinstance(line, bytes) else line
            if line[0] == "#":
                continue
            chrom, _, feature, start, end, _, strand, _, the_rest = line.strip().split("\t")
            if feature not in ("exon", "five_prime_utr", "three_prime_utr"):
                continue
            try:
                meta = dict(g.split(" ") for g in the_rest.split("; "))
            except ValueError as exc:
                print("Error in GTF: " + line, file=sys.stderr)
                raise exc

            meta = {key: value.strip('"') for key, value in meta.items()}
            gene_id = meta["gene_id"]
            exon = Interval(chr=chrom, strand=strand, start=int(start), end=int(end), gene=gene_id)
            if exon.start == exon.end:
                continue

            transcript_id = meta["transcript_id"]
            if transcript_id not in transcripts:
                transcripts[transcript_id] = set()
            transcripts[transcript_id].add(exon)
            genes[transcript_id] = gene_id
            exons.add(exon)
    return exons, transcripts, genes


def get_introns(transcripts, genes):
    """Construct intron intervals implied by a transcript-to-exon mapping."""
    introns = set()
    for transcript_id, transcript in transcripts.items():
        start_ends = {}
        for exon in transcript:
            start_ends[exon.start] = max(exon.end, start_ends.get(exon.start, 0))
        end_starts = {}
        for start, end in start_ends.items():
            end_starts[end] = min(start, end_starts.get(end, start))
        start_ends = {value: key for key, value in end_starts.items()}

        starts = sorted(start_ends.keys())
        ends = [start_ends[value] for value in starts]
        an_exon = next(iter(transcript))
        chrom = an_exon.chr
        strand = an_exon.strand

        for exon_index in range(len(starts) - 1):
            intron = Interval(
                chr=chrom,
                strand=strand,
                start=ends[exon_index],
                end=starts[exon_index + 1],
                gene=genes[transcript_id],
            )
            assert (intron.end - intron.start) > 0
            introns.add(intron)
    return introns


def load_gtf(gtf_filename: str | Path):
    """Load exon and intron interval sets from a GTF file."""
    exons, transcripts, genes = get_exons(gtf_filename)
    introns = get_introns(transcripts, genes)
    return exons, introns
