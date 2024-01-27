#!/usr/bin/env python


"""
    Output should look like this:
    seqid1	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
    seqid2	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__

    fasta output should be:
    >seqid1
    ATGCGGGCTAGAGTAGCGAT...
"""

import pandas as pd
from Bio.SeqIO import parse
from argparse import ArgumentParser

def refmt(col, last=None):
    suffix = ";"
    p = (col.name)[0]
    if last is not None and col.name == last:
        suffix = ""
    return f"{p}__"+col+suffix

def main(args):
    df = pd.read_csv(args.tsv, sep="\t", index_col=0)
    seqs = {}
    for record in parse(args.fasta, "fasta"):
        seqid = (record.id).split(";")[0]
        seqs[seqid] = record
    common_seqs = list(set(seqs.keys()).intersection(list(df.index)))
    df = df.loc[common_seqs]
    df = df.loc[:, args.ranks]
    df = df.apply(refmt, args=(args.ranks[-1],), axis=0)
    df.to_csv(args.tsv_output, sep="\t")
    with open(args.fasta_output, 'w') as fhout:
        for seqid, record in seqs.items():
            fhout.write(f">{seqid}\n{str(record.seq)}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("tsv", help="Tab separated file with taxonomic rank labels (columns) and sequence ids (rows)")
    parser.add_argument("fasta", help="Fasta file with sequences")
    parser.add_argument("tsv_output", help="Output file name")
    parser.add_argument("fasta_output", help="Output file name")
    parser.add_argument("--ranks", nargs="+", default=["kingdom","phylum","class","order","family","genus","species"], help="Taxonomic ranks included in input")
    args = parser.parse_args()
    main(args)