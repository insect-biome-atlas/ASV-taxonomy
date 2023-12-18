#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse
import re
import pandas as pd

def main(args):
    suff = {}
    if args.suffix_from:
        df = pd.read_csv(args.suffix_from, sep="\t", header=None, index_col=0, names=["seqid","taxid"])
        suff = df.to_dict(orient="index")
    seqs = {}
    for record in parse(args.input, "fasta"):
        seq_id = re.sub("ORF\d+_(.+):\d+:\d+$", r"\g<1>", (record.id).split("|")[1])
        l = len(record.seq)
        if seq_id not in seqs.keys():
            seqs[seq_id] = {"seq": str(record.seq), "id": record.id}
        else:
            if l > len(seqs[seq_id]["seq"]):
                seqs[seq_id] = {"seq": str(record.seq), "id": record.id}
    with open(args.faa_out, "w") as fhout, open(args.txt_out, "w") as fhout_txt:
        for seq_id, d in seqs.items():
            try:
                seq_id = f"{seq_id}_{suff[seq_id]['taxid']}"
            except KeyError:
                pass
            fhout.write(f">{seq_id}\n{d['seq']}\n")
            fhout_txt.write(f"{seq_id}\t{len(d['seq'])}\n")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input", type=str, help="Input fasta file")
    parser.add_argument("faa_out", type=str, help="Fasta file output")
    parser.add_argument("txt_out", type=str, help="Text output with orf lengths")
    parser.add_argument("-s", "--suffix_from", type=str, help="Suffix sequence ids with taxids given in file (tab separated)")
    args = parser.parse_args()
    main(args)