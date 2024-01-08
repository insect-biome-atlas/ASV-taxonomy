#!/usr/bin/env python

from Bio.SeqIO import write as write_fasta
from Bio.AlignIO import parse
from argparse import ArgumentParser
import tqdm


def main(args):
    records = {}
    for fasta in tqdm.tqdm(args.fasta, desc="reading fasta files", unit=" files"):
        for record in parse(fasta, "fasta"):
            if record.id not in records.keys():
                records[record.id] = record
    with open(args.outfile, 'w') as fhout:
        write_fasta(records.values(), fhout, "fasta")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", nargs="+", help="fasta files to concatenate")
    parser.add_argument("-o", "--outfile", help="output file")
    args = parser.parse_args()
    main(args)