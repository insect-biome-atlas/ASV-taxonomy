#!/usr/bin/env python

from argparse import ArgumentParser
from Bio import AlignIO
import sys

def main(args):
    counts = AlignIO.convert(args.infile, "stockholm", args.outfile, "fasta")
    sys.stderr.write(f"Converted {args.infile} to {args.outfile}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Input file in Stockholm format")
    parser.add_argument("outfile", type=str,
                        help="Output file in fasta format")
    args = parser.parse_args()
    main(args)