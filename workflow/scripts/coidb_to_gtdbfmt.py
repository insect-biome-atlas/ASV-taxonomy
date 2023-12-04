#!/usr/bin/env python3

"""
This script takes the COIDB fasta file and generates a new fasta and taxonomy
file that conforms with the gtdb_to_taxdump.py script from the Metabuli software
https://github.com/steineggerlab/Metabuli.

1. Dummy accession numbers are suffixed to the sequence ids
2. The taxonomy is either extracted from the fasta file or from an optional info
   file and written to a taxonomy file matching that of the GTDB.

Usage:
    coidb_to_metabulifmt.py <fasta> <outdir>
"""

from argparse import ArgumentParser
from Bio import SeqIO
import re

p = re.compile('([dkpcofgs]):')

def fmt_taxstring(tax):
    """
    This function takes the SINTAX format fasta header and returns a taxonomy
    string matching the GTDB format. If it finds a kingdom-level assignment, e.g.
    'k:Chordata', it will use this in place of any domain-level assignment:

    Input:
    tax=d:Animalia,k:Arthropoda,p:Insecta,c:Lepidoptera,o:Gracillariidae,f:Conopomorpha,g:Conopomorpha cramerella,s:BOLD:AAA4000
    Output:
    'd__Arthropoda;p__Insecta;c__Lepidoptera;o__Gracillariidae;f__Conopomorpha;g__Conopomorpha cramerella;s__BOLD:AAA4000'
    """
    rename = {"d__":"d__", "k__":"p__", "p__":"c__", "c__":"o__", "o__":"f__", "f__":"g__", "g__":"s__"}
    tax = tax.split("=")[1].replace("s:", "")
    bold_id = tax.split(",")[-1]
    refmt = p.sub(r'\1__', ",".join(tax.split(",")[0:-1])).replace(",",";")
    for key in reversed(rename.keys()):
        refmt = refmt.replace(key, rename[key])
    return bold_id, refmt


def fmt_seqs(f, outdir):
    """
    Iterates each sequence in the input fasta file and reformats the header into a GTDB taxonomy style definition
    Also, the taxonomy is written to a tab-separated file with updated sequence id in the first column and the taxonomy
    string in the second column
    
    This function expects the fasta header to be in the SINTAX format, e.g.:
    >GBGC18291-19;tax=d:Animalia,k:Chordata,p:Mammalia,c:Primates,o:Hominidae,f:Homo,g:Homo neanderthalensis,s:BOLD:AAA0001

    which in the tsv file would be turned into:
    
    GBGC18291-19.1  d__Animalia;k__Chordata;p__Mammalia;c__Primates;o__Hominidae;f__Homo;g__Homo neanderthalensis;s__BOLD:AAA0001

    and in the fasta file:

    >GBGC18291-19.1
    ATAGACAGACAGAGACACA...
    """
    with open(f"{outdir}/taxonomy.tsv", "w") as fh_tax, open(f"{outdir}/seqs.fasta", "w") as fh_fasta:
        for record in SeqIO.parse(f, "fasta"):
            seq_id, tax = (record.description).split(";")
            tax_items = tax.split(",")
            refmt_seqid = f"{seq_id}.1"
            bold_id, refmt_tax = fmt_taxstring(tax)
            fh_tax.write(f"{bold_id}.1\t{refmt_tax}\n")
            fh_fasta.write(f">{bold_id}.1\n{record.seq}\n")


def main(args):
    fmt_seqs(args.fasta, args.outdir)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta", help="Fasta file to be converted. Needs to be in SINTAX format, i.e. have headers such as >LPRCE088-21;tax=d:Animalia,k:Annelida,p:Clitellata,c:Crassiclitellata,o:Diporodrilidae,f:Diporodrilus,g:Diporodrilus_X,s:BOLD:AEM4827 ")
    parser.add_argument("outdir", help="Output directory")
    args = parser.parse_args()
    main(args)