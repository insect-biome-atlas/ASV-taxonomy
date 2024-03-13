#/!usr/bin/env python
from Bio.SeqIO import parse, write as write_fasta
from argparse import ArgumentParser
import subprocess
from glob import glob
import tqdm
import pandas as pd
import sys
import shutil
import os


def read_fasta(fasta):
    """
    Read a fasta file and return a dictionary with seqid -> sequence mapping
    """
    seqs = {}
    for record in tqdm.tqdm(parse(fasta, "fasta"), unit=" records", desc="reading fasta", leave=False):
        seqid = (record.id).split(";")[0]
        seqs[seqid] = record
    return seqs


def splitbyrank(taxfile, seqsfile, rank="family"):
    outbase=f"{os.path.dirname(taxfile)}/_splits"
    os.makedirs(outbase, exist_ok=True)
    sys.stderr.write(f"Reading {taxfile}\n")
    df = pd.read_csv(taxfile, sep="\t", index_col=0)
    seqs = read_fasta(seqsfile)
    df = df.loc[seqs.keys()]
    rank_counts = df.groupby(rank).size()
    rank_counts = rank_counts.sort_values(ascending=False)
    for tax in tqdm.tqdm(rank_counts.loc[rank_counts>1].index.tolist(), desc=f"Splitting {seqsfile} by {rank}", unit=f" {rank}s"):
        f=f"{outbase}/{tax}.fasta"
        if os.path.exists(f) and os.path.getsize(f) > 0:
            continue
        rank_seqs = [seqs[x] for x in df.loc[df[rank]==tax].index.tolist()]
        with open(f, 'w') as fhout:
            write_fasta(rank_seqs, fhout, "fasta")
    if rank_counts.loc[rank_counts==1].shape[0] > 0:
        f=f"{outbase}/singletons.fasta"
        rank_seqs = [seqs[x] for x in df.loc[df[rank].isin(rank_counts.loc[rank_counts==1].index.tolist())].index.tolist()]
        with open(f, 'w') as fhout:
            write_fasta(rank_seqs, fhout, "fasta")
        rank_counts = rank_counts._append(pd.Series({"singletons": len(rank_seqs)}))
    return rank_counts.loc[rank_counts>1]

def main(args):
    tmpdir = f"{os.path.dirname(args.taxfile)}/_splits"
    rank_counts = splitbyrank(args.taxfile, args.seqsfile, args.rank)
    file_order = [f"{tmpdir}/{x}.fasta" for x in rank_counts.sort_values(ascending=False).index.tolist()]
    # create first profile
    seedfile = file_order.pop(0)
    profile = f"{tmpdir}/profile.fna"
    logfile = f"{tmpdir}/clustalo.log"
    sys.stderr.write(f"Aligning fasta files, check {logfile} for progress\n")
    sys.stderr.write(f"Creating profile from {seedfile}\n")
    subprocess.run(["clustalo","-i",seedfile,"-o",profile,"--threads",args.threads,"-t","DNA","--outfmt=a2m", "--log", logfile, "-v"])
    for file in tqdm.tqdm(file_order, desc="aligning", unit=" files"):
        # align each file to the seed profile, then make the new alignment the profile for the next
        outfile = f"{profile}.tmp"
        subprocess.run(["clustalo","-i",file,"-o",outfile,"--threads",args.threads,"-t","DNA","--outfmt=a2m","--profile1",profile, "--log", logfile, "-v"])
        shutil.move(outfile, profile)
    if args.asvs:
        sys.stderr.write("Aligning ASVs to final profile\n")
        # align the ASVs to the final profile
        outfile = f"{profile}.tmp"
        subprocess.run(["clustalo","-i",args.asvs,"-o",outfile,"--threads",args.threads,"-t","DNA","--outfmt=a2m","--profile1",profile, "--log", logfile, "-v"])
        shutil.move(outfile, profile)
    shutil.move(profile, args.outfile)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--taxfile", "-t", help="Taxonomy file", required=True)
    parser.add_argument("--seqsfile", "-s", help="Sequences file", required=True)
    parser.add_argument("--rank", "-r", help="Rank to split on", default="family")
    parser.add_argument("--outfile", "-o", help="Output file", required=True)
    parser.add_argument("--threads", "-T", help="Number of threads", default=8)
    parser.add_argument("--asvs", help="ASV sequences to align to the final alignment")
    args = parser.parse_args()
    main(args)