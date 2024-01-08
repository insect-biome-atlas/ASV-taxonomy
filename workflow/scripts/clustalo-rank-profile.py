#/!usr/bin/env python
from Bio.SeqIO import parse, write as write_fasta
from argparse import ArgumentParser
import subprocess
import tqdm
import pandas as pd
import sys
import shutil
import os
import numpy as np


def read_fasta(fasta):
    """
    Read a fasta file and return a dictionary with seqid -> sequence mapping
    """
    seqs = {}
    seqlens = {}
    for record in tqdm.tqdm(parse(fasta, "fasta"), unit=" records", desc=f"reading {fasta}", leave=False):
        seqid = (record.id).split(";")[0]
        seqs[seqid] = record
        seqlens[seqid] = len(record.seq)
    return seqs, seqlens


def profile_rank(taxfile, seqsfile, outdir, rank="family"):
    os.makedirs(outdir, exist_ok=True)
    sys.stderr.write(f"Reading {taxfile}\n")
    df = pd.read_csv(taxfile, sep="\t", index_col=0)
    seqs, seqlens = read_fasta(seqsfile)
    df = df.loc[seqs.keys()]
    # Add sequence lengths to the dataframe
    df = pd.merge(df, pd.DataFrame(seqlens, index=["bp"]).T, left_index=True, right_index=True)
    # Select 1 sequence per rank, the longest
    profile_df = df.loc[df.groupby(rank)["bp"].idxmax()]
    # Write sequences for the profile
    profile_records = [seqs[x] for x in profile_df.index.tolist()]
    sys.stderr.write(f"Writing {len(profile_records)} sequences to {outdir}/profile.fna\n")
    with open(f"{outdir}/profile.fna", 'w') as fhout:
        write_fasta(profile_records, fhout, "fasta")
    # Write the non-profile remaining sequences into files of size profile_df.shape[0]
    remaining_df = df.loc[~df.index.isin(profile_df.index.tolist())]
    # split remaining sequences into files of size profile_df.shape[0]
    chunks = np.array_split(remaining_df.index.tolist(), remaining_df.shape[0]//profile_df.shape[0])
    files = []
    for i, chunk in enumerate(tqdm.tqdm(chunks, desc=f"Writing remaining sequences to {len(chunks)} files", unit=" files"),start=1):
        chunk_records = [seqs[x] for x in chunk]
        f = f"{outdir}/file_{i}.fna"
        files.append(f)
        with open(f, 'w') as fhout:
            write_fasta(chunk_records, fhout, "fasta")
    return files
    

def main(args):
    files = profile_rank(args.taxfile, args.seqsfile, args.outdir, args.rank)
    # create first profile
    seedfile = f"{args.outdir}/profile.fna"
    logfile = f"{args.outdir}/clustalo.log"
    sys.stderr.write(f"Creating profile {args.profile} from {seedfile}, check {logfile} for progress\n")
    subprocess.run(["clustalo","-i",seedfile,"-o",args.profile,"--threads",args.threads,"-t","DNA","--outfmt=a2m", "--log", logfile, "-v"])
    # write filelist
    sys.stderr.write(f"Writing list of files to align to {args.outdir}/fileList.txt\n")
    with open(f"{args.outdir}/fileList.txt", 'w') as fhout:
        fhout.write("\n".join(files))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--taxfile", "-t", help="Taxonomy file", required=True)
    parser.add_argument("--seqsfile", "-s", help="Sequences file", required=True)
    parser.add_argument("--rank", "-r", help="Rank to split on", default="family")
    parser.add_argument("--profile", "-p", help="Output profile file", required=True)
    parser.add_argument("--outdir", "-o", help="Output directory to store files to align", required=True)
    parser.add_argument("--threads", "-T", help="Number of threads", default=8)
    args = parser.parse_args()
    main(args)