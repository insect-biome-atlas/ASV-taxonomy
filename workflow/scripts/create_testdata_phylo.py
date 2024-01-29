#!/usr/bin/env python

import random
import os
from pathlib import Path
import pandas as pd
import subprocess
from Bio.SeqIO import parse, write as write_fasta
from tempfile import NamedTemporaryFile
import sys
from tqdm import tqdm
from argparse import ArgumentParser
from create_testdata import cluster_records

def read_data(db_taxfile, tree_taxfile, db_fasta, tree_fasta):
    sys.stderr.write(f"Reading {db_taxfile}\n")
    db_taxa = pd.read_csv(db_taxfile, sep="\t", header=0, index_col=0)
    sys.stderr.write(f"Reading {tree_taxfile}\n")
    tree_taxa = pd.read_csv(tree_taxfile, sep="\t", header=0, index_col=0)
    db_seqs = {}
    for record in tqdm(parse(db_fasta, "fasta"), desc="Reading db fasta", unit=" seqs", leave=False):
        seqid = (record.id).split(";")[0]
        db_seqs[seqid] = record
    tree_seqs = {}
    for record in tqdm(parse(tree_fasta, "fasta"), desc="Reading tree fasta", unit=" seqs", leave=False):
        seqid = (record.id).split(";")[0]
        record.seq = (record.seq).replace("-", "")
        tree_seqs[seqid] = record
    db_taxa = db_taxa.loc[db_seqs.keys()]
    sys.stderr.write(f"Proceeding with {db_taxa.shape[0]} records with sequences\n")
    sys.stderr.write(f"Removing sequences with ambiguous taxonomy\n")
    for r in ["family","genus","species"]:
        db_taxa = db_taxa.loc[~db_taxa[r].str.contains("_X+$")]
    tree_taxa = tree_taxa.loc[tree_seqs.keys()]
    return db_taxa, tree_taxa, db_seqs, tree_seqs


def case1_sample_keep_species_remove_identical(args):
    """
    Sample from k species that are present in the tree and db taxonomy, but skip identical sequences.
    """
    db_fasta = args.input_db_fasta
    tree_fasta = args.input_tree_fasta
    db_taxfile = args.input_db_taxfile
    tree_taxfile = args.input_tree_taxfile
    output_dir = args.output_dir
    k = args.k
    seed = args.seed
    threads = args.threads
    db_taxa, tree_taxa, db_seqs, tree_seqs = read_data(db_taxfile, tree_taxfile, db_fasta, tree_fasta)
    common_species = list(set(tree_taxa.species).intersection(db_taxa.species))
    sampled = []
    pbar = tqdm(common_species, desc="Sampling species", unit=" species", ncols=120, leave=False)
    for sp in pbar:
        pbar.set_postfix(
            {
                "sampled": len(sampled),
            }
        )
        db_sp_ids = db_taxa.loc[db_taxa.species == sp].index
        tree_sp_ids = tree_taxa.loc[tree_taxa.species == sp].index
        db_sp_records = [db_seqs[seqid] for seqid in db_sp_ids]
        tree_sp_records = [tree_seqs[seqid] for seqid in tree_sp_ids]
        uc_hits = pd.DataFrame()
        for tree_sp_record in tree_sp_records:
            clustered_records, _uc_hits = cluster_records([tree_sp_record]+db_sp_records, pid=1, threads=threads, usersort=True, clust_method="cluster_smallmem")
            uc_hits = pd.concat([uc_hits, _uc_hits])
        # remove any db sequences that are identical to tree sequences
        db_sp_ids_to_remove_idx = uc_hits.loc[[x for x in tree_sp_ids if x in uc_hits.index], 8].values
        db_sp_ids_to_remove_idx = [x.split(";")[0] for x in db_sp_ids_to_remove_idx]
        db_sp_ids_to_remove_val = list(uc_hits.loc[uc_hits[8].isin(tree_sp_ids)].index)
        db_sp_ids_to_remove_val = [x.split(";")[0] for x in db_sp_ids_to_remove_val]
        db_sp_ids_to_remove = list(set(db_sp_ids_to_remove_idx + db_sp_ids_to_remove_val))
        # these remaining ids are the ones that are not identical to any tree sequence and which we can sample from
        db_sp_ids_to_keep = list(set(db_sp_ids).difference(db_sp_ids_to_remove))
        if len(db_sp_ids_to_keep) == 0:
            continue
        # sample 1 sequence
        random.seed(seed)
        sampled+=random.sample(db_sp_ids_to_keep, 1)
        if len(sampled) == k:
            break
    sampled_seqs = [db_seqs[seqid] for seqid in sampled]
    with open(f"{output_dir}/test.fasta", "w") as f:
        write_fasta(sampled_seqs, f, "fasta")
    db_taxa.loc[sampled].to_csv(f"{output_dir}/test.tsv", sep="\t", header=True, index=True)

        
def main(args):
    if args.case == "1":
        case1_sample_keep_species_remove_identical(args)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input_db_fasta", type=str, required=True, help="Input database fasta file")
    parser.add_argument("--input_tree_fasta", type=str, required=True, help="Input tree fasta file")
    parser.add_argument("--input_db_taxfile", type=str, required=True, help="Input database taxonomy file")
    parser.add_argument("--input_tree_taxfile", type=str, required=True, help="Input tree taxonomy file")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser.add_argument("-k", type=int, default=100, help="Number of sequences to sample")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--case", type=str, default="lower", choices=["1", "2", "3", "4", "5"])
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    args = parser.parse_args()
    main(args)