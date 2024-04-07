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

def read_data(db_taxfile, tree_taxfile, db_fasta, tree_fasta, db_filter_rank=None, db_filter_taxa=[]):
    """
    Reads the database taxonomy, tree taxonomy and fasta files and returns
    dataframes with sequences identified in the fasta files
    """
    sys.stderr.write(f"Reading {db_taxfile}\n")
    db_taxa = pd.read_csv(db_taxfile, sep="\t", header=0, index_col=0)
    sys.stderr.write(f"Read {db_taxa.shape[0]} records\n")
    if db_filter_rank is not None and len(db_filter_taxa)>0:
        sys.stderr.write(f"Filtering database taxonomy for {db_filter_rank} in {','.join(db_filter_taxa)}\n")
        db_taxa = db_taxa.loc[db_taxa[db_filter_rank].isin(db_filter_taxa)]
        sys.stderr.write(f"Proceeding with {db_taxa.shape[0]} records\n")
    sys.stderr.write(f"Reading {tree_taxfile}\n")
    tree_taxa = pd.read_csv(tree_taxfile, sep="\t", header=None, index_col=0, names=["leaf","lineage"])
    db_seqs = {}
    for record in tqdm(parse(db_fasta, "fasta"), desc="Reading db fasta", unit=" seqs", leave=False):
        seqid = (record.id).split(";")[0]
        db_seqs[seqid] = record
    tree_seqs = {}
    for record in tqdm(parse(tree_fasta, "fasta"), desc="Reading tree fasta", unit=" seqs", leave=False):
        seqid = (record.id).split(";")[0]
        record.seq = (record.seq).replace("-", "")
        tree_seqs[seqid] = record
    db_taxa = db_taxa.loc[list(set(db_seqs.keys()).intersection(db_taxa.index))]
    sys.stderr.write(f"Proceeding with {db_taxa.shape[0]} records with sequences\n")
    sys.stderr.write(f"Removing sequences with ambiguous taxonomy\n")
    for r in ["family","genus","species"]:
        db_taxa = db_taxa.loc[~db_taxa[r].str.contains("_X+$")]
    tree_taxa = tree_taxa.loc[tree_seqs.keys()]
    return db_taxa, tree_taxa, db_seqs, tree_seqs

def extract_ranks(df, ranks):
    d = {}
    for row in df.iterrows():
        leaf, lineage = row
        lineage = lineage.values[0]
        _d = dict(zip(ranks, lineage.split(";")))
        d[leaf] = _d
    return pd.DataFrame(d).T


def case1_sample_keep_species_remove_identical(args):
    """
    Sample from k species that are present in the tree and db taxonomy, but skip identical sequences.
    """
    db_fasta = args.input_db_fasta
    tree_fasta = args.input_tree_fasta
    db_taxfile = args.input_db_taxfile
    db_filter_rank = args.db_filter_rank
    db_filter_taxa = args.db_filter_taxa
    tree_taxfile = args.input_tree_taxfile
    ranks = args.ranks
    output_dir = args.output_dir
    k = args.k
    seed = args.seed
    threads = args.threads
    db_taxa, tree_taxa, db_seqs, tree_seqs = read_data(db_taxfile, tree_taxfile, db_fasta, tree_fasta, db_filter_rank, db_filter_taxa)
    tree_taxa = extract_ranks(tree_taxa, ranks)
    common_species = list(set(tree_taxa.species).intersection(db_taxa.species))
    sampled = []
    # Loop through the common species (present in both the sequence database and
    # the tree)
    pbar = tqdm(common_species, desc="Sampling species", unit=" species", ncols=120, leave=False)
    for sp in pbar:
        pbar.set_postfix(
            {
                "sampled": len(sampled),
            }
        )
        # Get sequences in the reference taxonomy matching the species
        db_sp_ids = db_taxa.loc[db_taxa.species == sp].index
        # Get sequences in the tree fasta matching the species
        tree_sp_ids = tree_taxa.loc[tree_taxa.species == sp].index
        # Get all sequence records
        db_sp_records = [db_seqs[seqid] for seqid in db_sp_ids]
        tree_sp_records = [tree_seqs[seqid] for seqid in tree_sp_ids]
        # Iterate the records from the tree, and cluster each separately against all 
        # records from the reference
        uc_hits = pd.DataFrame()
        for tree_sp_record in tree_sp_records:
            clustered_records, _uc_hits = cluster_records([tree_sp_record]+db_sp_records, pid=1, threads=threads, clust_method="cluster_fast")
            uc_hits = pd.concat([uc_hits, _uc_hits])
        # remove any db sequences that are identical to tree sequences
        db_sp_ids_to_remove_idx = uc_hits.loc[[x for x in tree_sp_ids if x in uc_hits.index], 8].values
        db_sp_ids_to_remove_idx = [x.split(";")[0] for x in db_sp_ids_to_remove_idx]
        db_sp_ids_to_remove_val = list(uc_hits.loc[uc_hits[8].isin(tree_sp_ids)].index)
        db_sp_ids_to_remove_val = [x.split(";")[0] for x in db_sp_ids_to_remove_val]
        db_sp_ids_to_remove = list(set(db_sp_ids_to_remove_idx + db_sp_ids_to_remove_val))
        # these remaining ids are the ones that are not identical to any tree sequence and which we can sample from
        db_sp_ids_to_keep = list(set(db_sp_ids).difference(db_sp_ids_to_remove))
        # if no db sequences remain, continue to next species
        if len(db_sp_ids_to_keep) == 0:
            continue
        # else sample 1 sequence
        random.seed(seed)
        sampled+=random.sample(db_sp_ids_to_keep, 1)
        if len(sampled) == k:
            break
    sampled_seqs = [db_seqs[seqid] for seqid in sampled]
    with open(f"{output_dir}/test.fasta", "w") as f:
        write_fasta(sampled_seqs, f, "fasta")
    db_taxa.loc[sampled].to_csv(f"{output_dir}/test.tsv", sep="\t", header=True, index=True)


def case2_sample_keep_genus_remove_species(args):
    """
    Sample from k species that are present in the tree and db taxonomy, remove species but keep genus.
    """
    db_fasta = args.input_db_fasta
    tree_fasta = args.input_tree_fasta
    db_taxfile = args.input_db_taxfile
    db_filter_rank = args.db_filter_rank
    db_filter_taxa = args.db_filter_taxa
    tree_taxfile = args.input_tree_taxfile
    output_dir = args.output_dir
    k = args.k
    seed = args.seed
    ranks = args.ranks
    db_taxa, tree_taxa, db_seqs, tree_seqs = read_data(db_taxfile, tree_taxfile, db_fasta, tree_fasta, db_filter_rank, db_filter_taxa)
    tree_taxa = extract_ranks(tree_taxa, ranks)
    # get list of genera present in both reference db and the tree
    common_genera = list(set(tree_taxa.genus).intersection(db_taxa.genus))
    # get species in those genera in the reference database
    db_species = db_taxa.loc[db_taxa.genus.isin(common_genera)].species.unique().tolist()
    tree_species = tree_taxa.species.unique().tolist()
    # get species not in tree
    db_species_diff = list(set(db_species).difference(tree_species))
    sampled = []
    pbar = tqdm(db_species_diff, desc="Sampling from species not in tree", unit=" species", ncols=120, leave=False)
    for sp in pbar:
        pbar.set_postfix(
            {
                "sampled": len(sampled),
            }
        )
        # extract ids in reference database assigned to genus g
        db_sp_ids = db_taxa.loc[db_taxa.species == sp].index.tolist()
        if len(db_sp_ids) == 0:
            continue
        # sample 1 sequence
        random.seed(seed)
        sampled+=random.sample(db_sp_ids, 1)
        if len(sampled) == k:
            break
    for seqid in sampled:
        #if db_taxa.loc[seqid, "species"] in list(tree_taxa.species):
            #print(f"Species {seqid} is in tree species")
        #if db_taxa.loc[seqid, "genus"] not in list(tree_taxa.genus):
            #print(f"Genus {seqid} is not in tree genus")
        assert db_taxa.loc[seqid, "species"] not in tree_taxa.species.tolist()
        assert db_taxa.loc[seqid, "genus"] in tree_taxa.genus.tolist()
    sampled_seqs = [db_seqs[seqid] for seqid in sampled]
    with open(f"{output_dir}/test.fasta", "w") as f:
        write_fasta(sampled_seqs, f, "fasta")
    db_taxa.loc[sampled].to_csv(f"{output_dir}/test.tsv", sep="\t", header=True, index=True)

def case3_sample_keep_family_remove_genus(args):
    """
    Sample from k species that are present in the tree and db taxonomy, remove genus but keep family.
    """
    db_fasta = args.input_db_fasta
    tree_fasta = args.input_tree_fasta
    db_taxfile = args.input_db_taxfile
    db_filter_rank = args.db_filter_rank
    db_filter_taxa = args.db_filter_taxa
    tree_taxfile = args.input_tree_taxfile
    output_dir = args.output_dir
    k = args.k
    seed = args.seed
    ranks = args.ranks
    db_taxa, tree_taxa, db_seqs, tree_seqs = read_data(db_taxfile, tree_taxfile, db_fasta, tree_fasta, db_filter_rank, db_filter_taxa)
    tree_taxa = extract_ranks(tree_taxa, ranks)
    # get list of families present in both reference db and the tree
    common_families = list(set(tree_taxa.family).intersection(db_taxa.family))
    # get species in those families in the reference database
    db_species = db_taxa.loc[db_taxa.family.isin(common_families)].species.unique().tolist()
    db_genera = db_taxa.loc[db_taxa.family.isin(common_families)].genus.unique().tolist()
    tree_species = tree_taxa.species.unique().tolist()
    tree_genera = tree_taxa.genus.unique().tolist()
    # get species not in tree
    db_species_diff = list(set(db_species).difference(tree_species))
    # get genera not in tree
    db_genera_diff = list(set(db_genera).difference(tree_genera))
    # get intersection of species not in tree 
    # and species in db_genera_diff
    # and species matching common_families
    species = set(
                    db_species_diff
                    ).intersection(
                    db_taxa.loc[db_taxa.genus.isin(db_genera_diff)].species.unique().tolist()
                    ).intersection(
                    db_taxa.loc[db_taxa.family.isin(common_families)].species.unique().tolist()
                    )
    # species should now contain species not in the tree, and their genera should also not be in the tree
    sampled = []
    # iterate species found for db_genera_diff
    pbar = tqdm(species, desc="Sampling from species", unit=" species", ncols=120, leave=False)
    for sp in pbar:
        pbar.set_postfix(
            {
                "sampled": len(sampled),
            }
        )
        # extract ids in reference database assigned to genus g
        db_sp_ids = db_taxa.loc[db_taxa.species == sp].index.tolist()
        if len(db_sp_ids) == 0:
            continue
        # sample 1 sequence
        random.seed(seed)
        sampled+=random.sample(db_sp_ids, 1)
        if len(sampled) == k:
            break
    for seqid in sampled:
        #if db_taxa.loc[seqid, "species"] in list(tree_taxa.species):
            #print(f"Species {seqid} is in tree species")
        #if db_taxa.loc[seqid, "genus"] not in list(tree_taxa.genus):
            #print(f"Genus {seqid} is not in tree genus")
        assert db_taxa.loc[seqid, "species"] not in tree_taxa.species.tolist()
        assert db_taxa.loc[seqid, "genus"] not in tree_taxa.genus.tolist()
        assert db_taxa.loc[seqid, "family"] in tree_taxa.family.tolist()
    sampled_seqs = [db_seqs[seqid] for seqid in sampled]
    with open(f"{output_dir}/test.fasta", "w") as f:
        write_fasta(sampled_seqs, f, "fasta")
    db_taxa.loc[sampled].to_csv(f"{output_dir}/test.tsv", sep="\t", header=True, index=True)

def case4_sample_keep_order_remove_family(args):
    """
    Sample from k species that are present in the tree and db taxonomy, remove family but keep order.
    """
    db_fasta = args.input_db_fasta
    tree_fasta = args.input_tree_fasta
    db_taxfile = args.input_db_taxfile
    db_filter_rank = args.db_filter_rank
    db_filter_taxa = args.db_filter_taxa
    tree_taxfile = args.input_tree_taxfile
    output_dir = args.output_dir
    k = args.k
    seed = args.seed
    ranks = args.ranks
    db_taxa, tree_taxa, db_seqs, tree_seqs = read_data(db_taxfile, tree_taxfile, db_fasta, tree_fasta, db_filter_rank, db_filter_taxa)
    tree_taxa = extract_ranks(tree_taxa, ranks)
    # get list of orders present in both reference db and the tree
    common_orders = list(set(tree_taxa["order"]).intersection(db_taxa["order"]))
    # get species in those orders in the reference database
    db_species = db_taxa.loc[db_taxa["order"].isin(common_orders)].species.unique().tolist()
    db_genera = db_taxa.loc[db_taxa["order"].isin(common_orders)].genus.unique().tolist()
    db_families = db_taxa.loc[db_taxa["order"].isin(common_orders)].family.unique().tolist()
    tree_species = tree_taxa.species.unique().tolist()
    tree_genera = tree_taxa.genus.unique().tolist()
    tree_families = tree_taxa.family.unique().tolist()
    # get species not in tree
    db_species_diff = list(set(db_species).difference(tree_species))
    # get genera not in tree
    db_genera_diff = list(set(db_genera).difference(tree_genera))
    # get families not in tree
    db_families_diff = list(set(db_families).difference(tree_families))
    # get intersection of species not in tree 
    # and species in db_genera_diff
    # and species in db_families_diff
    # and species matching common_orders
    species = set(
                    db_species_diff
                    ).intersection(
                    db_taxa.loc[db_taxa.genus.isin(db_genera_diff)].species.unique().tolist()
                    ).intersection(
                    db_taxa.loc[db_taxa.family.isin(db_families_diff)].species.unique().tolist()
                    ).intersection(
                    db_taxa.loc[db_taxa["order"].isin(common_orders)].species.unique().tolist()
                    )
    # species should now contain species not in the tree, and their genera should also not be in the tree
    sampled = []
    # iterate species found for db_genera_diff
    pbar = tqdm(species, desc="Sampling from species", unit=" species", ncols=120, leave=False)
    for sp in pbar:
        pbar.set_postfix(
            {
                "sampled": len(sampled),
            }
        )
        # extract ids in reference database assigned to genus g
        db_sp_ids = db_taxa.loc[db_taxa.species == sp].index.tolist()
        if len(db_sp_ids) == 0:
            continue
        # sample 1 sequence
        random.seed(seed)
        sampled+=random.sample(db_sp_ids, 1)
        if len(sampled) == k:
            break
    for seqid in sampled:
        #if db_taxa.loc[seqid, "species"] in list(tree_taxa.species):
            #print(f"Species {seqid} is in tree species")
        #if db_taxa.loc[seqid, "genus"] not in list(tree_taxa.genus):
            #print(f"Genus {seqid} is not in tree genus")
        assert db_taxa.loc[seqid, "species"] not in tree_taxa.species.tolist()
        assert db_taxa.loc[seqid, "genus"] not in tree_taxa.genus.tolist()
        assert db_taxa.loc[seqid, "family"] not in tree_taxa.family.tolist()
        assert db_taxa.loc[seqid, "order"] in tree_taxa["order"].tolist()
    sampled_seqs = [db_seqs[seqid] for seqid in sampled]
    with open(f"{output_dir}/test.fasta", "w") as f:
        write_fasta(sampled_seqs, f, "fasta")
    db_taxa.loc[sampled].to_csv(f"{output_dir}/test.tsv", sep="\t", header=True, index=True)


def main(args):
    if args.case == "1":
        case1_sample_keep_species_remove_identical(args)
    elif args.case == "2":
        case2_sample_keep_genus_remove_species(args)
    elif args.case == "3":
        case3_sample_keep_family_remove_genus(args)
    elif args.case == "4":
        case4_sample_keep_order_remove_family(args)


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
    parser.add_argument("--ranks", nargs="+", default=["kingdom","phylum","class","order","family","genus","species"])
    parser.add_argument("--db_filter_rank", type=str, default=None, help="Rank to filter on (default: None)")
    parser.add_argument("--db_filter_taxa", nargs="+", default=[], help="Taxa to filter on (default: [])")
    args = parser.parse_args()
    main(args)