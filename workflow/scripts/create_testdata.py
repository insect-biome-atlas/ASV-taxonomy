#!/usr/bin/env python

import random
import os
from pathlib import Path
import pandas as pd
import subprocess
from Bio.SeqIO import parse, write as write_fasta
from tempfile import NamedTemporaryFile
import sys
import tqdm
from argparse import ArgumentParser


def cluster_records(records, pid, threads, clust_method="cluster_fast"):
    """
    Takes a list of sequence records, writes to a temporary file and
    clusters them with vsearch.

    :param records: Sequence records
    :param pid: Percent id to cluster by
    :return:
    """
    clustered_records = []
    uc_res = {}
    if len(records) == 1:
        return records
    # Write to file
    f = NamedTemporaryFile(mode="w", delete=False)
    with open(f.name, "w") as fh:
        write_fasta(records, fh, "fasta")
    f_null = open(os.devnull, "w")
    cons_out = NamedTemporaryFile(mode="w", delete=False)
    uc_out = NamedTemporaryFile(mode="w", delete=False)
    # Run vsearch on tempfile
    cmd = [
            "vsearch",
            f"--{clust_method}",
            f.name,
            "--id",
            str(pid),
            "--consout",
            cons_out.name,
            "--uc",
            uc_out.name,
            "--notrunclabels",
            "--threads",
            str(threads),
        ]
    subprocess.call(
        cmd,
        stdout=f_null,
        stderr=f_null,
    )
    # Read file with consensus sequences
    for record in parse(cons_out.name, "fasta"):
        seqid=(record.id).split(";")[0]
        record.id = seqid.replace("centroid=", "")
        clustered_records.append(record)
    cons_out.close()
    # Read file with cluster output
    try:
        uc_res = pd.read_csv(uc_out.name, sep="\t", header=None)
    except:
        return False, False
    uc_out.close()
    # Only keep hits and create dataframe mapping cluster reps to clustered seqs
    uc_hits = uc_res.loc[uc_res[0]=="H", [8,9]].set_index(9).rename(index = lambda x: x.split(" ")[0])
    uc_hits[8] = [x.split(" ")[0] for x in uc_hits[8]]
    f.close()
    os.remove(f.name)
    os.remove(uc_out.name)
    os.remove(cons_out.name)
    return clustered_records, uc_hits

def read_fasta(fasta):
    """
    Read a fasta file and return a dictionary with seqid -> sequence mapping
    """
    seqs = {}
    for record in tqdm.tqdm(parse(fasta, "fasta"), unit=" records", desc="reading fasta", leave=False):
        seqid = (record.id).split(";")[0]
        seqs[seqid] = record
    return seqs

def sample_dereplicated(species, seqs, df, k=100, seed=42):
    sampled = []
    to_remove = []
    for sp in tqdm.tqdm(species, unit=" species", desc="sampling", total=k, leave=False):
        # get all sequences for the species
        records = [seqs[x] for x in df.loc[df.species==sp].index.tolist()]
        # dereplicate the sequences so that derep_species_seqs[sp] only has unique sequences
        _derep, clusttable = cluster_records(records, 1.0, 1)
        if _derep == False and clusttable == False:
            sys.stderr.write(f"Error clustering {sp} ({len(records)} records)\n")
            continue
        # sample 1 sequence from the dereplicated set of sequences
        random.seed(seed)
        _sampled = random.sample(_derep, 1)
        # add the sampled sequence to the ones to remove
        _to_remove = [_sampled[0].id]
        # remove any identical sequences from the database            
        try:
            _to_remove += clusttable.loc[_sampled[0].id].values
        except KeyError:
            pass
        # check that the species still exists after removing the sampled + identical seqs
        if df.loc[(~df.index.isin(_to_remove))&(df.species==sp)].shape[0] > 0:
            sampled.append(_sampled[0].id)
            to_remove += _to_remove
        if len(sampled) == k:
            break
    return sampled, to_remove


def sample_keep_genus(df, _df, species, k=100, seed=42):
    to_remove = []
    sampled = []
    for sp in tqdm.tqdm(species, unit=" species", desc="sampling", total=k, leave=False):
        g = _df.loc[_df.species==sp, "genus"].iloc[0]
        # sample 1 sequence from the species
        _df_sp = _df.loc[_df.species==sp]
        random.seed(seed)
        _sampled = random.sample(_df_sp.index.tolist(), 1)
        _to_remove = [_sampled[0]]
        # remove all other sequences from the species
        _to_remove += df.loc[df.species==sp].index.tolist()
        # check that the genus still exists even after removing the species
        if df.loc[(~df.index.isin(_to_remove+to_remove))&(df.genus==g)].shape[0] > 0:
            to_remove += _to_remove
            sampled.append(_sampled[0])
        if len(sampled) == k:
            break
    return sampled, to_remove


def sample_keep_family(df, _df, genera, k=100, seed=42):
    to_remove = []
    sampled = []
    for g in tqdm.tqdm(genera, unit=" genera", desc="sampling", total=k, leave=False):
        f = _df.loc[_df.genus==g, "family"].iloc[0]
        # sample 1 sequence from the genus
        _df_g = _df.loc[_df.genus==g]
        random.seed(seed)
        _sampled = random.sample(_df_g.index.tolist(), 1)
        _to_remove = [_sampled[0]]
        # remove all other sequences from the genus
        _to_remove += df.loc[df.genus==g].index.tolist()
        # check that the family still exists even after removing the genus
        if df.loc[(~df.index.isin(_to_remove+to_remove))&(df.family==f)].shape[0] > 0:
            to_remove += _to_remove
            sampled.append(_sampled[0])
        if len(sampled) == k:
            break
    return sampled, to_remove

def case1_sample_keep_species_in_db(args):
    db_fasta_abspath = os.path.abspath(args.input_fasta)
    train_fasta_abspath = os.path.abspath(args.output_dir+"/train.fasta")
    df = pd.read_csv(args.input_taxfile, sep="\t", index_col=0)
    df.index.name="seqid"
    sys.stderr.write(f"Read {df.shape[0]} records from {args.input_taxfile}\n")
    # Read the fasta file
    seqs = read_fasta(args.input_fasta)
    sys.stderr.write(f"Read {len(seqs)} sequences from {args.input_fasta}\n")
    # Limit taxonomy to sequences in fasta file
    df = df.loc[df.index.isin(seqs.keys())]
    _df = df.copy()
    sys.stderr.write(f"Proceeding with {df.shape[0]} records with sequences\n")
    # Remove sequences with ambigous species, genus and family assignments
    sys.stderr.write(f"Removing sequences with ambiguous taxonomy\n")
    for r in ["family","genus","species"]:
        _df = _df.loc[~_df[r].str.contains("_X+$")]
    sys.stderr.write(f"{_df.shape[0]} records remaining\n")
    species = list(_df.species.unique())
    sys.stderr.write(f"Sampling {args.k} species\n")
    # sample k species as queries
    random.seed(args.seed)
    species_sample = random.sample(species, k=args.k)
    # get subset with only sampled species
    df_subset = _df.loc[_df.species.isin(species_sample)]
    # For each species select one sequence at random
    sys.stderr.write(f"Sampling one sequence from each species\n")
    queries = df_subset.groupby("species").apply(lambda x: x.sample(1, random_state=args.seed))
    queries = queries.drop("species", axis=1).reset_index().set_index("seqid")
    sys.stderr.write(f"Writing taxonomy for test set to {args.output_dir}/test.tsv\n")
    queries.to_csv(args.output_dir+"/test.tsv", sep="\t")
    query_seqids = queries.index.tolist()
    query_records = [seqs[x] for x in query_seqids]
    # output the queries
    sys.stderr.write(f"Writing test set to {args.output_dir}/test.fasta\n")
    with open(args.output_dir+"/test.fasta", "w") as fh:
        write_fasta(query_records, fh, "fasta")
    # write/symlink reference
    sys.stderr.write(f"Symlinking reference fasta {db_fasta_abspath} to {args.output_dir}/train.fasta\n")
    Path(train_fasta_abspath).symlink_to(db_fasta_abspath)
    sys.stderr.write(f"Writing taxonomy for train set to {args.output_dir}/train.tsv\n")
    df.to_csv(args.output_dir+"/train.tsv", sep="\t")


def case2_sample_keep_species_remove_identical(args):
    """
    Hard positive case: keep species in db but remove identical sequences.
    Sample 100 species, for each species dereplicate to unique sequences, then sample 1 sequence 
    from the dereplicated sequences and remove any identical sequences from the db.
    """
    df = pd.read_csv(args.input_taxfile, sep="\t", index_col=0)
    df.index.name="seqid"
    sys.stderr.write(f"Read {df.shape[0]} records from {args.input_taxfile}\n")
    # Read the fasta file
    seqs = read_fasta(args.input_fasta)
    sys.stderr.write(f"Read {len(seqs)} sequences from {args.input_fasta}\n")
    # Limit taxonomy to sequences in fasta file
    df = df.loc[df.index.isin(seqs.keys())]
    _df = df.copy()
    sys.stderr.write(f"Proceeding with {df.shape[0]} records with sequences\n")
    # Remove sequences with ambigous species, genus and family assignments
    sys.stderr.write(f"Removing sequences with ambiguous taxonomy\n")
    for r in ["family","genus","species"]:
        _df = _df.loc[~_df[r].str.contains("_X+$")]
    sys.stderr.write(f"{_df.shape[0]} records remaining\n")
    # Get list of species with at least 2 sequences
    species_counts = _df.groupby("species").size()
    species = species_counts.loc[species_counts>1].index.tolist()
    sys.stderr.write(f"{len(species)} species with at least 2 sequences\n")
    # Iterate the species and get all sequences for each species
    sampled, to_remove = sample_dereplicated(species, seqs, _df, k=args.k, seed=args.seed)
    df = df.loc[~df.index.isin(to_remove)]
    sys.stderr.write(f"Writing taxonomy for train set to {args.output_dir}/train.tsv\n")
    df.to_csv(args.output_dir+"/train.tsv", sep="\t")
    sys.stderr.write(f"Writing taxonomy for test set to {args.output_dir}/test.tsv\n")
    _df.loc[sampled].to_csv(args.output_dir+"/test.tsv", sep="\t")
    query_records = [seqs[x] for x in sampled]
    sys.stderr.write(f"Writing test set to {args.output_dir}/test.fasta\n")
    with open(args.output_dir+"/test.fasta", "w") as fh:
        write_fasta(query_records, fh, "fasta")
    train_records = [seqs[x] for x in df.index.tolist()]
    sys.stderr.write(f"Writing train set to {args.output_dir}/train.fasta\n")
    with open(args.output_dir+"/train.fasta", "w") as fh:
        write_fasta(train_records, fh, "fasta")


def case3_remove_species_keep_genus(args):
    """
    Negative case: remove species in db, but make sure at least 1 species for the genus is present
    """
    df = pd.read_csv(args.input_taxfile, sep="\t", index_col=0)
    df.index.name="seqid"
    sys.stderr.write(f"Read {df.shape[0]} records from {args.input_taxfile}\n")
    # Read the fasta file
    seqs = read_fasta(args.input_fasta)
    sys.stderr.write(f"Read {len(seqs)} sequences from {args.input_fasta}\n")
    # Limit taxonomy to sequences in fasta file
    df = df.loc[df.index.isin(seqs.keys())]
    _df = df.copy()
    sys.stderr.write(f"Proceeding with {df.shape[0]} records with sequences\n")
    # Remove sequences with ambigous species, genus and family assignments
    sys.stderr.write(f"Removing sequences with ambiguous taxonomy\n")
    for r in ["family","genus","species"]:
        _df = _df.loc[~_df[r].str.contains("_X+$")]
    sys.stderr.write(f"{_df.shape[0]} records remaining\n")
    # Sample from genera with at least 2 sequences
    genera_counts = _df.groupby("genus").size()
    genera = genera_counts.loc[genera_counts>1].index.tolist()
    species = _df.loc[_df.genus.isin(genera), "species"].unique().tolist()
    sys.stderr.write(f"Sampling {args.k} sequences from {len(species)} species in {len(genera)} genera\n")
    sampled, to_remove = sample_keep_genus(df, _df, species, k=args.k, seed=args.seed)
    sys.stderr.write(f"Will remove {len(to_remove)} sequences from database\n")
    df_out = df.loc[~df.index.isin(to_remove)]
    _df_out = _df.loc[sampled]
    for s in sampled:
        assert _df_out.loc[s].genus in df_out.genus.unique()
        assert _df_out.loc[s].species not in df_out.species.unique()
    sys.stderr.write(f"Writing taxonomy for train set to {args.output_dir}/train.tsv\n")
    df_out.to_csv(args.output_dir+"/train.tsv", sep="\t")
    sys.stderr.write(f"Writing taxonomy for test set to {args.output_dir}/test.tsv\n")
    _df_out.loc[sampled].to_csv(args.output_dir+"/test.tsv", sep="\t")
    query_records = [seqs[x] for x in sampled]
    sys.stderr.write(f"Writing test set to {args.output_dir}/test.fasta\n")
    with open(args.output_dir+"/test.fasta", "w") as fh:
        write_fasta(query_records, fh, "fasta")
    train_records = [seqs[x] for x in df_out.index.tolist()]
    sys.stderr.write(f"Writing train set to {args.output_dir}/train.fasta\n")
    with open(args.output_dir+"/train.fasta", "w") as fh:
        write_fasta(train_records, fh, "fasta")


def case4_remove_genus_keep_family(args):
    """
    Negative case: remove genus in db, but make sure at least 1 genus for the family is present
    """
    df = pd.read_csv(args.input_taxfile, sep="\t", index_col=0)
    df.index.name="seqid"
    sys.stderr.write(f"Read {df.shape[0]} records from {args.input_taxfile}\n")
    # Read the fasta file
    seqs = read_fasta(args.input_fasta)
    sys.stderr.write(f"Read {len(seqs)} sequences from {args.input_fasta}\n")
    # Limit taxonomy to sequences in fasta file
    df = df.loc[df.index.isin(seqs.keys())]
    _df = df.copy()
    sys.stderr.write(f"Proceeding with {df.shape[0]} records with sequences\n")
    # Remove sequences with ambigous species, genus and family assignments
    sys.stderr.write(f"Removing sequences with ambiguous taxonomy\n")
    for r in ["family","genus","species"]:
        _df = _df.loc[~_df[r].str.contains("_X+$")]
    sys.stderr.write(f"{_df.shape[0]} records remaining\n")
    # Sample from families with at least 2 sequences
    family_counts = _df.groupby("family").size()
    families = family_counts.loc[family_counts>1].index.tolist()
    genera = _df.loc[_df.family.isin(families), "genus"].unique().tolist()
    sys.stderr.write(f"Sampling {args.k} sequences from {len(genera)} genera in {len(families)} families\n")
    sampled, to_remove = sample_keep_family(df, _df, genera, k=args.k, seed=args.seed)
    sys.stderr.write(f"Will remove {len(to_remove)} sequences from database\n")
    df_out = df.loc[~df.index.isin(to_remove)]
    _df_out = _df.loc[sampled]
    for s in sampled:
        assert _df_out.loc[s].family in df_out.family.unique()
        assert _df_out.loc[s].genus not in df_out.genus.unique()
    sys.stderr.write(f"Writing taxonomy for train set to {args.output_dir}/train.tsv\n")
    df_out.to_csv(args.output_dir+"/train.tsv", sep="\t")
    sys.stderr.write(f"Writing taxonomy for test set to {args.output_dir}/test.tsv\n")
    _df_out.loc[sampled].to_csv(args.output_dir+"/test.tsv", sep="\t")
    query_records = [seqs[x] for x in sampled]
    sys.stderr.write(f"Writing test set to {args.output_dir}/test.fasta\n")
    with open(args.output_dir+"/test.fasta", "w") as fh:
        write_fasta(query_records, fh, "fasta")
    train_records = [seqs[x] for x in df_out.index.tolist()]
    sys.stderr.write(f"Writing train set to {args.output_dir}/train.fasta\n")
    with open(args.output_dir+"/train.fasta", "w") as fh:
        write_fasta(train_records, fh, "fasta")


def case5_remove_family(args):
    """
    Negative case: remove all sequences for the family
    """
    df = pd.read_csv(args.input_taxfile, sep="\t", index_col=0)
    df.index.name="seqid"
    sys.stderr.write(f"Read {df.shape[0]} records from {args.input_taxfile}\n")
    # Read the fasta file
    seqs = read_fasta(args.input_fasta)
    sys.stderr.write(f"Read {len(seqs)} sequences from {args.input_fasta}\n")
    # Limit taxonomy to sequences in fasta file
    df = df.loc[df.index.isin(seqs.keys())]
    _df = df.copy()
    sys.stderr.write(f"Proceeding with {df.shape[0]} records with sequences\n")
    # Remove sequences with ambigous species, genus and family assignments
    sys.stderr.write(f"Removing sequences with ambiguous taxonomy\n")
    for r in ["family","genus","species"]:
        _df = _df.loc[~_df[r].str.contains("_X+$")]
    sys.stderr.write(f"{_df.shape[0]} records remaining\n")
    families = _df.family.unique().tolist()
    sys.stderr.write(f"Sampling {args.k} sequences from {len(families)} families\n")
    # Sample k sequences
    random.seed(args.seed)
    sampled = random.sample(_df.index.tolist(), args.k)
    desc = _df.loc[sampled].describe().loc["unique",].to_dict()
    sys.stderr.write(f"Sampled {len(sampled)} sequences in {desc['species']} species, {desc['genus']} genera and {desc['family']} families\n")
    _df_out = _df.loc[sampled]
    # Remove all sequences for the sampled families
    to_remove = df.loc[df.family.isin( _df_out.family.unique())].index.tolist()
    sys.stderr.write(f"Will remove {len(to_remove)} sequences from database\n")
    df_out = df.loc[~df.index.isin(to_remove)]
    for s in sampled:
        assert _df_out.loc[s].family not in df_out.family.unique()
    sys.stderr.write(f"Writing taxonomy for train set to {args.output_dir}/train.tsv\n")
    df_out.to_csv(args.output_dir+"/train.tsv", sep="\t")
    sys.stderr.write(f"Writing taxonomy for test set to {args.output_dir}/test.tsv\n")
    _df_out.loc[sampled].to_csv(args.output_dir+"/test.tsv", sep="\t")
    query_records = [seqs[x] for x in sampled]
    sys.stderr.write(f"Writing test set to {args.output_dir}/test.fasta\n")
    with open(args.output_dir+"/test.fasta", "w") as fh:
        write_fasta(query_records, fh, "fasta")
    train_records = [seqs[x] for x in df_out.index.tolist()]
    sys.stderr.write(f"Writing train set to {args.output_dir}/train.fasta\n")
    with open(args.output_dir+"/train.fasta", "w") as fh:
        write_fasta(train_records, fh, "fasta")


def main(args):
    if args.case == "1":
        case1_sample_keep_species_in_db(args)
    elif args.case == "2":
        case2_sample_keep_species_remove_identical(args)
    elif args.case == "3":
        case3_remove_species_keep_genus(args)
    elif args.case == "4":
        case4_remove_genus_keep_family(args)
    elif args.case == "5":
        case5_remove_family(args)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input_fasta", type=str, required=True)
    parser.add_argument("--input_taxfile", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("-k", type=int, default=100)
    parser.add_argument("-s", "--seed", type=int, default=42)
    parser.add_argument("--case", type=str, default="lower", choices=["1", "2", "3", "4", "5", "6"])
    args = parser.parse_args()
    main(args)