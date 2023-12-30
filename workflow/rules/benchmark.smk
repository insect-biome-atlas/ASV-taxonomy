import random
import os
from pathlib import Path
import pandas as pd
import subprocess
from Bio.SeqIO import parse, write as write_fasta
from tempfile import NamedTemporaryFile
from collections import defaultdict
import tqdm

localrules:
    sample_keep_species_in_db,
    sample_keep_species_remove_identical

def cluster_records(records, pid, threads):
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
    subprocess.call(
        [
            "vsearch",
            "--cluster_fast",
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
        ],
        stdout=f_null,
        stderr=f_null,
    )
    # Read file with consensus sequences
    for record in parse(cons_out.name, "fasta"):
        record.id = (record.id).replace("centroid=", "")
        clustered_records.append(record)
    cons_out.close()
    # Read file with cluster output
    uc_res = pd.read_csv(uc_out.name, sep="\t", header=None)
    uc_out.close()
    # Only keep hits and create dataframe mapping cluster reps to clustered seqs
    uc_hits = uc_res.loc[uc_res[0]=="H", [8,9]].set_index(9).rename(index = lambda x: x.split(" ")[0])
    uc_hits[8] = [x.split(" ")[0] for x in uc_hits[8]]
    f.close()
    return clustered_records, uc_hits

rule sample_keep_species_in_db:
    """
    For this case, the train fasta is the same as the original database file, so use symlink
    """
    output:
        train_tsv="benchmark/{db}/case1-keep-species-in-db/train.tsv",
        train_fasta="benchmark/{db}/case1-keep-species-in-db/train.fasta",
        test_tsv="benchmark/{db}/case1-keep-species-in-db/test.tsv",
        test_fasta="benchmark/{db}/case1-keep-species-in-db/test.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        tax = lambda wildcards: config["benchmark"][wildcards.db]["taxonomy"],
    params:
        k=100,
        seed=42
    run:
        db_fasta_abspath = os.path.abspath(input.fasta)
        db_tsv_abspath = os.path.abspath(input.tsv)
        train_fasta_abspath = os.path.abspath(output.train_fasta)
        train_tsv_abspath = os.path.abspath(output.train_tsv)
        df = pd.read_csv(input.tax, sep="\t", index_col=0)
        df.index.name="seqid"
        df = df.loc[~df["species"].str.contains("_X+$")]
        species = list(df.species.unique())
        # sample params.k species as queries
        random.seed(params.seed)
        species_sample = random.sample(species, k=params.k)
        # get subset with only sampled species
        df_subset = df.loc[df.species.isin(species_sample)]
        # For each species select one sequence at random
        queries = df_subset.groupby("species").apply(lambda x: x.sample(1, random_state=params.seed))
        queries = queries.drop("species", axis=1).reset_index().set_index("seqid")
        queries.to_csv(output.test_tsv, sep="\t")
        query_seqids = queries.index.tolist()
        # output the test seqs
        with open(output.test_fasta, 'w') as fhout:
            for record in SeqIO.parse(input.fasta, "fasta"):
                if record.id in query_seqids:
                    fhout.write(f">{record.description}\n{record.seq}\n")
        # symlink reference
        Path(train_fasta_abspath).symlink_to(db_fasta_abspath)
        Path(train_tsv_abspath).symlink_to(db_tsv_abspath)

rule sample_keep_species_remove_identical:
    """
    Hard positive case: keep species in db but remove identical sequences.
    Sample 100 species, for each species dereplicate to unique sequences, then sample 1 sequence 
    from the dereplicated sequences and remove any identical sequences from the db.
    """
    output:
        train_tsv="benchmark/{db}/case2-keep-species-remove-identical/train.tsv",
        train_fasta="benchmark/{db}/case2-keep-species-remove-identical/train.fasta",
        test_tsv="benchmark/{db}/case2-keep-species-remove-identical/test.tsv",
        test_fasta="benchmark/{db}/case2-keep-species-remove-identical/test.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        tax = lambda wildcards: config["benchmark"][wildcards.db]["taxonomy"],
    params:
        k=100,
        seed=42
    run:
        df = pd.read_csv(input.tax, sep="\t", index_col=0)
        df.index.name="seqid"
        _df = df.copy()
        # Remove sequences with ambigous species, genus and family assignments
        for r in ["family","genus","species"]:
            _df = _df.loc[~_df[rank].str.contains("_X+$")]
        # Get list of species with at least 2 sequences
        species_counts = _df.groupby("species").size()
        species = species_counts.loc[species_counts>1].index.tolist()
        # Read in sequences
        # create seqid -> species mapping
        record2species = df.loc[df.species.isin(species_sample), "species"].to_dict()
        # create dictionary to hold sequences for each species
        species_seqs = {} #defaultdict.fromkeys(species_sample, [])
        # Read the fasta file and store the records for each species sampled
        seqs = {}
        for record in tqdm.tqdm(parse("data/benchmark/coidb.221216.fasta", "fasta"), unit=" records", desc="reading fasta"):
            try:
                _ = df.loc[record.id]
                seqs[record.id] = record
            except KeyError:
                continue
            try:
                sp = record2species[record.id]
            except KeyError:
                continue
            try:
                species_seqs[sp].append(record)
            except KeyError:
                species_seqs[sp] = [record]
        df = df.loc[seqs.keys()]

        # sample 5xparams.k species
        # do this to be able to later sample params.k species with 1 sequence remaining after dereplication
        random.seed(params.seed)
        species_sample = random.sample(species, k=5*params.k)
        
        derep_species_seqs = defaultdict.fromkeys(species_sample, [])
        to_remove = []
        sampled = []
        # iterate the species and dereplicate with vsearch
        for sp in tqdm.tqdm(species_sample, unit=" species", desc="dereplicating"):
            # dereplicate the sequences so that derep_species_seqs[sp] only has unique sequences
            _derep, clusttable = cluster_records(species_seqs[sp], 1.0, 1)
            # sample 1 sequence from the dereplicated set of sequences
            random.seed(42)
            _sampled = random.sample(_derep, 1)
            _to_remove = [_sampled[0].id]
            # remove any identical sequences from the database            
            try:
                _to_remove += clusttable.loc[_sampled[0].id].values
            except KeyError:
                pass
            # check that the species still exists after removing the sampled + identical seqs
            df.loc[(~df.index.isin(_to_remove))&(df.species==sp)]
        with open(output.test_fasta, "w") as fh:
            write_fasta(records, fh, "fasta")
        with open(output.train_fasta, "w") as fh:
            write_fasta(seqs.values, fh, "fasta")
        df.loc[list(seqs.keys())].to_csv(output.train_tsv, sep="\t")
        df.loc[]

            
            

            
