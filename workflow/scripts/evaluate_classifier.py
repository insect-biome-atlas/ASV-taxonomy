#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys

def extract_taxa(row, taxcol="Taxon"):
    rank_translate = {"k__": "kingdom", "p__": "phylum", "c__": "class", "o__": "order", "f__": "family", "g__": "genus", "s__": "species"}
    data = {}
    items = row[taxcol].split("; ")
    for item in items:
        try:
            data[rank_translate[item[0:3]]] = item[3:]
        except KeyError:
            continue
    #if "species" in data.keys() and data["species"] != "" and data["genus"] != "":
    #    data["species"] = data["genus"] + " " + data["species"]
    return pd.Series(data)

def parse_qiime2(df, taxcol="Taxon", ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]):
    """
    Parses the output from QIIME2 with this format:
                                                                                    Taxon  Confidence/Consensus
    Feature ID                                                                                     
    1d17561d16d803f652b0d6bdd671c1cc  k__Bacteria; p__Firmicutes; c__Clostridia; o__...    0.999651
    ac809fd715cced98911f73f1dfb1ffb9  k__Bacteria; p__Proteobacteria; c__Alphaproteo...    0.999861
    ca96b8fd01aa6c508305bacac37da6c2  k__Bacteria; p__Fusobacteria; c__Fusobacteriia...    0.893995
    e1b32002ee69c4e3210a2a257511f96a  k__Bacteria; p__Bacteroidetes; c__Bacteroidia;...    0.992977
    d32d15407d86a2044eead0d72cd3f9e7  k__Bacteria; p__Firmicutes; c__Clostridia; o__...    0.999981

    To a dataframe with separate columns for each rank:
                                    kingdom          phylum                class            order            family          genus species
    Feature ID                                                                                                                               
    1d17561d16d803f652b0d6bdd671c1cc  Bacteria      Firmicutes           Clostridia    Clostridiales   Ruminococcaceae   Oscillospira        
    ac809fd715cced98911f73f1dfb1ffb9  Bacteria  Proteobacteria  Alphaproteobacteria  Rhodobacterales  Rhodobacteraceae     Paracoccus        
    ca96b8fd01aa6c508305bacac37da6c2  Bacteria    Fusobacteria        Fusobacteriia  Fusobacteriales  Fusobacteriaceae  Fusobacterium        
    e1b32002ee69c4e3210a2a257511f96a  Bacteria   Bacteroidetes          Bacteroidia    Bacteroidales                                         
    d32d15407d86a2044eead0d72cd3f9e7  Bacteria      Firmicutes           Clostridia    Clostridiales   Veillonellaceae                       
    """
    parsed = df.apply(extract_taxa, args=(taxcol,), axis=1).loc[:, ranks].fillna("")
    parsed.fillna("", inplace=True)
    # add confidence/consensus column
    if "Consensus" in df.columns:
        confcol = "Consensus"
    elif "Confidence" in df.columns:
        confcol = "Confidence"
    parsed = pd.merge(parsed, df.loc[:, confcol], left_index=True, right_index=True)
    return parsed

def count_assigned(df):
    eval = {}
    for col in df.columns:
            eval[col] = df.loc[(df[col]!="")&(df[col]!="unassigned")&(~df[col].str.startswith("unclassified."))].shape[0]
    eval = pd.DataFrame(eval, index=["assigned"]).T
    eval = eval.assign(n = pd.Series(df.shape[0], index=eval.index))
    eval = eval.assign(percent = eval["assigned"]/eval["n"]*100)
    return eval

def evaluate(res, taxfile, _ranks):
    tax = pd.read_csv(taxfile, sep="\t", header=0, index_col=0)
    ranks = [r for r in _ranks if r in list(set(list(tax.columns)).intersection(list(res.columns)))]
    df = pd.merge(res, tax, left_index=True, right_index=True, suffixes=["_res", "_truth"])
    eval = {}
    for rank in ranks:
        eval[rank] = {}
        eval[rank]["correct"] = df.loc[df[rank+"_res"] == df[rank+"_truth"]].shape[0]
        eval[rank]["unassigned"] = df.loc[(df[rank+"_res"]=="")|(df[rank+"_res"]=="unassigned")|(df[rank+"_res"].str.startswith("unclassified."))].shape[0]
        eval[rank]["incorrect"] = df.loc[(df[rank+"_res"]!=df[rank+"_truth"])&(df[rank+"_res"]!="")&(df[rank+"_res"]!="unassigned")&(~df[rank+"_res"].str.startswith("unclassified."))].shape[0]
    eval = pd.DataFrame(eval).T
    eval = eval.assign(n = pd.Series(df.shape[0], index=eval.index))
    return eval

def main(args):
    _ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    resultsfile = args.results
    taxfile = args.taxonomy
    classifier = args.classifier
    res = pd.read_csv(resultsfile, sep="\t", header=0, index_col=0, dtype=str)
    res.fillna("", inplace=True)
    res.rename(index = lambda x: x.split(";")[0], inplace=True)
    if res.shape[1] == 2 and res.index.name=="Feature ID":
        res = parse_qiime2(res)
    eval = {}
    if not args.taxonomy:
        eval = count_assigned(res)
    else:
        eval = evaluate(res, taxfile, _ranks)
    eval = eval.assign(classifier = pd.Series([classifier]*eval.shape[0], index=eval.index))
    eval.index.name = "rank"
    if args.output:
        eval.to_csv(args.output, sep="\t")
    else:
        eval.to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("results", type=str, help="Path to the results file")
    parser.add_argument("--taxonomy", type=str, help="Path to known taxonomy file")
    parser.add_argument("--classifier", type=str, required=True, help="Which classifier was used to generate the results")
    parser.add_argument("--output", type=str, help="Path to output file")
    args = parser.parse_args()
    main(args)