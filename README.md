# Taxonomic annotations workflow

## Overview

This workflow runs a number of tools to assign taxonomy to ASV sequences. The tools used are:

- [sintax](https://github.com/torognes/vsearch) (implemented in vsearch)
- [dada2](https://github.com/benjjneb/dada2)
- [Naive Bayes classifier](https://github.com/qiime2/qiime2) (implemented in qiime2)
- [Metabuli](https://github.com/steineggerlab/Metabuli/)
- [Kaiju](https://github.com/bioinformatics-centre/kaiju)
- [EPA-ng](https://github.com/pierrebarbera/epa-ng)

## Installation

The workflow is implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and uses [conda](https://docs.conda.io/en/latest/) and [apptainer](https://apptainer.org/) to manage software dependencies. To install the workflow, clone the repository and create a conda environment with the required dependencies:

```bash
git clone git@github.com:johnne/IBA-taxonomy.git
cd IBA-taxonomy
conda env create -f environment.yml
conda activate IBA-taxonomy
```

### Installation of Rackham

On the UPPMAX cluster Rackham, instead of using the `environment.yml` file to set up required software you can use the module system to load snakemake and biopython:

```bash
module load bioinfo-tools snakemake biopython
```

Then you can use the flag `--profile slurm` to run the workflow on the cluster. But first make sure that your NAISS project id is entered into `slurm/config.yaml` at the `default-resources:` parameter, *e.g.*:

```yaml
default-resources: "slurm_account=naiss2023-5-209"
```

## Configuration

Most tools used are configured using this config entry structure:

```yaml
<tool>:
  ref:
    <ref-name>:
        fasta: "<path-to-reference>.fasta"
        taxonomy: "<path-to-reference-taxonomy>.tsv"
        ranks: ["kingdom","phylum","class","order","family","genus","species"]
  query:
    <query-name>: "<path-to-query>.fasta"
```

## Benchmark setup

This workflow can create train and test datasets from a set of reference sequences. The rules file `workflow/rules/benchmark.smk` contains 
rules to generate 5 sets of train/test datasets. It requires a `benchmark:` entry in the configuration file that should look like this:

```yaml
benchmark:
    coidb.221216.7lvl:
        fasta: "data/benchmark/bold_clustered.sintax.7levels.fasta"
        taxonomy: "data/benchmark/coidb.221216.taxonomy.tsv.gz"
        ranks: ["kingdom","phylum","class","order","family","genus","species"]
        extract_subseq: False
        subseq_min_len: 300
        hmm_from: 78
        hmm_to: 202
        asv_seq: "AATAAATAACATAAGATTTTGA"
```

In the example above, `coidb.221216.7lvl` is a name given to this reference database. The `fasta:` entry should point to a fasta file with the
reference sequences. `taxonomy:` should point to a tab-separated file with sequence ids in the first column and taxonomic rank labels in 
subsequent columns. Using the fasta and taxonomy file 5 pairs of train/test datasets are produced. For each set, 100 species are sampled and 1 sequence picked from each and put into the `test` set. In cases 2-5, the test sequences are also removed from the `train` set.

| benchmark dataset name | description | output directory | 
| ---------------------- | ----------- | ---------------- |
| sample_keep_species_in_db | sample 100 species into test, keep all seqs in train | "benchmark/coidb.221216.7lvl/case1-keep-species-in-db/" |
| sample_keep_species_remove_identical | sample 100 species into test, remove test and identical sequences from train | "benchmark/coidb.221216.7lvl/case2-sample_keep_species_remove_identical/" |
| sample_keep_genus | sample 100 species into test, remove all sequences for sampled species from train but ensure that at least one sequence for the same genus is kept | "benchmark/coidb.221216.7lvl/case3-sample_keep_genus/" |
| sample_keep_family | sample 100 species into test, remove all sequences for sampled genera from train but ensure at least one sequence for each family is kept | "benchmark/coidb.221216.7lvl/case4-remove-genus-keep-family/" |
| sample_remove_family | sample 100 species into test, remove all sequences for sampled families from train | "benchmark/coidb.221216.7lvl/case5-remove-family/" |

If the config param `subseq:` is set to `True` the reference sequences will be trimmed to match the ASV region being studied. The region is matched using the COX1 hmm as a proxy. First you need to know the location on the HMM to which your ASVs match. This can be done by translating an ASV sequence into it's amino acid sequence, then using it as a query with `hmmscan` *e.g.* on the [HMMER website](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan). Enter the hmm start and stop coordinates in the config as `hmm_from:` and `hmm_to:`, respectively. In addition, supply the nucleotide sequence for one of the ASVs in the study using the `asv_seq:` config entry. Using this information, subsequences will be extracted from the references following these steps:

1. translate seqs to the longest open reading frame (using the Invertebrate Mitochondrial Code table)
2. hmmsearch of translated amino acid sequences against the COX1 HMM
3. reference seqs that overlap the `hmm_from:` and `hmm_to:` region are extracted and converted back to nucleotide sequences
4. multiple sequence alignments with matched nucleotide sequences together with the ASV sequence from `asv_seq:`
5. trimming of aligned sequences to only the start and stop of the ASV sequence, followed by gap removal

To benchmark `sintax`, `dada2` and the `naive-bayes` and `vsearch` classifiers, run:

```bash
snakemake --configfile config/benchmark-config.yml --profile slurm eval_all_dada2 eval_all_qiime2 eval_all_sintax
```