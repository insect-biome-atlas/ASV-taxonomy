# Taxonomic annotations workflow

## Overview

This workflow runs a number of tools to assign taxonomy to ASV sequences. The tools used are:

- [sintax](https://github.com/torognes/vsearch) (implemented in vsearch)
- [Metabuli](https://github.com/steineggerlab/Metabuli/)
- [Kaiju](https://github.com/bioinformatics-centre/kaiju)
- [EPA-ng](https://github.com/pierrebarbera/epa-ng)

## Installation

The workflow is implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and uses [conda](https://docs.conda.io/en/latest/) to manage software dependencies. To install the workflow, clone the repository and create a conda environment with the required dependencies:

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

