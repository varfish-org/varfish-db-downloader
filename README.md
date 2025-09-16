[![CI](https://github.com/bihealth/varfish-db-downloader/actions/workflows/main.yml/badge.svg)](https://github.com/bihealth/varfish-db-downloader/actions/workflows/main.yml)
[![check-urls](https://github.com/bihealth/varfish-db-downloader/actions/workflows/check-urls.yml/badge.svg)](https://github.com/bihealth/varfish-db-downloader/actions/workflows/check-urls.yml)

# VarFish DB Downloader

The purpose of this repository is to collect various data form public sources that is eventually used in VarFish for annotation and display to the user.
This repository contains a Snakemake workflow with supporting code for downloading the data and preparing it for being used with VarFish.

**Quick Facts**

- License: MIT
- Programming Language: Python / Snakemake

## Running

Use the utility rule `help` to get a list of all available rules:

```
# snakemake --cores=1 help
```

Run them all with `all`:

```
# snakemake --cores=1 all
```

Note that this will take a long time, use a lot of disk space, and download a lot of data.

To run on a Slurm cluster, you can use the Snakemake `--slurm` option.
See `run-slurm.sh` for an example.

## Development Setup

### Prerequisites: Install `mamba` for Conda Package Management

Install conda, ideally via [miniforge](https://github.com/conda-forge/miniforge).
A quickstart:

```
# wget -O /tmp/Mambaforge-Linux-x86_64.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
# bash /tmp/Mambaforge-Linux-x86_64.sh -b -p ~/mambaforge3 -s
# source ~/mambaforge3/bin/activate
```

### Clone Project

```
# git clone git@github.com:bihealth/varfish-db-downloader.git
# cd varfish-db-downloader
```

### Setup Environment and Install Tools

This will setup the conda environment:


```
# mamba env create --file environment.yml
# conda activate varfish-db-downloader
```

This will install the `varfish-db-downloader` tools:

```
# pip install -e .
```

## Developer Rules

### Download Commands

We use `wget` and `aria2c` only and not `curl`.
The rationale is that for the **test mode**, we are overriding the two executables with helper commands.

### Development Subsets

Besides the full output, we also build a subset of the data suitable for development.
At the moment of writing, the subset is to the BRCA1 gene only.
The rationale is that this gene and its variants are heavily annotated as breast cancer predisposition screening is a common task and users/data is plenty.

### Running in Test Mode

TODO: This part is currently disabled in CI (`Run-Test-Mode` in `.github/workflows/main.yml`).
TODO: pre-mehari workflow files are missing.

The download of files can be disabled to enable a test mode.
Instead, the files in `excerpt-data` are used when `CI=true` is set in the environment.

This is done by overriding the download executables `wget` and `aria2` in the Snakemake file when `CI=true` has been set.
This again is done by overriding the `PATH` environment variable.

The files can be updated by calling

```
# varfish-db-downloader wget urls-download
```

The known URLs are managed in `download_urls.yml`.

#### Manually maintained excerpts

Some files have to be put manually in the excerpt-data folder.
Some need to be reduced in size.

- dbNSFP(a): https://usf.box.com/shared/static/2hzcx5s6p1xui7oen16xqzndfrkt8l9l
- dbNSFP(c): https://usf.box.com/shared/static/03xsrpna0nzgrytfo2pzk326t8jad4oc
- dbscSNV: https://usf.box.com/shared/static/ffwlywsat3q5ijypvunno3rg6steqfs8
- https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2025-04-01.tsv
- https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2025-04-01.json
- https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
- https://github.com/varfish-org/clinvar-data-jsonl/releases/download/clinvar-weekly-20250410/clinvar-data-extract-vars-20250410+0.18.5.tar.gz
  - Extract the tar.gz, take the first 1000 lines of the containing files and zip them again.
- https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/gene2xml/linux64.gene2xml.gz
- https://ftp.ncbi.nih.gov/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz 

### Managing GitHub Project with Terraform

```
# export GITHUB_OWNER=bihealth
# export GITHUB_TOKEN=ghp_<thetoken>

# cd utils/terraform
# terraform init
# terraform import github_repository.varfish-db-downloader varfish-db-downloader

# terraform validate
# terraform fmt
# terraform plan
# terraform apply
```

### Uploading Data to S3

For example, as follows

```
# s5cmd --dry-run --profile ext-varfish-public \
    --endpoint-url https://ceph-s3-ext.cubi.bihealth.org \
    sync \
    'output/full/mehari/genes-txs-grch3*' \
    s3://varfish-public/public/
```

### Semantic Commits

Generally, follow [Semantic Commits v1](https://www.conventionalcommits.org/en/v1.0.0/#specification), also see [examples](https://www.conventionalcommits.org/en/v1.0.0/#examples).

Here is a list of the commit message prefixes that we use:

| prefix | description |
| ------ | ----------- |
| feat | Features |
| fix | Bug Fixes |
| perf | Performance Improvements |
| deps | Dependencies |
| revert | Reverts |
| docs | Documentation |
| style | Styles |
| chore | Miscellaneous Chores |
| refactor | Code Refactoring |
| test | Tests |
| build | Build System |
| ci | Continuous Integration |
