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

The download of files can be disabled to enable a test mode.
Instead, the files in `excerpt-data` are used when `CI=true` is set in the environment.

This is done by overriding the download executables `wget` and `aria2` in the Snakemake file when `CI=true` has been set.
This again is done by overriding the `PATH` environment variable.

The files can be updated by calling

```
# varfish-db-downloader wget urls-download
```

The known URLs are managed in `download_urls.yml`.

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
