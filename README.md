# VarFish DB Downloader

The purpose of this project is to collect various annotation database files
required by VarFish Web UI and convert them into the format that is required for
the import into VarFish Web UI.

## Requirements

- [conda](https://conda.io/miniconda.html)

## Installation

### Clone project

```
git clone git@cubi-gitlab.bihealth.org:CUBI_Engineering/VarFish/varfish-db-downloader
cd varfish-db-downloader
```

### Setup environment

```
conda env create -n varfish-db-downloader -f environment.yaml
conda activate varfish-db-downloader
pip install -r requirements.txt
```

## Usage

```
snakemake
```
