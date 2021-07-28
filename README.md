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
conda env create -f environment.frozen-2020-10-23.yaml
conda activate varfish-db-downloader
pip install -r requirements.txt
```

## Prepare Data Release

This command will download and process all necessary files as well as link all files in a folder to build the final packages.

```
snakemake
```

The following variables can be adjusted on the commandline (either export them as variables or preceed the command).

| Variable           | Default                 |
|--------------------|-------------------------|
| `RELEASE_PATH`     | `releases`              |
| `DATA_RELEASE`     | current date (YYYYMMDD) |
| `CLINVAR_RELEASE`  | current date (YYYYMMDD) |
| `JANNOVAR_RELEASE` | current date (YYYYMMDD) |

The output can be found in the folders `GRCh37`, `GRCh38` and `releases` (default setting).

## Pack Data Release

Use the `Makefile` to pack the output of the snakemake run to finalize the data release.

```
make pack_server
make pack_annotator
make pack_jannovar
```

The following variables can be adjusted on the commandline (either export them as variables, preceed the command or edit directly in the `Makefile`).

| Variable           | Default                 |
|--------------------|-------------------------|
| `RELEASE_PATH`     | `releases`              |
| `DATA_RELEASE`     | current date (YYYYMMDD) |
| `JANNOVAR_RELEASE` | current date (YYYYMMDD) |

The output can be found in the folder `releases`:

```
releases/jannovar-db-YYYYMMDD.tar.gz.sha256
releases/jannovar-db-YYYYMMDD.tar.gz
releases/varfish-annotator-db-YYYYMMDD.tar.gz.sha256
releases/varfish-annotator-db-YYYYMMDD.tar.gz
releases/varfish-server-background-db-YYYYMMDD.tar.gz.sha256
releases/varfish-server-background-db-YYYYMMDD.tar.gz
```

## VarFish Server Compatibility Table

Information about which VarFish DB Downloader version corresponds to which data
release version and can be used with which VarFish Server version:

| VarFish DB Downloader | Data Release | VarFish Server |
|-----------------------|--------------|----------------|
| v0.2                  | 20201006     | <= v0.22.1     |


## Developer Info

- Use `wget` only and not `curl`.
  The rationale is that for the "test mode", we are overriding a single command with a helper command.
- Output files either go to `GRCh37/...`, or `GRCh38/...` or `noref/...`.
- All `*.smk` files in `snakefiles/*` are automatically included and the output files are dynamically computed from this.
  This implies we have to follow some conventions in for rule names and output file list.
    - When using the genome build as a wildcards, it must be called `{genome_build}`.
    - Lists generating output files to be included in data release must be called `result_`.
    - You can use `{chrom}` for the chromosomes `1, 2, ..., 22, X, Y`.
    - Any other wildcard in the output file must be an entry in the configuration read from `configfile: "config.yaml"`.
      The canonical example is `download_date`.
