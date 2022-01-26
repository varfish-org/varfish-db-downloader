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
cp config.yaml.example config.yaml
```

Adjust the variables in `config.yaml` to your requirements.

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
      Use `{chrom_no_y}` for `1, 1, ..., 22, X` without `Y`.
    - Any other wildcard in the output file must be an entry in the configuration read from `configfile: "config.yaml"`.
      The canonical example is `download_date`.

## Building Specific Tables

Some data is downloaded from sources that have a rolling release, i.e. the
files are updated without any version number. This release system tends to
update the data more frequently. This brings two issues with it: Firstly, we
can't rely on a version number to distinguish releases, and secondly, as
releases tend to happen more often, VarFish users expect the database to be
up-to-date.

The first problem is solved by replacing what would resemble a version number
in a versioned release with the date downloaded (this is why you can provide
`download_date` in the `config.yaml` file). The second issue is solved by
building your own data as we will described now on the example of HPO and OMIM.

The instructions assume that you followed the installation and have set up
the conda environment, installed the requirements and copied the `config.yaml`
file.

We do build our own tables by specifying the output files expected by Snakemake.
One can specify rule names in Snakemake, but then it is expected that they do
not contain any wildcards, which is happens not to be the case in
VarFish DB Downloader.

The variable in our Snakemake workflow to specify the download date is
`download_date`. Note that specifying this variable in the `config.yaml` is
only required when running the full workflow and has no impact when building
only specific tables.

### Update HPO and OMIM Tables

The target files to build HPO and OMIM tables are:

```
noref/mim2gene/{download_date}/Mim2geneMedgen.tsv
noref/mim2gene/{download_date}/Mim2geneMedgen.release_info
noref/hpo/{download_date}/Hpo.tsv
noref/hpo/{download_date}/Hpo.release_info
```

Replace the variable `download_date` with the current date and pass it to the
Snakemake command:

```
$ snakemake -p \
    noref/mim2gene/20220126/Mim2geneMedgen.{tsv,release_info} \
    noref/hpo/20220126/Hpo.{tsv,release_info} \
    noref/hpo/20220126/HpoName.{tsv,release_info}
```

Create a backup copy of the `import_tables.tsv` file:

```
$ cp import_tables.tsv{,.bak}
```

Replace the content of the `import_tables.tsv` file (should be tab-separated!):

```
build   table_group version
noref   hpo 20220126
noref   mim2gene  20220126
```

Import the files into VarFish:

```
$ python manage.py import_tables --tables-path /path/to/varfish-db-downloader
```
