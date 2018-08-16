# Varfish Anno

The purpose of this project is to convert databases that are required by
varfish into a format that can be easily imported (i.e. TSV files with a header
containing the column names of the corresponding varfish database table).

## Requirements

* bcftools

## Installation

### Requirements

We recommend the installation of the requirements via [conda](https://conda.io/miniconda.html):

```
conda install bcftools
```

### Clone project

```
git clone git@github.com:bihealth/varfish-anno.git
cd varfish-anno
```

## Usage

### Initialize folder structure

```
make init
```

This step creates the folder structure in `databases/`.

### Download databases

```
make download
make dl_ref
```

The downloads will be stored in `databases/<database_name>/download/`.
The reference is placed in `downloads/`.

Please note that the download routine is not sophisticated. You might want to
double check the process, especially in case something breaks. It is thought
as extensive instructions to download the required databases. If the files are
already available to you, you can place them in the corresponding download
folder and omit this step. Note that in this case the conversion scripts might
need some adaption to match the correct file name (see next section).

The download links are defined in `downloads/Makefile` and the variable names
are prefixed with `URL_`. Those variables are safe to change (if the downloaded
file contains the expected format).

The KEGG database is not automatically downloadable. Instructions are printed
to obtain the required files. They need to be placed in
`databases/kegg/downloads`.

The case files are in `.ped` format and are individual depending on your
project. You need to place them in `databases/case/download`.

Note that ExAC, gnomAD and dbSNP databases are rather large files and will take
time to download.

### Convert databases

```
make convert
```

Every script defines a `HEADER`, `INPUT` and `OUTPUT` variable, and, if needed,
a `REF` variable. The names should be self-explanatory. They are preset to the
downloaded files. You can change the `INPUT` and `REF`, if needed.

Note that ExAC, gnomAD and dbSNP databases are rather large files and
especially dbSNP will take time to convert.