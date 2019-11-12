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

## Realize data freeze

Prepare a new environment:

```
conda create -n varfish-annotator varfish-annotator-cli jannovar-cli
conda activate varfish-annotator
```

Adapt `DATA_RELEASE` and `ANNOTATOR_VERSION` to the current values.

```
DATA_RELEASE=20190820
ANNOTATOR_VERSION=0.9
```

```
# Pack server background databases
tar chzvf \
    varfish-server-background-db-$DATA_RELEASE.tar.gz \
    varfish-server-background-db-$DATA_RELEASE/
sha256sum \
    varfish-server-background-db-$DATA_RELEASE.tar.gz \
    > varfish-server-background-db-$DATA_RELEASE.tar.gz.sha256

# Prepare & pack Jannovar DB
jannovar download \
    -d hg19/refseq_curated \
    --download-dir varfish-annotator-transcripts-$DATA_RELEASE
jannovar download \
    -d hg19/ensembl \
    --download-dir varfish-annotator-transcripts-$DATA_RELEASE
tar czvf \
    varfish-annotator-transcripts-$DATA_RELEASE.tar.gz \
    varfish-annotator-transcripts-$DATA_RELEASE/*.ser
sha256sum \
    varfish-annotator-transcripts-$DATA_RELEASE.tar.gz \
    > varfish-annotator-transcripts-$DATA_RELEASE.tar.gz.sha256

# Prepare Varfish Annotator DB
varfish-annotator init-db \
    --db-release-info "varfish-annotator:v$ANNOTATOR_VERSION" \
    --db-release-info "varfish-annotator-db:r$DATA_RELEASE" \
    \
    --ref-path $DOWNLOAD/GRCh37/reference/hs37d5/hs37d5.fa \
    \
    --db-release-info "clinvar:2019-06-22" \
    --clinvar-path $DOWNLOAD/GRCh37/clinvar/latest/clinvar_tsv_main/output/clinvar_allele_trait_pairs.single.b37.tsv.gz \
    --clinvar-path $DOWNLOAD/GRCh37/clinvar/latest/clinvar_tsv_main/output/clinvar_allele_trait_pairs.multi.b37.tsv.gz \
    \
    --db-path ./varfish-annotator-db-$DATA_RELEASE \
    \
    --db-release-info "exac:r1.0" \
    --exac-path $DOWNLOAD/GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz \
    \
    --db-release-info "gnomad_exomes:r2.1" \
    $(for path in $DOWNLOAD/GRCh37/gnomAD_exomes/r2.1/download/gnomad.exomes.r2.1.sites.chr*.normalized.vcf.bgz; do \
        echo --gnomad-exomes-path $path; \
    done) \
    \
    --db-release-info "gnomad_genomes:r2.1" \
    $(for path in $DOWNLOAD/GRCh37/gnomAD_genomes/r2.1/download/gnomad.genomes.r2.1.sites.chr*.normalized.vcf.bgz; do \
        echo --gnomad-genomes-path $path; \
    done) \
    \
    --db-release-info "thousand_genomes:v3.20101123" \
    --thousand-genomes-path $DOWNLOAD/GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz \
    \
    --db-release-info "hgmd_public:ensembl_r75" \
    --hgmd-public $DOWNLOAD/GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.tsv
gzip -c \
    varfish-annotator-db-${DATA_RELEASE}.db.h2 \
    > varfish-annotator-db-${DATA_RELEASE}.db.h2.gz
sha256sum \
    varfish-annotator-db-${DATA_RELEASE}.h2.db.gz \
    > varfish-annotator-db-${DATA_RELEASE}.h2.db.gz.sha256
```
