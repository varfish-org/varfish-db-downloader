# Main Snakefile for the varfish-db-downloader.
#
# This Snakemake workflow allows for downlaoding the background data needed by VarFish.  The data
# is downloaded from public sources and transformed as needed.  Eventually, the data is converted
# into RocksDB databases or protobuf binary files that can be used directly by
# ``varfish-server-worker`` and is used in the backend for filtering and/or exposed to the
# user via a REST API.

from varfish_db_downloader.data_versions import DATA_VERSIONS as DV

# The prefix to use for all shell commands.
SHELL_PREFIX = "export LC_ALL=C; set -x -euo pipefail;"
# Setup the shell prefix by default.
shell.prefix(SHELL_PREFIX)


# ===============================================================================================
# Test Mode
# ===============================================================================================

import os

# Activate test mode by prepending the path to the "test-mode-bin" directory to the PATH.
if os.environ.get("CI", "false").lower() == "true":
    cwd = os.getcwd()
    old_path = os.environ["PATH"]
    os.environ["PATH"] = f"{cwd}/test-mode-bin:{old_path}"


# ===============================================================================================
# Top-Level Rules
# ===============================================================================================


## help -- print this help
rule help:
    input:
        "Snakefile",
    run:
        shell.prefix("")  # no ``set -x`` for this rule
        shell(
            r"""
        echo
        echo "=== Available Rules ==="
        echo
        for f in Snakefile $(find rules/* -name '*.smk' | sort); do
            echo "--- $f ---"
            echo
            grep '^##' $f
            echo
            grep -e '^rule' $f
            echo
        done
        """
        )


## all -- run all rules
# rule all:
#     input:
# Gene-Related Information
# "work/genes/hgnc/hgnc_info.jsonl",
# "genes/ncbi/gene_info.jsonl",
# "genes/dbnsfp/genes.tsv.gz",
# "genes/xlink/ensembl.tsv",
# "genes/xlink/hgnc.tsv",
# "genes/mim2gene/mim2gene.tsv",
# # Per-Reference Variant Annotations
# "annos/grch37/cadd/.done",
# "annos/grch37/dbnsfp-4.4a/.done",
# "annos/grch37/dbnsfp-4.4c/.done",
# "annos/grch37/dbscsnv/.done",
# "annos/grch37/helixmtdb/helixmtdb.vcf.gz",
# "annos/grch37/gnomad_mtdna/gnomad_mtdna.vcf.gz",
# "annos/grch37/ucsc_conservation/ucsc_conservation.tsv",
# "annos/grch37/dbsnp/dbsnp.vcf.gz",
# "annos/grch37/gnomad_exomes/.done",
# "annos/grch37/gnomad_genomes/.done",
# "annos/grch38/cadd/.done",
# "annos/grch38/dbnsfp-4.4a/.done",
# "annos/grch38/dbnsfp-4.4c/.done",
# "annos/grch38/gnomad_exomes/.done",
# "annos/grch38/gnomad_genomes/.done",
# "annos/grch38/gnomad_mtdna/gnomad_mtdna.vcf.gz",
# "annos/grch38/helixmtdb/helixmtdb.vcf.gz",
# # Per-Reference "Features"
# "features/grch37/tads/imr90.bed",
# "features/grch37/tads/hesc.bed",
# "features/grch37/gene_regions/refseq.bed.gz",
# "features/grch37/gene_regions/ensembl.bed.gz",
# "features/grch37/masked/repeat.bed.gz",
# "features/grch37/masked/segdup.bed.gz",
# "tracks/grch37/ucsc_genomicSuperDups.bed.gz",
# "tracks/grch37/ucsc_rmsk.bed.gz",
# "tracks/grch37/ucsc_fixSeqLiftOverPsl.bed.gz",
# "tracks/grch37/ucsc_altSeqLiftOverPsl.bed.gz",
# # "vardbs/grch37/strucvar/clinvar.bed.gz",
# "vardbs/grch37/strucvar/dbvar.bed.gz",
# "vardbs/grch37/strucvar/dgv.bed.gz",
# "vardbs/grch37/strucvar/dgv_gs.bed.gz",
# "vardbs/grch37/strucvar/g1k.bed.gz",
# "vardbs/grch37/strucvar/gnomad_sv.bed.gz",
# "vardbs/grch37/strucvar/exac.bed.gz",


# ===============================================================================================
# Modular Snakefile Includes
# ===============================================================================================


# Gene-related information.
include: "rules/genes/dbnsfp.smk"
include: "rules/genes/ensembl.smk"
include: "rules/genes/gnomad.smk"
include: "rules/genes/hgnc.smk"
include: "rules/genes/ncbi.smk"
# Refernece sequence--related information.
include: "rules/reference/human.smk"


# include: "rules/annos.smk"
# include: "rules/features.smk"
# include: "rules/vardbs-grch37-strucvars.smk"
# include: "rules/tracks-grch37.smk"
