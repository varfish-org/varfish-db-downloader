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
rule all:
    input:
        "work/genes/dbnsfp/4.4/genes.tsv.gz",
        "work/genes/dbnsfp/4.4/genes.tsv.gz.md5",
        "work/genes/ensembl/ensembl_xlink.tsv",
        "work/genes/ensembl/ensembl_xlink.tsv.md5",
        "work/genes/enst_ensg/grch37/download/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
        "work/genes/enst_ensg/grch37/download/knowntoEnsembl.txt.gz",
        "work/genes/enst_ensg/grch37/enst_ensg.tsv",
        "work/genes/enst_ensg/grch37/enst_ensg.tsv.md5",
        "work/genes/entrez/gene_info.jsonl",
        "work/genes/entrez/gene_info.jsonl.md5",
        "work/genes/gnomad/gnomad_constraints.tsv",
        "work/genes/gnomad/gnomad_constraints.tsv.md5",
        "work/genes/hgnc/hgnc_info.jsonl",
        "work/genes/hgnc/hgnc_info.jsonl.md5",
        "work/genes/hgnc/hgnc_xlink.tsv",
        "work/genes/hgnc/hgnc_xlink.tsv.md5",
        "work/genes/mim2gene/mim2gene.tsv",
        "work/genes/mim2gene/mim2gene.tsv.md5",


# ===============================================================================================
# Modular Snakefile Includes
# ===============================================================================================


# Gene-related information.
include: "rules/genes/dbnsfp.smk"
include: "rules/genes/ensembl.smk"
include: "rules/genes/gnomad.smk"
include: "rules/genes/hgnc.smk"
include: "rules/genes/ncbi.smk"
# Reference sequence--related information.
include: "rules/reference/human.smk"


# include: "rules/annos.smk"
# include: "rules/features.smk"
# include: "rules/vardbs-grch37-strucvars.smk"
# include: "rules/tracks-grch37.smk"
