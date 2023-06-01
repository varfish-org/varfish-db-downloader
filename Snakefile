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
        # genes
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
        # reference-specific annotations
        # -- background/population sequence variants and annotations thereof
        # ---- GRCh37
        "work/download/annos/grch37/seqvars/cadd/whole_genome_SNVs_inclAnno.tsv.gz",
        "work/download/annos/grch37/seqvars/cadd/InDels_inclAnno.tsv.gz",
        f"work/download/annos/grch37/seqvars/dbnsfp/{DV.dbnsfp}a/LICENSE.txt",
        f"work/download/annos/grch37/seqvars/dbnsfp/{DV.dbnsfp}c/LICENSE.txt",
        f"work/download/annos/grch37/seqvars/dbscsnv/{DV.dbscsnv}/dbscSNV{DV.dbscsnv}.chr1",
        "work/download/annos/grch37/seqvars/dbsnp/dbsnp.vcf.gz",
        "work/annos/grch37/seqvars/helixmtdb/helixmtdb.vcf.gz"
        "work/annos/grch37/seqvars/gnomad_mtdna/gnomad_mtdna.vcf.gz",
        "work/annos/grch37/seqvars/gnomad_exomes/.done",
        "work/annos/grch37/seqvars/gnomad_genomes/.done",
        # ---- GRCh38
        "work/download/annos/grch38/seqvars/cadd/whole_genome_SNVs_inclAnno.tsv.gz",
        "work/download/annos/grch38/seqvars/cadd/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz",
        # NB: dbNSFP is dual reference (for download)
        # NB: dbscSNV is dual reference (for download)
        # TODO: "work/download/annos/grch38/seqvars/dbsnp/dbsnp.vcf.gz",
        "work/annos/grch38/seqvars/helixmtdb/helixmtdb.vcf.gz"
        "work/annos/grch38/seqvars/gnomad_mtdna/gnomad_mtdna.vcf.gz",
        "work/annos/grch38/seqvars/gnomad_exomes/.done",
        "work/annos/grch38/seqvars/gnomad_genomes/.done",
        ###### "work/annos/grch37/seqvars",
        # -- background/population structural variants and annoations thereof
        # ---- GRCh37
        "work/annos/grch37/strucvars/dbvar/dbvar.bed.gz",
        "work/annos/grch37/strucvars/dgv/dgv.bed.gz",
        "work/annos/grch37/strucvars/dgv_gs/dgv_gs.bed.gz",
        "work/annos/grch37/strucvars/exac/exac.bed.gz",
        "work/annos/grch37/strucvars/g1k/g1k.bed.gz"
        "work/annos/grch37/strucvars/gnomad/gnomad_sv.bed.gz",
        # ---- GRCh38
        # TODO: "work/annos/grch37/strucvars/dbvar/dbvar.bed.gz",
        # TODO: "work/annos/grch37/strucvars/gnomad/gnomad_sv.bed.gz",
        # -- genome browser "features" (position-specific)
        # ---- GRCh37
        # ---- GRCh38
        # TODO


# ===============================================================================================
# Modular Snakefile Includes
# ===============================================================================================


# Gene-related rules.
include: "rules/genes/dbnsfp.smk"
include: "rules/genes/ensembl.smk"
include: "rules/genes/gnomad.smk"
include: "rules/genes/hgnc.smk"
include: "rules/genes/ncbi.smk"
# Reference sequence--related rules.
include: "rules/reference/human.smk"
# Features (position and not variant specific).
include: "rules/annos/features/cons.smk"
include: "rules/annos/features/ensembl.smk"
include: "rules/annos/features/refseq.smk"
include: "rules/annos/features/tads.smk"
include: "rules/annos/features/ucsc.smk"
# Sequence variants and annotations.
include: "rules/annos/seqvars/cadd.smk"
include: "rules/annos/seqvars/dbnsfp.smk"
include: "rules/annos/seqvars/dbscsnv.smk"
include: "rules/annos/seqvars/dbsnp.smk"
include: "rules/annos/seqvars/gnomad_mtdna.smk"
include: "rules/annos/seqvars/gnomad_nuclear.smk"
include: "rules/annos/seqvars/helix.smk"
# Structural variant related.
include: "rules/annos/strucvars/dbvar.smk"
include: "rules/annos/strucvars/dgv.smk"
include: "rules/annos/strucvars/exac.smk"
include: "rules/annos/strucvars/g1k.smk"
include: "rules/annos/strucvars/gnomad.smk"
