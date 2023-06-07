# Main Snakefile for the varfish-db-downloader.
#
# This Snakemake workflow allows for downlaoding the background data needed by VarFish.  The data
# is downloaded from public sources and transformed as needed.  Eventually, the data is converted
# into RocksDB databases or protobuf binary files that can be used directly by
# ``varfish-server-worker`` and is used in the backend for filtering and/or exposed to the
# user via a REST API.

from varfish_db_downloader.versions import DATA_VERSIONS as DV, PACKAGE_VERSIONS as PV

# The prefix to use for all shell commands.
SHELL_PREFIX = "export LC_ALL=C; set -x -euo pipefail;"
# Setup the shell prefix by default.
shell.prefix(SHELL_PREFIX)

# Regular expression for genome release.
RE_GENOME = r"grch(37|38)"
# Regular expression for versions.
RE_VERSION = r"\d+(\.\d+)*"

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
        # == work directory =====================================================================
        #
        # genes
        f"work/genes/dbnsfp/{DV.dbnsfp}/genes.tsv.gz",
        f"work/genes/ensembl/{DV.ensembl}/ensembl_xlink.tsv",
        f"work/genes/enst_ensg/grch37/{DV.ensembl_37}/enst_ensg.tsv",
        f"work/genes/entrez/{DV.today}/gene_info.jsonl",
        f"work/genes/gnomad/{DV.gnomad_constraints}/gnomad_constraints.tsv",
        f"work/genes/hgnc/{DV.today}/hgnc_info.jsonl",
        f"work/genes/hgnc/{DV.today}/hgnc_xlink.tsv",
        f"work/genes/mim2gene/{DV.today}/mim2gene.tsv",
        # reference-specific annotations
        # -- background/population sequence variants and annotations thereof
        # ---- GRCh37
        f"work/download/annos/grch37/seqvars/cadd/{DV.cadd}/whole_genome_SNVs_inclAnno.tsv.gz",
        f"work/download/annos/grch37/seqvars/cadd/{DV.cadd}/InDels_inclAnno.tsv.gz",
        f"work/download/annos/grch37/seqvars/dbnsfp/{DV.dbnsfp}a/LICENSE.txt",
        f"work/download/annos/grch37/seqvars/dbnsfp/{DV.dbnsfp}c/LICENSE.txt",
        f"work/download/annos/grch37/seqvars/dbscsnv/{DV.dbscsnv}/dbscSNV{DV.dbscsnv}.chr1",
        f"work/download/annos/grch37/seqvars/dbsnp/{DV.dbsnp}/dbsnp.vcf.gz",
        f"work/annos/grch37/seqvars/helixmtdb/{DV.helixmtdb}/helixmtdb.vcf.gz",
        f"work/annos/grch37/seqvars/gnomad_mtdna/{DV.gnomad_mtdna}/gnomad_mtdna.vcf.gz",
        f"work/download/annos/grch37/seqvars/gnomad_exomes/{DV.gnomad_v2}/.done",
        f"work/download/annos/grch37/seqvars/gnomad_genomes/{DV.gnomad_v2}/.done",
        # ---- GRCh38
        f"work/download/annos/grch38/seqvars/cadd/{DV.cadd}/whole_genome_SNVs_inclAnno.tsv.gz",
        f"work/download/annos/grch38/seqvars/cadd/{DV.cadd}/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz",
        # NB: dbNSFP is dual reference (for download)
        # NB: dbscSNV is dual reference (for download)
        f"work/download/annos/grch37/seqvars/dbsnp/{DV.dbsnp}/dbsnp.vcf.gz",
        f"work/annos/grch38/seqvars/helixmtdb/{DV.helixmtdb}/helixmtdb.vcf.gz",
        f"work/annos/grch38/seqvars/gnomad_mtdna/{DV.gnomad_mtdna}/gnomad_mtdna.vcf.gz",
        f"work/download/annos/grch38/seqvars/gnomad_exomes/{DV.gnomad_v2}/.done",
        f"work/download/annos/grch38/seqvars/gnomad_genomes/{DV.gnomad_v3}/.done",
        # -- background/population structural variants and annoations thereof
        # ---- GRCh37
        f"work/annos/grch37/strucvars/dbvar/{DV.dbvar}/dbvar.bed.gz",
        f"work/annos/grch37/strucvars/dgv/{DV.dgv}/dgv.bed.gz",
        f"work/annos/grch37/strucvars/dgv_gs/{DV.dgv_gs}/dgv_gs.bed.gz",
        f"work/annos/grch37/strucvars/exac/{DV.exac_cnv}/exac.bed.gz",
        f"work/annos/grch37/strucvars/g1k/{DV.g1k_svs}/g1k.bed.gz",
        f"work/annos/grch37/strucvars/gnomad/{DV.gnomad_sv}/gnomad_sv.bed.gz",
        # ---- GRCh38
        f"work/annos/grch38/strucvars/dbvar/{DV.dbvar}/dbvar.bed.gz",
        f"work/annos/grch38/strucvars/dgv/{DV.dgv}/dgv.bed.gz",
        f"work/annos/grch38/strucvars/dgv_gs/{DV.dgv_gs}/dgv_gs.bed.gz",
        # NB: gnomAD-SV GRCh38 was announced end of 2020 but not released yet
        # -- genome browser "features" (position-specific)
        # ---- GRCh37
        f"work/annos/grch37/features/cons/{DV.ucsc_cons_37}/ucsc_conservation.tsv",
        f"work/annos/grch37/features/ensembl/{DV.ensembl_37}/ensembl_genes.bed.gz",
        f"work/annos/grch37/features/refseq/{DV.refseq_37}/refseq_genes.bed.gz",
        "work/annos/grch37/features/tads/dixon2015/hesc.bed",
        f"work/annos/grch37/features/ucsc/{DV.ucsc_genomic_super_dups_37}/genomicSuperDups.bed.gz",
        f"work/annos/grch37/features/ucsc/{DV.ucsc_rmsk_37}/rmsk.bed.gz",
        f"work/annos/grch37/features/ucsc/{DV.ucsc_alt_seq_liftover_37}/altSeqLiftOverPsl.bed.gz",
        f"work/annos/grch37/features/ucsc/{DV.ucsc_fix_seq_liftover_37}/fixSeqLiftOverPsl.bed.gz",
        # ---- GRCh38
        f"work/annos/grch38/features/cons/{DV.ucsc_cons_38}/ucsc_conservation.tsv",
        f"work/annos/grch38/features/ensembl/{DV.ensembl_38}/ensembl_genes.bed.gz",
        f"work/annos/grch38/features/refseq/{DV.refseq_38}/refseq_genes.bed.gz",
        "work/annos/grch38/features/tads/dixon2015/hesc.bed",
        f"work/annos/grch38/features/ucsc/{DV.ucsc_genomic_super_dups_38}/genomicSuperDups.bed.gz",
        f"work/annos/grch38/features/ucsc/{DV.ucsc_rmsk_38}/rmsk.bed.gz",
        f"work/annos/grch38/features/ucsc/{DV.ucsc_alt_seq_liftover_38}/altSeqLiftOverPsl.bed.gz",
        f"work/annos/grch38/features/ucsc/{DV.ucsc_fix_seq_liftover_38}/fixSeqLiftOverPsl.bed.gz",
        #
        # == output directory ===================================================================
        #
        # -- mehari data
        # ---- frequencies (via annonars)
        f"output/mehari/freqs-grch37-{DV.gnomad_v2}+{DV.gnomad_v2}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/mehari/freqs-grch38-{DV.gnomad_v3}+{DV.gnomad_v2}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        # -- varfish-server-worker data
        # ---- CADD
        f"output/worker/annos/seqvars/cadd-grch37-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/worker/annos/seqvars/cadd-grch38-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- dbSNP
        # f"output/worker/annos/seqvars/dbsnp-grch37-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/dbsnp-grch38-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- dbNSFP
        # f"output/worker/annos/seqvars/dbnsfp-grch37-{DV.dbnsfp}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/dbnsfp-grch38-{DV.dbnsfp}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- dbscSNV
        # f"output/worker/annos/seqvars/dbscsnv-grch37-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/dbscsnv-grch38-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- gnomAD mtDNA
        # f"output/worker/annos/seqvars/gnomad-mtdna-grch37-{DV.gnomad_mtdna}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/gnomad-mtdna-grch38-{DV.gnomad_mtdna}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- gnomAD exomes
        # f"output/worker/annos/seqvars/gnomad-exomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/gnomad-exomes-grch38-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- gnomAD genomes
        # f"output/worker/annos/seqvars/gnomad-genomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/gnomad-genomes-grch38-{DV.gnomad_v3}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- HelixMtDb
        # f"output/worker/annos/seqvars/helixmtdb-grch37-{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/helixmtdb-grch38-{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        # # ---- UCSC conservation
        # f"output/worker/annos/seqvars/cons-grch37-{DV.ucsc_cons_37}+{PV.annonars}/rocksdb/IDENTITY",
        # f"output/worker/annos/seqvars/cons-grch38-{DV.ucsc_cons_38}+{PV.annonars}/rocksdb/IDENTITY",


# ===============================================================================================
# Modular Snakefile Includes
# ===============================================================================================


# -- work directory -----------------------------------------------------------------------------
# Gene-related rules.
include: "rules/work/genes/dbnsfp.smk"
include: "rules/work/genes/ensembl.smk"
include: "rules/work/genes/gnomad.smk"
include: "rules/work/genes/hgnc.smk"
include: "rules/work/genes/ncbi.smk"
# Reference sequence--related rules.
include: "rules/work/reference/human.smk"
# Features (position and not variant specific).
include: "rules/work/annos/features/cons.smk"
include: "rules/work/annos/features/ensembl.smk"
include: "rules/work/annos/features/refseq.smk"
include: "rules/work/annos/features/tads.smk"
include: "rules/work/annos/features/ucsc.smk"
# Sequence variants and annotations.
include: "rules/work/annos/seqvars/cadd.smk"
include: "rules/work/annos/seqvars/dbnsfp.smk"
include: "rules/work/annos/seqvars/dbscsnv.smk"
include: "rules/work/annos/seqvars/dbsnp.smk"
include: "rules/work/annos/seqvars/gnomad_mtdna.smk"
include: "rules/work/annos/seqvars/gnomad_nuclear.smk"
include: "rules/work/annos/seqvars/helix.smk"
# Structural variant related.
include: "rules/work/annos/strucvars/dbvar.smk"
include: "rules/work/annos/strucvars/dgv.smk"
include: "rules/work/annos/strucvars/exac.smk"
include: "rules/work/annos/strucvars/g1k.smk"
include: "rules/work/annos/strucvars/gnomad.smk"
# -- output directory ---------------------------------------------------------------------------
include: "rules/output/mehari/freqs.smk"
include: "rules/output/worker/cadd.smk"
