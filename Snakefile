# Main Snakefile for the varfish-db-downloader.
#
# This Snakemake workflow allows for downlaoding the background data needed by VarFish.  The data
# is downloaded from public sources and transformed as needed.  Eventually, the data is converted
# into RocksDB databases or protobuf binary files that can be used directly by
# ``varfish-server-worker`` and is used in the backend for filtering and/or exposed to the
# user via a REST API.

from varfish_db_downloader.versions import DATA_VERSIONS as DV, PACKAGE_VERSIONS as PV, TODAY

# The prefix to use for all shell commands.
SHELL_PREFIX = "export LC_ALL=C; set -x -euo pipefail;"
# Setup the shell prefix by default.
shell.prefix(SHELL_PREFIX)

# Regular expression for genome release.
RE_GENOME = r"grch(37|38)"
# Regular expression for versions.
RE_VERSION = r"\w+(\.\w+)*"

# ===============================================================================================
# Test Mode
# ===============================================================================================

import os

# Activate test mode by prepending the path to the "test-mode-bin" directory to the PATH.
if os.environ.get("CI", "false").lower() == "true":
    cwd = os.getcwd()
    old_path = os.environ["PATH"]
    os.environ["PATH"] = f"{cwd}/test-mode-bin:{old_path}"
    RUNS_IN_CI = True
else:
    RUNS_IN_CI = False


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
        # NB: gnomAD-SV GRCh38 was announced end of 2020 but not released yet
        # -- genome browser "features" (position-specific)
        # ---- GRCh37
        f"work/annos/grch37/features/cons/{DV.ucsc_cons_37}/ucsc_conservation.tsv",
        f"work/annos/grch37/features/ensembl/{DV.ensembl_37}/ensembl_genes.bed.gz",
        f"work/annos/grch37/features/refseq/{DV.refseq_37}/refseq_genes.bed.gz",
        # ---- GRCh38
        f"work/annos/grch38/features/cons/{DV.ucsc_cons_38}/ucsc_conservation.tsv",
        f"work/annos/grch38/features/ensembl/{DV.ensembl_38}/ensembl_genes.bed.gz",
        f"work/annos/grch38/features/refseq/{DV.refseq_38}/refseq_genes.bed.gz",
        #
        # == output directory ===================================================================
        #
        # -- mehari data
        # ---- frequencies (via annonars)
        f"output/mehari/freqs-grch37-{DV.gnomad_v2}+{DV.gnomad_v2}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/mehari/freqs-grch38-{DV.gnomad_v3}+{DV.gnomad_v2}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        # ---- annonars data
        f"output/annonars/cadd-grch37-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/cadd-grch38-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbsnp-grch37-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbsnp-grch38-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbnsfp-grch37-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbnsfp-grch38-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbnsfp-grch37-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbnsfp-grch38-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbscsnv-grch37-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/dbscsnv-grch38-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/gnomad-mtdna-grch37-{DV.gnomad_mtdna}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/gnomad-mtdna-grch38-{DV.gnomad_mtdna}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/gnomad-exomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/gnomad-exomes-grch38-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/gnomad-genomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/gnomad-genomes-grch38-{DV.gnomad_v3}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/helixmtdb-grch37-{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/helixmtdb-grch38-{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/cons-grch37-{DV.ucsc_cons_37}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/annonars/cons-grch38-{DV.ucsc_cons_38}+{PV.annonars}/rocksdb/IDENTITY",
        # ----- Genes
        f"output/worker/genes-{DV.acmg_sf}+{DV.gnomad_constraints}+{DV.dbnsfp}+{DV.today}+{PV.worker}/rocksdb/IDENTITY",
        f"output/worker/genes-xlink-{DV.today}/genes-xlink.tsv",
        f"output/worker/genes-txs-grch37-{DV.mehari_tx}/mehari-data-txs-grch37-{DV.mehari_tx}.bin.zst",
        f"output/worker/genes-txs-grch38-{DV.mehari_tx}/mehari-data-txs-grch38-{DV.mehari_tx}.bin.zst",
        # ----- HPO
        f"output/viguno/hpo-{DV.hpo}+{PV.viguno}/hp.obo",
        f"output/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype.hpoa",
        f"output/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype_to_genes.txt",
        f"output/viguno/hpo-{DV.hpo}+{PV.viguno}/hpo.bin",
        # ----- background/population structural variants and annotations thereof
        f"output/worker/annos/strucvars/dbvar-grch37-{DV.dbvar}/dbvar.bed.gz",
        f"output/worker/annos/strucvars/dbvar-grch38-{DV.dbvar}/dbvar.bed.gz",
        f"output/worker/annos/strucvars/dgv-grch37-{DV.dgv}/dgv.bed.gz",
        f"output/worker/annos/strucvars/dgv-grch38-{DV.dgv}/dgv.bed.gz",
        f"output/worker/annos/strucvars/dgv-gs-grch37-{DV.dgv_gs}/dgv-gs.bed.gz",
        f"output/worker/annos/strucvars/dgv-gs-grch38-{DV.dgv_gs}/dgv-gs.bed.gz",
        f"output/worker/annos/strucvars/exac-grch37-{DV.exac_cnv}/exac.bed.gz",
        f"output/worker/annos/strucvars/g1k-grch37-{DV.g1k_svs}/g1k.bed.gz",
        f"output/worker/annos/strucvars/gnomad-grch37-{DV.gnomad_sv}/gnomad.bed.gz",
        # ----- known pathogenic MMS
        f"output/worker/annos/strucvars/patho-mms-grch37-{DV.patho_mms}/patho-mms.bed",
        f"output/worker/annos/strucvars/patho-mms-grch38-{DV.patho_mms}/patho-mms.bed",
        # ----- problematic regions (rmsk, genomicSuperDups, altSeqLiftOverPsl, fixSeqLiftOverPsl)
        f"output/worker/annos/features/ucsc-genomicsuperdups-grch37-{DV.ucsc_genomic_super_dups_37}/genomicSuperDups.bed.gz",
        f"output/worker/annos/features/ucsc-genomicsuperdups-grch38-{DV.ucsc_genomic_super_dups_38}/genomicSuperDups.bed.gz",
        f"output/worker/annos/features/ucsc-rmsk-grch37-{DV.ucsc_rmsk_37}/rmsk.bed.gz",
        f"output/worker/annos/features/ucsc-rmsk-grch38-{DV.ucsc_rmsk_38}/rmsk.bed.gz",
        f"output/worker/annos/features/ucsc-altseqliftoverpsl-grch37-{DV.ucsc_alt_seq_liftover_37}/altSeqLiftOverPsl.bed.gz",
        f"output/worker/annos/features/ucsc-altseqliftoverpsl-grch38-{DV.ucsc_alt_seq_liftover_38}/altSeqLiftOverPsl.bed.gz",
        f"output/worker/annos/features/ucsc-fixseqliftoverpsl-grch37-{DV.ucsc_fix_seq_liftover_37}/fixSeqLiftOverPsl.bed.gz",
        f"output/worker/annos/features/ucsc-fixseqliftoverpsl-grch38-{DV.ucsc_fix_seq_liftover_38}/fixSeqLiftOverPsl.bed.gz",
        # ----- tads
        "output/worker/annos/strucvars/tads-grch37-dixon2015/hesc.bed",
        "output/worker/annos/strucvars/tads-grch38-dixon2015/hesc.bed",


# ===============================================================================================
# Modular Snakefile Includes
# ===============================================================================================


# -- work directory -----------------------------------------------------------------------------
# Misc rules.
include: "rules/work/misc/hpo.smk"
# Gene-related rules.
include: "rules/work/genes/dbnsfp.smk"
include: "rules/work/genes/ensembl.smk"
include: "rules/work/genes/gnomad.smk"
include: "rules/work/genes/hgnc.smk"
include: "rules/work/genes/mehari_data_tx.smk"
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
# ---- mehari
include: "rules/output/mehari/freqs.smk"
# ---- viguno
include: "rules/output/viguno/hpo.smk"
# ------ annonars
include: "rules/output/annonars/cadd.smk"
include: "rules/output/annonars/cons.smk"
include: "rules/output/annonars/dbnsfp.smk"
include: "rules/output/annonars/dbscsnv.smk"
include: "rules/output/annonars/dbsnp.smk"
include: "rules/output/annonars/gnomad_exomes.smk"
include: "rules/output/annonars/gnomad_genomes.smk"
include: "rules/output/annonars/gnomad_mtdna.smk"
include: "rules/output/annonars/helix.smk"
# ---- worker
# ------ global
include: "rules/output/worker/genes.smk"
include: "rules/output/worker/patho_mms.smk"
