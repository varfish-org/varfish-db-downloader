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
        "work/annos/grch37/seqvars/helixmtdb/20200327/helixmtdb.vcf.gz",
        f"work/annos/grch37/seqvars/gnomad_mtdna/{DV.gnomad_mtdna}/gnomad_mtdna.vcf.gz",
        f"work/annos/grch37/seqvars/gnomad_exomes/{DV.gnomad_v2}/.done",
        f"work/annos/grch37/seqvars/gnomad_genomes/{DV.gnomad_v2}/.done",
        # ---- GRCh38
        f"work/download/annos/grch38/seqvars/cadd/{DV.cadd}/whole_genome_SNVs_inclAnno.tsv.gz",
        f"work/download/annos/grch38/seqvars/cadd/{DV.cadd}/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz",
        # NB: dbNSFP is dual reference (for download)
        # NB: dbscSNV is dual reference (for download)
        f"work/download/annos/grch37/seqvars/dbsnp/{DV.dbsnp}/dbsnp.vcf.gz",
        "work/annos/grch38/seqvars/helixmtdb/20200327/helixmtdb.vcf.gz",
        f"work/annos/grch38/seqvars/gnomad_mtdna/{DV.gnomad_mtdna}/gnomad_mtdna.vcf.gz",
        f"work/annos/grch38/seqvars/gnomad_exomes/{DV.gnomad_v2}/.done",
        f"work/annos/grch38/seqvars/gnomad_genomes/{DV.gnomad_v3}/.done",
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
