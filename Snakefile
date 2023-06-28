# Main Snakefile for the varfish-db-downloader.
#
# This Snakemake workflow allows for downlaoding the background data needed by VarFish.  The data
# is downloaded from public sources and transformed as needed.  Eventually, the data is converted
# into RocksDB databases or protobuf binary files that can be used directly by
# ``varfish-server-worker`` and is used in the backend for filtering and/or exposed to the
# user via a REST API.

from varfish_db_downloader.versions import (
    DATA_VERSIONS as DV,
    PACKAGE_VERSIONS as PV,
    FORCE_TODAY,
    TODAY,
    RUNS_IN_CI,
)

# The prefix to use for all shell commands.
SHELL_PREFIX = "export LC_ALL=C; set -x -euo pipefail;"
# Setup the shell prefix by default.
shell.prefix(SHELL_PREFIX)

# Regular expression for genome release.
RE_GENOME = r"grch(37|38)"
# Regular expression for a "word".
RE_NAME = r"[\w-]+"
# Regular expression for versions.
RE_VERSION = r"\w+(\.\w+)*"
# Regular expression for versions.
RE_VERSION_MULTI = r"\w+(\.\w+)*(\+\w+(\.\w+)*)*"

# Gene symbols for the "dev" subset.
DEV_GENE_SYMBOLS = "|".join(
    [
        "BRCA1",  # commonly tested, well-annotated
        "TTN",  # commonly tested, well-annotated
        "OMA1",  # commonly tested, well-annotated
        "SAMD11",  # very small coordinates => found when CI=true
    ]
)
# Padding to add to exons in dev and exons mode.
EXON_PADDING = 200

# ===============================================================================================
# Test Mode
# ===============================================================================================

# Activate test mode by prepending the path to the "test-mode-bin" directory to the PATH.
if RUNS_IN_CI:
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
        f"output/full/mehari/freqs-grch37-{DV.gnomad_v2}+{DV.gnomad_v2}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/mehari/freqs-grch38-{DV.gnomad_v3}+{DV.gnomad_v2}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        # -- annonars data
        # ----- sequence variant annotations
        f"output/full/annonars/cadd-grch37-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/cadd-grch38-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbsnp-grch37-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbsnp-grch38-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbnsfp-grch37-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbnsfp-grch38-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbnsfp-grch37-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbnsfp-grch38-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbscsnv-grch37-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/dbscsnv-grch38-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/gnomad-mtdna-grch37-{DV.gnomad_mtdna}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/gnomad-mtdna-grch38-{DV.gnomad_mtdna}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/gnomad-exomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/gnomad-exomes-grch38-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/gnomad-genomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/gnomad-genomes-grch38-{DV.gnomad_v3}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/helixmtdb-grch37-{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/helixmtdb-grch38-{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
        # ----- conservation
        f"output/full/annonars/cons-grch37-{DV.ucsc_cons_37}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/full/annonars/cons-grch38-{DV.ucsc_cons_38}+{PV.annonars}/rocksdb/IDENTITY",
        # ----- genes
        f"output/full/annonars/genes-{DV.acmg_sf}+{DV.gnomad_constraints}+{DV.dbnsfp}+{DV.today}+{PV.worker}/rocksdb/IDENTITY",
        # -- worker data
        f"output/full/worker/genes-regions-grch37-{DV.refseq_37}/refseq_genes.bin",
        f"output/full/worker/genes-regions-grch37-{DV.ensembl_37}/ensembl_genes.bin",
        f"output/full/worker/genes-regions-grch38-{DV.refseq_38}/refseq_genes.bin",
        f"output/full/worker/genes-regions-grch38-{DV.ensembl_38}/ensembl_genes.bin",
        f"output/full/worker/genes-xlink-{DV.today}/genes-xlink.bin",
        f"output/full/worker/acmg-sf-{DV.acmg_sf}/acmg_sf.tsv",
        f"output/full/worker/mim2gene/{DV.today}/mim2gene.tsv",
        f"output/full/worker/masked-repeat-grch37-{DV.ucsc_rmsk_37}/masked-repeat.bin",
        f"output/full/worker/masked-repeat-grch38-{DV.ucsc_rmsk_38}/masked-repeat.bin",
        f"output/full/worker/masked-segdup-grch37-{DV.ucsc_genomic_super_dups_37}/masked-segdup.bin",
        f"output/full/worker/masked-segdup-grch38-{DV.ucsc_genomic_super_dups_38}/masked-segdup.bin",
        f"output/full/worker/bgdb-dbvar-grch37-{DV.dbvar}/bgdb-dbvar.bin",
        f"output/full/worker/bgdb-dbvar-grch38-{DV.dbvar}/bgdb-dbvar.bin",
        f"output/full/worker/bgdb-dgv-grch37-{DV.dgv}/bgdb-dgv.bin",
        f"output/full/worker/bgdb-dgv-grch38-{DV.dgv}/bgdb-dgv.bin",
        f"output/full/worker/bgdb-dgv-gs-grch37-{DV.dgv}/bgdb-dgv-gs.bin",
        f"output/full/worker/bgdb-dgv-gs-grch38-{DV.dgv}/bgdb-dgv-gs.bin",
        f"output/full/worker/bgdb-gnomad-grch37-{DV.gnomad_sv}/bgdb-gnomad.bin",
        f"output/full/worker/bgdb-exac-grch37-{DV.exac_cnv}/bgdb-exac.bin",
        f"output/full/worker/bgdb-g1k-grch37-{DV.g1k_svs}/bgdb-g1k.bin",
        # -- mehari data
        f"output/full/mehari/genes-xlink-{DV.today}/genes-xlink.tsv",
        f"output/full/mehari/genes-txs-grch37-{DV.mehari_tx}/mehari-data-txs-grch37-{DV.mehari_tx}.bin.zst",
        f"output/full/mehari/genes-txs-grch38-{DV.mehari_tx}/mehari-data-txs-grch38-{DV.mehari_tx}.bin.zst",
        # ----- HPO
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/hp.obo",
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype.hpoa",
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype_to_genes.txt",
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/hpo.bin",
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/scores-fun-sim-avg-resnik-gene/IDENTITY",
        # ----- background/population structural variants and annotations thereof
        f"output/full/worker/track-strucvars-dbvar-grch37-{DV.dbvar}/dbvar.bed.gz",
        f"output/full/worker/track-strucvars-dbvar-grch38-{DV.dbvar}/dbvar.bed.gz",
        f"output/full/worker/track-strucvars-dgv-grch37-{DV.dgv}/dgv.bed.gz",
        f"output/full/worker/track-strucvars-dgv-grch38-{DV.dgv}/dgv.bed.gz",
        f"output/full/worker/track-strucvars-dgv-gs-grch37-{DV.dgv_gs}/dgv-gs.bed.gz",
        f"output/full/worker/track-strucvars-dgv-gs-grch38-{DV.dgv_gs}/dgv-gs.bed.gz",
        f"output/full/worker/track-strucvars-exac-grch37-{DV.exac_cnv}/exac.bed.gz",
        f"output/full/worker/track-strucvars-g1k-grch37-{DV.g1k_svs}/g1k.bed.gz",
        f"output/full/worker/track-strucvars-gnomad-grch37-{DV.gnomad_sv}/gnomad.bed.gz",
        # ----- ClinVar SVs
        f"output/full/worker/track-strucvars-clinvar-grch37-{DV.clinvar_version}/clinvar_strucvars.bin",
        f"output/full/worker/track-strucvars-clinvar-grch38-{DV.clinvar_version}/clinvar_strucvars.bin",
        # ----- known pathogenic MMS
        f"output/full/worker/track-strucvars-patho-mms-grch37-{DV.patho_mms}/patho-mms.bed",
        f"output/full/worker/track-strucvars-patho-mms-grch38-{DV.patho_mms}/patho-mms.bed",
        # ----- problematic regions (rmsk, genomicSuperDups, altSeqLiftOverPsl, fixSeqLiftOverPsl)
        f"output/full/worker/track-features-ucsc-genomicsuperdups-grch37-{DV.ucsc_genomic_super_dups_37}/genomicSuperDups.bed.gz",
        f"output/full/worker/track-features-ucsc-genomicsuperdups-grch38-{DV.ucsc_genomic_super_dups_38}/genomicSuperDups.bed.gz",
        f"output/full/worker/track-features-ucsc-rmsk-grch37-{DV.ucsc_rmsk_37}/rmsk.bed.gz",
        f"output/full/worker/track-features-ucsc-rmsk-grch38-{DV.ucsc_rmsk_38}/rmsk.bed.gz",
        f"output/full/worker/track-features-ucsc-altseqliftoverpsl-grch37-{DV.ucsc_alt_seq_liftover_37}/altSeqLiftOverPsl.bed.gz",
        f"output/full/worker/track-features-ucsc-altseqliftoverpsl-grch38-{DV.ucsc_alt_seq_liftover_38}/altSeqLiftOverPsl.bed.gz",
        f"output/full/worker/track-features-ucsc-fixseqliftoverpsl-grch37-{DV.ucsc_fix_seq_liftover_37}/fixSeqLiftOverPsl.bed.gz",
        f"output/full/worker/track-features-ucsc-fixseqliftoverpsl-grch38-{DV.ucsc_fix_seq_liftover_38}/fixSeqLiftOverPsl.bed.gz",
        # ----- tads
        "output/full/worker/track-tads-grch37-dixon2015/hesc.bed",
        "output/full/worker/track-tads-grch38-dixon2015/hesc.bed",
        #
        # == development (reduced data) directories =============================================
        #
        # -- targets
        f"output/reduced-dev/targets/grch37/refseq/{DV.refseq_37}/refseq_target_exons.bed.gz",
        f"output/reduced-dev/targets/grch38/refseq/{DV.refseq_38}/refseq_target_exons.bed.gz",
        # -- viguno
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/hp.obo",
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype.hpoa",
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype_to_genes.txt",
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/hpo.bin",
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/scores-fun-sim-avg-resnik-gene/IDENTITY",
        # -- annonars
        f"output/reduced-dev/annonars/cadd-grch37-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/cadd-grch38-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/cons-grch37-{DV.ucsc_cons_37}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/cons-grch38-{DV.ucsc_cons_38}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbsnp-grch37-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbsnp-grch38-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbnsfp-grch37-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbnsfp-grch38-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbnsfp-grch37-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbnsfp-grch38-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbscsnv-grch37-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/dbscsnv-grch38-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/gnomad-exomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/gnomad-exomes-grch38-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/gnomad-genomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-dev/annonars/gnomad-genomes-grch38-{DV.gnomad_v3}+{PV.annonars}/rocksdb/IDENTITY",
        #
        # == exomes (reduced data) directories ==================================================
        #
        # -- targets
        f"output/reduced-exomes/targets/grch37/refseq/{DV.refseq_37}/refseq_target_exons.bed.gz",
        f"output/reduced-exomes/targets/grch38/refseq/{DV.refseq_38}/refseq_target_exons.bed.gz",
        # -- viguno
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/hp.obo",
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype.hpoa",
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype_to_genes.txt",
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/hpo.bin",
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/scores-fun-sim-avg-resnik-gene/IDENTITY",
        # -- annonars
        f"output/reduced-exomes/annonars/cadd-grch37-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/cadd-grch38-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/cons-grch37-{DV.ucsc_cons_37}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/cons-grch38-{DV.ucsc_cons_38}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbsnp-grch37-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbsnp-grch38-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbnsfp-grch37-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbnsfp-grch38-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbnsfp-grch37-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbnsfp-grch38-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbscsnv-grch37-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/dbscsnv-grch38-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/gnomad-exomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/gnomad-exomes-grch38-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/gnomad-genomes-grch37-{DV.gnomad_v2}+{PV.annonars}/rocksdb/IDENTITY",
        f"output/reduced-exomes/annonars/gnomad-genomes-grch38-{DV.gnomad_v3}+{PV.annonars}/rocksdb/IDENTITY",


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
include: "rules/work/annos/strucvars/clinvar.smk"
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
include: "rules/output/annonars/genes.smk"
# ---- worker
include: "rules/output/worker/patho_mms.smk"
include: "rules/output/worker/clinvar.smk"
include: "rules/output/worker/genes_regions.smk"
include: "rules/output/worker/hgnc.smk"
include: "rules/output/worker/acmg.smk"
include: "rules/output/worker/mim2gene.smk"
include: "rules/output/worker/masked.smk"
include: "rules/output/worker/bgdb.smk"
# -- reduced output directory (dev/exomes) ------------------------------------------------------
# ---- bed file
include: "rules/reduced/annonars.smk"
include: "rules/reduced/hpo.smk"
include: "rules/reduced/targets.smk"
