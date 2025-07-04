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

GENOMEBUILDS = os.environ.get("GENOMEBUILDS", "grch37,grch38").split(",")

# Define number of threads for resource-intensive tasks
THREADS = int(os.environ.get("THREADS", "32"))

# Define amount of ram for resource-intesnive tasks (in mb)
MEMORY = int(os.environ.get("MEMORY", "196000"))

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

CHROMS = tuple(list(map(str, range(1, 23))) + ["X", "Y", "MT"])

def generate_input_files():
    """
    Helper function to return the input files for the rules.
    """
    genomebuilds = GENOMEBUILDS
    gnomad_versions = {
        "grch37": DV.gnomad_v2,
        "grch38": DV.gnomad_v4,
    }
    gnomad_cnv_versions = {
        "grch37": DV.exac_cnv,
        "grch38": DV.gnomad_cnv4,
    }
    gnomad_sv_versions = {
        "grch37": DV.gnomad_sv,
        "grch38": DV.gnomad_sv4,
    }
    refseq_versions = {
        "grch37": DV.refseq_37,
        "grch38": DV.refseq_38,
    }
    ensembl_versions = {
        "grch37": DV.ensembl_37,
        "grch38": DV.ensembl_38,
    }
    cons_versions = {
        "grch37": DV.ucsc_cons_37,
        "grch38": DV.ucsc_cons_38,
    }
    rmask_versions = {
        "grch37": DV.ucsc_rmsk_37,
        "grch38": DV.ucsc_rmsk_38,
    }
    genomic_super_dups_versions = {
        "grch37": DV.ucsc_genomic_super_dups_37,
        "grch38": DV.ucsc_genomic_super_dups_38,
    }
    lift_over_versions = {
        "grch37": DV.ucsc_alt_seq_liftover_37,
        "grch38": DV.ucsc_alt_seq_liftover_38,
    }
    lift_over_fix_versions = {
        "grch37": DV.ucsc_fix_seq_liftover_37,
        "grch38": DV.ucsc_fix_seq_liftover_38,
    }
    genomebuild_cap = {
        "grch37": "GRCh37",
        "grch38": "GRCh38",
    }

    input_files = []
    if "grch38" in genomebuilds:
        input_files += [
            f"work/genes/clingen/{DV.today}/ClinGen_gene_curation_list_GRCh38.tsv",
            "work/download/genes/alphamissense/1/AlphaMissense_gene_hg38.tsv.gz",
            f"work/download/annos/grch38/seqvars/cadd/{DV.cadd}/whole_genome_SNVs_inclAnno.tsv.gz",
            f"work/download/annos/grch38/seqvars/cadd/{DV.cadd}/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz",
            f"work/annos/grch38/seqvars/helixmtdb/{DV.helixmtdb}/helixmtdb.vcf.gz",
            f"work/annos/grch38/seqvars/gnomad_mtdna/{DV.gnomad_mtdna}/gnomad_mtdna.vcf.gz",
            f"work/download/annos/grch38/seqvars/gnomad_exomes/{DV.gnomad_v4}/.done",
            f"work/download/annos/grch38/seqvars/gnomad_genomes/{DV.gnomad_v4}/.done",
            f"work/annos/grch38/features/cons/{DV.ucsc_cons_38}/ucsc_conservation.tsv",
            f"work/annos/grch38/features/refseq/{DV.refseq_38}/refseq_genes.bed.gz",
            f"output/full/worker/bgdb-gnomad-exomes-cnv-grch38-{DV.gnomad_sv4}+{PV.worker}/bgdb-gnomad-exomes-cnv-grch38.bin",
            f"output/full/worker/bgdb-gnomad-genomes-sv-grch38-{DV.gnomad_sv4}+{PV.worker}/bgdb-gnomad-genomes-sv-grch38.bin",
            # ----- background/population structural variants and annotations thereof
            f"output/full/tracks/track-strucvars-gnomad-sv-grch38-{DV.gnomad_sv4}+{DV.tracks}/gnomad-sv.bed.gz",
            f"output/full/tracks/track-strucvars-gnomad-cnv-grch38-{DV.gnomad_cnv4}+{DV.tracks}/gnomad-cnv.bed.gz",
            # -- pre-mehari
            "output/pre-mehari/GRCh38/gnomAD_constraints/v2.1.1/GnomadConstraints.tsv",
            "output/pre-mehari/GRCh38/gnomAD_constraints/v2.1.1/GnomadConstraints.release_info",
        ]
    if "grch37" in genomebuilds:
        input_files += [
            f"work/genes/enst_ensg/grch37/{DV.ensembl_37}/enst_ensg.tsv",
            f"work/genes/clingen/{DV.today}/ClinGen_gene_curation_list_GRCh37.tsv",
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
            f"work/annos/grch37/features/cons/{DV.ucsc_cons_37}/ucsc_conservation.tsv",
            f"work/annos/grch37/features/refseq/{DV.refseq_37}/refseq_genes.bed.gz",
            f"output/full/worker/bgdb-gnomad-grch37-{DV.gnomad_sv}+{PV.worker}/bgdb-gnomad.bin",
            f"output/full/worker/bgdb-exac-grch37-{DV.exac_cnv}+{PV.worker}/bgdb-exac.bin",
            f"output/full/worker/bgdb-g1k-grch37-{DV.g1k_svs}+{PV.worker}/bgdb-g1k.bin",
            f"output/full/tracks/track-strucvars-exac-grch37-{DV.exac_cnv}+{DV.tracks}/exac.bed.gz",
            f"output/full/tracks/track-strucvars-g1k-grch37-{DV.g1k_svs}+{DV.tracks}/g1k.bed.gz",
            f"output/full/tracks/track-strucvars-gnomad-grch37-{DV.gnomad_sv}+{DV.tracks}/gnomad.bed.gz",
        ]
    for genomebuild in genomebuilds:
        input_files += [
            # -- work
            f"work/annos/{genomebuild}/features/ensembl/{ensembl_versions[genomebuild]}/ensembl_genes.bed.gz",
            # -- mehari data
            # ---- frequencies (via annonars)
            f"output/full/mehari/freqs-{genomebuild}-{gnomad_versions[genomebuild]}+{gnomad_versions[genomebuild]}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/alphamissense-{genomebuild}-{DV.alphamissense}+{PV.annonars}/rocksdb/IDENTITY",
            # -- annonars data
            # ----- sequence variant annotations
            f"output/full/annonars/cadd-{genomebuild}-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/dbsnp-{genomebuild}-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/dbscsnv-{genomebuild}-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/gnomad-mtdna-{genomebuild}-{DV.gnomad_mtdna}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/gnomad-exomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/gnomad-genomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/helixmtdb-{genomebuild}-{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/gnomad-sv-exomes-{genomebuild}-{gnomad_cnv_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/gnomad-sv-genomes-{genomebuild}-{gnomad_sv_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            # ----- sequence annotation
            f"output/full/annonars/functional-{genomebuild}-{refseq_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/full/annonars/regions-{genomebuild}-{DV.today}+{PV.annonars}/rocksdb/IDENTITY",
            # ----- conservation
            f"output/full/annonars/cons-{genomebuild}-{cons_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            # ----- worker
            f"output/full/worker/masked-repeat-{genomebuild}-{rmask_versions[genomebuild]}+{PV.worker}/masked-repeat.bin",
            f"output/full/worker/masked-segdup-{genomebuild}-{genomic_super_dups_versions[genomebuild]}+{PV.worker}/masked-segdup.bin",
            f"output/full/worker/bgdb-dbvar-{genomebuild}-{DV.dbvar}+{PV.worker}/bgdb-dbvar.bin",
            f"output/full/worker/bgdb-dgv-{genomebuild}-{DV.dgv}+{PV.worker}/bgdb-dgv.bin",
            f"output/full/worker/bgdb-dgv-gs-{genomebuild}-{DV.dgv}+{PV.worker}/bgdb-dgv-gs.bin",
            f"output/full/worker/clinvar-strucvars-{genomebuild}-{DV.clinvar_version}+{PV.worker}/clinvar-strucvars.bin",
            f"output/full/worker/patho-mms-{genomebuild}-{DV.patho_mms}+{PV.worker}/patho-mms.bed",
            f"output/full/worker/tads-{genomebuild}-dixon2015/hesc.bed",
            # -- mehari data
            f"output/full/mehari/genes-txs-{genomebuild}-{DV.mehari_tx}/mehari-data-txs-{genomebuild}-ensembl-{DV.mehari_tx}.bin.zst",
            f"output/full/mehari/genes-txs-{genomebuild}-{DV.mehari_tx}/mehari-data-txs-{genomebuild}-refseq-{DV.mehari_tx}.bin.zst",
            f"output/full/mehari/genes-txs-{genomebuild}-{DV.mehari_tx}/mehari-data-txs-{genomebuild}-ensembl-and-refseq-{DV.mehari_tx}.bin.zst",
            f"output/full/tracks/track-strucvars-dbvar-{genomebuild}-{DV.dbvar}+{DV.tracks}/dbvar.bed.gz",
            f"output/full/tracks/track-strucvars-dgv-{genomebuild}-{DV.dgv}+{DV.tracks}/dgv.bed.gz",
            f"output/full/tracks/track-strucvars-dgv-gs-{genomebuild}-{DV.dgv_gs}+{DV.tracks}/dgv-gs.bed.gz",
            # ----- known pathogenic MMS
            f"output/full/tracks/track-strucvars-patho-mms-{genomebuild}-{DV.patho_mms}+{DV.tracks}/patho-mms.bed",
            # ----- problematic regions (rmsk, genomicSuperDups, altSeqLiftOverPsl, fixSeqLiftOverPsl)
            f"output/full/tracks/track-features-ucsc-genomicsuperdups-{genomebuild}-{genomic_super_dups_versions[genomebuild]}+{DV.tracks}/genomicSuperDups.bed.gz",
            f"output/full/tracks/track-features-ucsc-rmsk-{genomebuild}-{rmask_versions[genomebuild]}+{DV.tracks}/rmsk.bed.gz",
            f"output/full/tracks/track-features-ucsc-altseqliftoverpsl-{genomebuild}-{lift_over_versions[genomebuild]}+{DV.tracks}/altSeqLiftOverPsl.bed.gz",
            f"output/full/tracks/track-features-ucsc-fixseqliftoverpsl-{genomebuild}-{lift_over_fix_versions[genomebuild]}+{DV.tracks}/fixSeqLiftOverPsl.bed.gz",
            # ----- tads
            f"output/full/tracks/track-tads-{genomebuild}-dixon2015+{DV.tracks}/hesc.bed",
            # ----- probesets
            f"output/full/tracks/track-enrichment-probesets-targets/agilent-all-exon-v4-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/agilent-all-exon-v5-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/agilent-all-exon-v6-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/agilent-all-exon-v7-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/agilent-all-exon-v8-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/idt-xgen-exome-research-panel-v1-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/idt-xgen-exome-research-panel-v2-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/twist-comprehensive-exome-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/twist-core-exome-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/twist-exome-v2_0-{genomebuild}-{DV.tracks}.bed.gz",
            f"output/full/tracks/track-enrichment-probesets-targets/twist-refseq-exome-{genomebuild}-{DV.tracks}.bed.gz",
            # -- targets
            f"output/reduced-dev/targets/{genomebuild}/refseq/{refseq_versions[genomebuild]}/refseq_target_exons.bed.gz",
            #
            # == development (reduced data) directories =============================================
            #
            # -- annonars
            f"output/reduced-dev/annonars/cadd-{genomebuild}-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/cadd-{genomebuild}-{DV.cadd}+{PV.annonars}/spec.yaml",
            f"output/reduced-dev/annonars/cons-{genomebuild}-{cons_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/cons-{genomebuild}-{cons_versions[genomebuild]}+{PV.annonars}/spec.yaml",
            f"output/reduced-dev/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}a+{PV.annonars}/spec.yaml",
            f"output/reduced-dev/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}c+{PV.annonars}/spec.yaml",
            f"output/reduced-dev/annonars/dbscsnv-{genomebuild}-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/dbscsnv-{genomebuild}-{DV.dbscsnv}+{PV.annonars}/spec.yaml",
            f"output/reduced-dev/annonars/dbsnp-{genomebuild}-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/dbsnp-{genomebuild}-{DV.dbsnp}+{PV.annonars}/spec.yaml",
            f"output/reduced-dev/annonars/gnomad-exomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/gnomad-exomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/spec.yaml",
            f"output/reduced-dev/annonars/gnomad-genomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/annonars/gnomad-genomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/spec.yaml",
            # -- mehari
            f"output/reduced-dev/mehari/freqs-{genomebuild}-{gnomad_versions[genomebuild]}+{gnomad_versions[genomebuild]}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-dev/mehari/freqs-{genomebuild}-{gnomad_versions[genomebuild]}+{gnomad_versions[genomebuild]}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/spec.yaml",
            #
            # == exomes (reduced data) directories ==================================================
            #
            # -- targets
            f"output/reduced-exomes/targets/{genomebuild}/refseq/{refseq_versions[genomebuild]}/refseq_target_exons.bed.gz",
            # -- annonars
            f"output/reduced-exomes/annonars/cadd-{genomebuild}-{DV.cadd}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/cadd-{genomebuild}-{DV.cadd}+{PV.annonars}/spec.yaml",
            f"output/reduced-exomes/annonars/cons-{genomebuild}-{cons_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/cons-{genomebuild}-{cons_versions[genomebuild]}+{PV.annonars}/spec.yaml",
            f"output/reduced-exomes/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}a+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}a+{PV.annonars}/spec.yaml",
            f"output/reduced-exomes/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}c+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/dbnsfp-{genomebuild}-{DV.dbnsfp}c+{PV.annonars}/spec.yaml",
            f"output/reduced-exomes/annonars/dbscsnv-{genomebuild}-{DV.dbscsnv}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/dbscsnv-{genomebuild}-{DV.dbscsnv}+{PV.annonars}/spec.yaml",
            f"output/reduced-exomes/annonars/dbsnp-{genomebuild}-{DV.dbsnp}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/dbsnp-{genomebuild}-{DV.dbsnp}+{PV.annonars}/spec.yaml",
            f"output/reduced-exomes/annonars/gnomad-exomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/gnomad-exomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/spec.yaml",
            f"output/reduced-exomes/annonars/gnomad-genomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/annonars/gnomad-genomes-{genomebuild}-{gnomad_versions[genomebuild]}+{PV.annonars}/spec.yaml",
            # -- mehari
            f"output/reduced-exomes/mehari/freqs-{genomebuild}-{gnomad_versions[genomebuild]}+{gnomad_versions[genomebuild]}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/rocksdb/IDENTITY",
            f"output/reduced-exomes/mehari/freqs-{genomebuild}-{gnomad_versions[genomebuild]}+{gnomad_versions[genomebuild]}+{DV.gnomad_mtdna}+{DV.helixmtdb}+{PV.annonars}/spec.yaml",
            # -- pre-mehari
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/hgnc/{DV.hgnc_quarterly}+{DV.cdot}+{refseq_versions[genomebuild]}/Hgnc.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/hgnc/{DV.hgnc_quarterly}+{DV.cdot}+{refseq_versions[genomebuild]}/Hgnc.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/hgnc/{DV.hgnc_quarterly}+{DV.cdot}+{refseq_versions[genomebuild]}/RefseqToHgnc.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/hgnc/{DV.hgnc_quarterly}+{DV.cdot}+{refseq_versions[genomebuild]}/RefseqToHgnc.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/clinvar/{DV.hgnc_quarterly}+{DV.clinvar_release}/Clinvar.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/clinvar/{DV.hgnc_quarterly}+{DV.clinvar_release}/Clinvar.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/HelixMtDb/{DV.helixmtdb}/HelixMtDb.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/HelixMtDb/{DV.helixmtdb}/HelixMtDb.release_info",
            expand("output/pre-mehari/{genomebuild}/dbSNP/{dbsnp}/Dbsnp.{chrom}.{file_ext}", genomebuild=genomebuild_cap[genomebuild], dbsnp=DV.dbsnp, chrom=CHROMS, file_ext=["tsv", "release_info"]),
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/ensembltogenesymbol/{ensembl_versions[genomebuild]}/EnsemblToGeneSymbol.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/ensembltogenesymbol/{ensembl_versions[genomebuild]}/EnsemblToGeneSymbol.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/MITOMAP/{DV.today}/Mitomap.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/MITOMAP/{DV.today}/Mitomap.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/mtDB/{DV.mtdb}/MtDb.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/mtDB/{DV.mtdb}/MtDb.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/extra_annos/{DV.cadd}/ExtraAnno.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/extra_annos/{DV.cadd}/ExtraAnno.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/extra_annos/{DV.cadd}/ExtraAnnoField.tsv",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/extra_annos/{DV.cadd}/ExtraAnnoField.release_info",
        ]
    # Files independent of genomebuild (or serving both)
    input_files += [
        # -- download
        f"work/download/genes/rcnv/2022/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz",
        f"work/download/genes/ctd/{DV.today}/CTD_diseases.tsv.gz",
        f"work/download/do/{DV.today}/omim-unmapped.csv",
        # NB: dbNSFP is dual reference (for download)
        # NB: dbscSNV is dual reference (for download)
        f"work/download/annos/grch37/seqvars/dbsnp/{DV.dbsnp}/dbsnp.vcf.gz",
        # -- genes
        f"work/genes/dbnsfp/{DV.dbnsfp}/genes.tsv.gz",
        "work/genes/decipher/v3/decipher_hi_prediction.tsv",
        f"work/genes/ensembl/{DV.ensembl_38}/ensembl_xlink.tsv",
        f"work/genes/entrez/{DV.today}/gene_info.jsonl",
        f"work/genes/gnomad/{DV.gnomad_constraints}/gnomad_constraints.tsv",
        f"work/genes/hgnc/{DV.hgnc_quarterly}/hgnc_info.jsonl",
        f"work/genes/omim/{DV.hpo}+{DV.today}+{DV.hgnc_quarterly}/omim_diseases.tsv",
        f"work/genes/orphadata/{DV.orphadata}/orphadata.jsonl",
        f"work/genes/mondo/{DV.today}/mondo.obo",
        "work/genes/rcnv/2022/rcnv_collins_2022.tsv",
        "work/genes/shet/2019/shet_weghorn_2019.tsv",
        "work/genes/domino/20190219/domino.tsv",
        # -- annonars
        f"output/full/annonars/genes-{DV.acmg_sf}+{DV.gnomad_constraints}+{DV.dbnsfp}+{DV.hpo}+{DV.today}+{DV.hgnc_quarterly}+{PV.annonars}/rocksdb/IDENTITY",
        # -- worker data
        f"output/full/worker/genes-xlink-{DV.hgnc_quarterly}+{PV.worker}/genes-xlink.bin",
        f"output/full/worker/acmg-sf-{DV.acmg_sf}+{PV.worker}/acmg_sf.tsv",  # ATTN! source file is placed manually in the data directory
        f"output/full/worker/mim2gene-{DV.today}+{PV.worker}/mim2gene.tsv",
        f"output/full/mehari/genes-xlink-{DV.hgnc_quarterly}/genes-xlink.tsv",
        # -- viguno
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/hp.obo",
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype.hpoa",
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype_to_genes.txt",
        f"output/full/viguno/hpo-{DV.hpo}+{PV.viguno}/hpo.bin",
        # -- viguno (reduced dev)
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/hp.obo",
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype.hpoa",
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype_to_genes.txt",
        f"output/reduced-dev/viguno/hpo-{DV.hpo}+{PV.viguno}/hpo.bin",
        # -- viguno (reduced exomes)
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/hp.obo",
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype.hpoa",
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/phenotype_to_genes.txt",
        f"output/reduced-exomes/viguno/hpo-{DV.hpo}+{PV.viguno}/hpo.bin",
        # -- pre-mehari
        f"output/pre-mehari/noref/hpo/{DV.hpo}/Hpo.tsv",
        f"output/pre-mehari/noref/hpo/{DV.hpo}/Hpo.release_info",
        f"output/pre-mehari/noref/hpo/{DV.hpo}/HpoName.tsv",
        f"output/pre-mehari/noref/hpo/{DV.hpo}/HpoName.release_info",
        # TODO double check the versionen of ensembl
        f"output/pre-mehari/noref/refseqtoensembl/{DV.ensembl_38}/RefseqToEnsembl.tsv",
        f"output/pre-mehari/noref/refseqtoensembl/{DV.ensembl_38}/RefseqToEnsembl.release_info",
        f"output/pre-mehari/noref/acmg/{DV.acmg_sf}/Acmg.tsv",  # ATTN! source file is placed manually in the data directory
        f"output/pre-mehari/noref/acmg/{DV.acmg_sf}/Acmg.release_info",
        f"output/pre-mehari/noref/mim2gene/{DV.today}/Mim2geneMedgen.tsv",
        f"output/pre-mehari/noref/mim2gene/{DV.today}/Mim2geneMedgen.release_info",
        f"output/pre-mehari/noref/refseqtogenesymbol/{DV.today}/RefseqToGeneSymbol.tsv",
        f"output/pre-mehari/noref/refseqtogenesymbol/{DV.today}/RefseqToGeneSymbol.release_info",
    ]
    return input_files


## all -- run all rules
rule all:
    input:
        generate_input_files(),


# ===============================================================================================
# Modular Snakefile Includes
# ===============================================================================================


# -- work directory -----------------------------------------------------------------------------
# Misc rules.
include: "rules/work/misc/hpo.smk"
include: "rules/work/misc/clinvar_download.smk"
include: "rules/work/misc/hgnc_download.smk"
# Gene-related rules.
include: "rules/work/genes/alphamissense.smk"
include: "rules/work/genes/dbnsfp.smk"
include: "rules/work/genes/clingen.smk"
include: "rules/work/genes/conditions.smk"
include: "rules/work/genes/ctd.smk"
include: "rules/work/genes/do.smk"
include: "rules/work/genes/decipher.smk"
include: "rules/work/genes/ensembl.smk"
include: "rules/work/genes/gnomad.smk"
include: "rules/work/genes/gtex.smk"
include: "rules/work/genes/hgnc.smk"
include: "rules/work/genes/mehari_data_tx.smk"
include: "rules/work/genes/ncbi.smk"
include: "rules/work/genes/omim.smk"
include: "rules/work/genes/panelapp.smk"
include: "rules/work/genes/mondo.smk"
include: "rules/work/genes/orphadata.smk"
include: "rules/work/genes/rcnv.smk"
include: "rules/work/genes/shet.smk"
include: "rules/work/genes/domino.smk"
# Reference sequence--related rules.
include: "rules/work/reference/human.smk"
# Features (position and not variant specific).
include: "rules/work/annos/features/cons.smk"
include: "rules/work/annos/features/ensembl.smk"
include: "rules/work/annos/features/refseq.smk"
include: "rules/work/annos/features/tads.smk"
include: "rules/work/annos/features/ucsc.smk"
# Sequence variants and annotations.
include: "rules/work/annos/seqvars/alphamissense.smk"
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
include: "rules/work/annos/strucvars/gnomad_sv4.smk"
# -- output directory ---------------------------------------------------------------------------
# ---- mehari
include: "rules/output/mehari/freqs.smk"
# ---- viguno
include: "rules/output/viguno/hpo.smk"
# ---- annonars
include: "rules/output/annonars/alphamissense.smk"
include: "rules/output/annonars/cadd.smk"
include: "rules/output/annonars/cons.smk"
include: "rules/output/annonars/dbnsfp.smk"
include: "rules/output/annonars/dbscsnv.smk"
include: "rules/output/annonars/dbsnp.smk"
include: "rules/output/annonars/gnomad_exomes.smk"
include: "rules/output/annonars/gnomad_genomes.smk"
include: "rules/output/annonars/gnomad_mtdna.smk"
include: "rules/output/annonars/gnomad_sv.smk"
include: "rules/output/annonars/helix.smk"
include: "rules/output/annonars/genes.smk"
include: "rules/output/annonars/functional.smk"
include: "rules/output/annonars/regions.smk"
# ---- worker
include: "rules/output/worker/patho_mms.smk"
include: "rules/output/worker/clinvar.smk"
include: "rules/output/worker/hgnc.smk"
include: "rules/output/worker/acmg.smk"
include: "rules/output/worker/mim2gene.smk"
include: "rules/output/worker/masked.smk"
include: "rules/output/worker/bgdb.smk"
include: "rules/output/worker/tads.smk"
# ---- tracks
include: "rules/output/tracks/exome_probesets.smk"
# -- reduced output directory (dev/exomes) ------------------------------------------------------
# ---- bed file
include: "rules/reduced/annonars.smk"
include: "rules/reduced/hpo.smk"
include: "rules/reduced/targets.smk"
include: "rules/reduced/mehari.smk"

# pre-mehari rules for tsv output
# -- both releases
include: "rules/pre-mehari/snakefiles/GRCh3_/hgnc.smk"
include: "rules/pre-mehari/snakefiles/GRCh3_/clinvar.smk"
include: "rules/pre-mehari/snakefiles/GRCh3_/dbsnp.smk"
include: "rules/pre-mehari/snakefiles/GRCh3_/helixdb.smk"
include: "rules/pre-mehari/snakefiles/GRCh3_/ensembltogenesymbol.smk"
include: "rules/pre-mehari/snakefiles/GRCh3_/mitomap.smk"
include: "rules/pre-mehari/snakefiles/GRCh3_/mtdb.smk"
include: "rules/pre-mehari/snakefiles/GRCh3_/extra_annos.smk"
# -- grch38
include: "rules/pre-mehari/snakefiles/GRCh38/gnomad.smk"
# -- noref
include: "rules/pre-mehari/snakefiles/noref/hpo.smk"
include: "rules/pre-mehari/snakefiles/noref/refseqtoensembl.smk"
include: "rules/pre-mehari/snakefiles/noref/acmg.smk"
include: "rules/pre-mehari/snakefiles/noref/mim2gene.smk"
include: "rules/pre-mehari/snakefiles/noref/refseqtogenesymbol.smk"