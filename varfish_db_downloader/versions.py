"""Declaration of data versions."""

import os
import subprocess
from datetime import datetime

import attrs

#: Value to use for "today".
TODAY = os.environ.get("TODAY", datetime.today().strftime("%Y%m%d"))
#: Allow to disable the today check.
FORCE_TODAY = os.environ.get("FORCE_TODAY", "false").lower() == "true"
#: The ClinVar release to use (includes annonars version used for building).
CLINVAR_RELEASE = os.environ.get("CLINVAR_RELEASE", "20250410+0.18.5")
#: The ClinVar version to use (part of the tag and release name).
CLINVAR_VERSION = CLINVAR_RELEASE.replace("-", "").split("+")[0]


#: Wether we run in CI/test mode.
RUNS_IN_CI = os.environ.get("CI", "false").lower() == "true"


@attrs.frozen()
class DataVersions:
    """Container with data versions."""

    #: String to use for AlphaMissense version.
    alphamissense: str
    #: String to use for ClinGen gene curation version.
    clingen_gene: str
    #: String to use for ClinGen variant curation version.
    clingen_variant: str
    #: String to use for GRCh37 ENSEMBL version.
    ensembl_37: str
    #: String to use for GRCh38 ENSEMBL version.
    ensembl_38: str
    #: String to use for ENSEMBL version.
    ensembl: str
    #: URL to use for ENSEMBL archive.
    ensembl_37_archive_url: str
    #: URL to use for ENSEMBL archive.
    ensembl_38_archive_url: str
    #: String to use for current date.
    today: str
    #: Version of dbNSFP.
    dbnsfp: str
    #: Version of dbscSNV.
    dbscsnv: str
    #: Version of CADD.
    cadd: str
    #: Version of gnomAD for constraints.
    gnomad_constraints: str
    #: Version of gnomAD mtDNA.
    gnomad_mtdna: str
    #: Version of gnomAD v2.
    gnomad_v2: str
    #: Version of gnomAD v3.
    gnomad_v3: str
    #: Version of gnomAD v4.
    gnomad_v4: str
    #: Version of gnomAD SV.
    gnomad_sv: str
    #: Version of gnomAD CNV v4.
    gnomad_cnv4: str
    #: Version of gnomAD SV v4.
    gnomad_sv4: str
    #: Version of dbVar.
    dbvar: str
    #: Version of DGV.
    dgv: str
    #: Version of DGV Gold Standard.
    dgv_gs: str
    #: ExAC CNVs.
    exac_cnv: str
    #: Thousand Genomes SVs.
    g1k_svs: str
    #: HelixMtDb
    helixmtdb: str
    #: UCSC conservation (GRCh37).
    ucsc_cons_37: str
    #: UCSC conservation (GRCh38).
    ucsc_cons_38: str
    #: UCSC repeat masker (GRCh37).
    ucsc_rmsk_37: str
    #: UCSC repeat masker (GRCh38).
    ucsc_rmsk_38: str
    #: UCSC genomicSuperDups (GRCh37).
    ucsc_genomic_super_dups_37: str
    #: UCSC genomicSuperDups (GRCh38).
    ucsc_genomic_super_dups_38: str
    #: UCSC genome browser altSeqLiftOverPsl (GRCh37).
    ucsc_alt_seq_liftover_37: str
    #: UCSC genome browser altSeqLiftOverPsl (GRCh38).
    ucsc_alt_seq_liftover_38: str
    #: UCSC genome browser fixSeqLiftOverPsl (GRCh37).
    ucsc_fix_seq_liftover_37: str
    #: UCSC genome browser fixSeqLiftOverPsl (GRCh38).
    ucsc_fix_seq_liftover_38: str
    #: RefSeq version (GRCh37).
    refseq_37: str
    #: RefSeq version (GRCh38).
    refseq_38: str
    #: dbSNP version.
    dbsnp: str
    #: dbSNP reference version (GRCh37).
    dbsnp_reference_37: str
    #: dbSNP reference version (GRCh38).
    dbsnp_reference_38: str
    #: dbSNP reference version (GRCh37) for reports.
    dbsnp_reference_37_report: str
    #: dbSNP reference version (GRCh38) for reports.
    dbsnp_reference_38_report: str
    #: dbSNP reference version (GRCh37).
    dbsnp_reference_37_ext: str
    #: dbSNP reference version (GRCh38).
    dbsnp_reference_38_ext: str
    #: ACMG secondary findings version.
    acmg_sf: str
    #: HPO
    hpo: str
    #: Orphadata
    orphadata: str
    #: Pathogenic MMS
    patho_mms: str
    #: Mehari transcript data.
    mehari_tx: str
    #: ClinVar release.
    clinvar_release: str
    #: ClinVar version.
    clinvar_version: str
    #: Marker file for the tracks version.  This allows us to update the
    #: tracks BED files later on.
    tracks: str
    #: RefSeq functional elements for GRCh37.
    refseq_fe_37: str
    #: RefSeq functional elements for GRCh38.
    refseq_fe_38: str
    #: CDOT version.
    cdot: str
    #: HGNC quarterly release date.
    hgnc_quarterly: str
    #: HGNC GFF for GRCh37.
    hgnc_gff_37: str
    #: HGNC GFF for GRCh38.
    hgnc_gff_38: str


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    alphamissense="1",
    clingen_gene=TODAY,
    clingen_variant=TODAY,
    ensembl_37="87",
    ensembl_38="112",  # keep at 112 for consistency with mehari
    ensembl="112",  # keep at 112 for consistency with mehari and ensembl_archive_url
    ensembl_37_archive_url="https://grch37.archive.ensembl.org",  # not possible to tag a specific version but according to them they didn't essenitally update sind 75 (2014)
    ensembl_38_archive_url="https://may2024.archive.ensembl.org",
    today=TODAY,
    dbnsfp="4.5",  # update to 4.9 or 5.1 ?
    dbscsnv="1.1",
    cadd="1.6",
    gnomad_constraints="4.1",
    gnomad_mtdna="3.1",
    gnomad_v2="2.1.1",
    gnomad_v3="3.1.2",
    gnomad_v4="4.1",
    gnomad_sv="2.1.1",
    gnomad_cnv4="4.1",
    gnomad_sv4="4.1",
    dbvar="20231030",
    dgv="20200225",
    dgv_gs="20160515",
    exac_cnv="0.3.1",
    g1k_svs="phase3v2",
    helixmtdb="20200327",
    ucsc_cons_37="20161007",
    ucsc_cons_38="20190906",
    ucsc_rmsk_37="20200322",
    ucsc_rmsk_38="20221018",
    ucsc_genomic_super_dups_37="20111025",
    ucsc_genomic_super_dups_38="20141019",
    ucsc_alt_seq_liftover_37="20200322",
    ucsc_alt_seq_liftover_38="20221103",
    ucsc_fix_seq_liftover_37="20200609",
    ucsc_fix_seq_liftover_38="20221103",
    refseq_37="105",
    refseq_38="GCF_000001405.40+RS_2023_03",
    dbsnp="b157",
    dbsnp_reference_37="GCF_000001405.25",
    dbsnp_reference_38="GCF_000001405.40",
    dbsnp_reference_37_report="GCF_000001405.25-RS_2024_09",
    dbsnp_reference_38_report="GCF_000001405.40-RS_2024_08",
    dbsnp_reference_37_ext="GCF_000001405.25_GRCh37.p13",
    dbsnp_reference_38_ext="GCF_000001405.40_GRCh38.p14",
    acmg_sf="3.1", # ATTN! source file is placed manually in the data directory
    hpo="20250303",
    orphadata=TODAY,
    patho_mms="20220730",
    clinvar_release=CLINVAR_RELEASE,
    clinvar_version=CLINVAR_VERSION,
    tracks="0",
    refseq_fe_37="105.20201022",
    refseq_fe_38="110",
    mehari_tx="0.10.3",
    # The following 4 lines/versions should be consistent with the mehari-data-tx release:
    # https://github.com/varfish-org/mehari-data-tx/blob/v{mehari_tx}/config/config.yaml
    cdot="0.2.27", # 1
    hgnc_quarterly="2025-04-01", # 2
    hgnc_gff_37="GCF_000001405.25_GRCh37.p13_genomic.105.20220307.gff", # 4
    hgnc_gff_38="GCF_000001405.40_GRCh38.p14_genomic.110.gff", # 3
)


@attrs.frozen()
class PackageVersions:
    """Container with package versions."""

    # NB: we do not need the version of mehari as transcripts are built with GitHub Actions.

    #: VarFish DB Downloader.
    downloader: str
    #: Version of ``annonars`` executable.
    annonars: str
    #: Version of ``viguno`` executable.
    viguno: str
    #: Version of ``varfish-server-worker`` executable.
    worker: str


def get_version(executable: str) -> str:
    """Return version of ``executable``."""
    tmp: str = subprocess.check_output([executable, "--version"], text=True)
    if executable == "varfish-db-downloader":
        _, _, version = tmp.strip().split(" ", 2)
    else:
        _, version = tmp.strip().split(" ", 1)
    return version


#: The package versions from environment.
PACKAGE_VERSIONS = PackageVersions(
    annonars=get_version("annonars"),
    viguno=get_version("viguno"),
    worker=get_version("varfish-server-worker"),
    downloader=get_version("varfish-db-downloader"),
)
