"""Declaration of data versions."""

import os
import subprocess
from datetime import datetime

import attrs

#: Value to use for "today".
TODAY = os.environ.get("TODAY", datetime.today().strftime("%Y%m%d"))
#: Allow to disable the today check.
FORCE_TODAY = os.environ.get("FORCE_TODAY", "false").lower() == "true"

# Clinvar source:
# https://github.com/varfish-org/clinvar-data-jsonl/releases
# https://github.com/varfish-org/clinvar-data-jsonl/releases/download/clinvar-weekly-20250706/clinvar-data-extract-vars-20250706+0.18.5.tar.gz
# Used by:
# - worker
# - pre-mehari release
CLINVAR_VERSION = "20250706"
# Clinvar-this version
CLINVAR_THIS = "0.18.5"
#: The ClinVar release to use (weekly clinvar release data + clinvar-this).
CLINVAR_RELEASE = os.environ.get("CLINVAR_RELEASE", f"{CLINVAR_VERSION}+{CLINVAR_THIS}")

# The following should be consistent with the mehari-data-tx release:
# https://github.com/varfish-org/mehari-data-tx/blob/main/config/config.yaml#L31
# https://github.com/varfish-org/mehari-data-tx/blob/main/config/config.yaml#L114
#: RefSeq release for GRCh37
REFSEQ_37 = "105.20220307"
#: RefSeq release for GRCh38
REFSEQ_38 = "110"
#:  RefSeq reference for GRCh37, corresponding to REFSEQ_37
REFSEQ_REF_37 = "GCF_000001405.25"
#:  RefSeq reference for GRCh38, corresponding to REFSEQ_38
REFSEQ_REF_38 = "GCF_000001405.40"
#: RefSeq reference build for GRCh37, corresponding to REFSEQ_REF_37
REFSEQ_REF_37_BUILD = "GRCh37.p13"
#: RefSeq reference build for GRCh38, corresponding to REFSEQ_REF_38
REFSEQ_REF_38_BUILD = "GRCh38.p14"
#: Ensembl release for GRCh37
ENSEMBL_37 = "87"
#: Ensembl release for GRCh38
ENSEMBL_38 = "112"

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
    #: URL to use for ENSEMBL archive.
    ensembl_37_archive_url: str
    #: URL to use for ENSEMBL archive.
    ensembl_38_archive_url: str
    #: URL to use for ENSEMBL archive FTP.
    ensembl_37_archive_ftp: str
    #: URL to use for ENSEMBL archive FTP.
    ensembl_38_archive_ftp: str
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
    #: Version of gnomAD v3.
    gnomad_v2: str
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
    #: Base URL for RefSeq releases.
    refseq_base_url: str
    #: Refseq release (GRCh37).
    refseq_37: str
    #: Refseq release (GRCh38).
    refseq_38: str
    #: Refseq reference version (GRCh37).
    refseq_ref_37: str
    #: Refseq reference version (GRCh38).
    refseq_ref_38: str
    #: Refseq assembly (refseq reference + GRCh37 build with patch).
    refseq_ref_37_assembly: str
    #: Refseq assembly (refseq reference + GRCh38 build with patch).
    refseq_ref_38_assembly: str
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
    #: CDOT version.
    cdot: str
    #: HGNC quarterly release date.
    hgnc_quarterly: str
    #: cdot refseq GFF for GRCh37.
    cdot_refseq_gff_json_37: str
    #: cdot refseq GFF for GRCh38.
    cdot_refseq_gff_json_38: str
    #: MtDb version
    mtdb: str
    #: Pre-Mehari release date.
    pre_mehari_release: str


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    alphamissense="1",
    clingen_gene=TODAY,
    clingen_variant=TODAY,
    ensembl_37_archive_url="https://grch37.archive.ensembl.org",  # not possible to tag a specific version but according to them they didn't update essential parts since release 75 (2014)
    ensembl_38_archive_url="https://may2024.archive.ensembl.org",
    ensembl_37_archive_ftp="https://ftp.ensembl.org/pub/grch37",
    ensembl_38_archive_ftp="https://ftp.ensembl.org/pub",
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
    dbsnp="b157",
    acmg_sf="3.1", # ATTN! source file is placed manually in the data directory
    orphadata=TODAY,
    patho_mms="20220730",
    clinvar_release=CLINVAR_RELEASE,
    clinvar_version=CLINVAR_VERSION,
    tracks="0",
    # refseq_37="105",
    # refseq_38="GCF_000001405.40+RS_2023_03",
    # this url should be consisent with where cdot gets its data from
    # https://github.com/SACGF/cdot/blob/main/generate_transcript_data/cdot_transcripts.yaml#L115
    refseq_base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606",
    # The lines/versions below the mehari_tx should be consistent with the mehari-data-tx release:
    # https://github.com/varfish-org/mehari-data-tx/blob/v{mehari_tx}/config/config.yaml
    mehari_tx="0.10.4",  # ATTN! version to be consistent with
    # --- 
    cdot="0.2.27", # line #L30 and others
    hgnc_quarterly="2025-04-01", # line #L239
    cdot_refseq_gff_json_37=f"{REFSEQ_REF_37}_{REFSEQ_REF_37_BUILD}_genomic.{REFSEQ_37}.gff", # line #L114
    cdot_refseq_gff_json_38=f"{REFSEQ_REF_38}_{REFSEQ_REF_38_BUILD}_genomic.{REFSEQ_38}.gff", # line #L31
    hpo="v2025-05-06", # line #L250
    ensembl_37=ENSEMBL_37, # line #L217
    ensembl_38=ENSEMBL_38, # line #L155
    refseq_37=REFSEQ_37,
    refseq_38=REFSEQ_38,
    refseq_ref_37=REFSEQ_REF_37,
    refseq_ref_38=REFSEQ_REF_38,
    refseq_ref_37_assembly=f"{REFSEQ_REF_37}_{REFSEQ_REF_37_BUILD}",
    refseq_ref_38_assembly=f"{REFSEQ_REF_38}_{REFSEQ_REF_38_BUILD}",
    # ---
    mtdb="20210728",  # Was manually downloaded at that date, database doesn't exist anymore
    pre_mehari_release=TODAY,
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
