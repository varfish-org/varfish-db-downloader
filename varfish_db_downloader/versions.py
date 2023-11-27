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
CLINVAR_RELEASE = os.environ.get("CLINVAR_RELEASE", "2023-0625+0.6.3")
#: The ClinVar version to use (part of the tag and release name).
CLINVAR_VERSION = CLINVAR_RELEASE.replace("-", "").split("+")[0]


#: Wether we run in CI/test mode.
RUNS_IN_CI = os.environ.get("CI", "false").lower() == "true"


@attrs.frozen()
class DataVersions:
    """Container with data versions."""

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
    #: ACMG secondary findings version.
    acmg_sf: str
    #: HPO
    hpo: str
    #: OrphaPacket
    orphapacket: str
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


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    clingen_gene=TODAY,
    clingen_variant=TODAY,
    ensembl_37="87",
    ensembl_38="109",
    ensembl="110",
    today=TODAY,
    dbnsfp="4.4",
    dbscsnv="1.1",
    cadd="1.6",
    gnomad_constraints="2.1.1",
    gnomad_mtdna="3.1",
    gnomad_v2="2.1.1",
    gnomad_v3="3.1.2",
    gnomad_sv="2.1.1",
    gnomad_cnv4="4.0",
    gnomad_sv4="4.0",
    dbvar="20230516",
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
    ucsc_fix_seq_liftover_37="20200524",
    ucsc_fix_seq_liftover_38="20221103",
    refseq_37="105",
    refseq_38="GCF_000001405.40+RS_2023_03",
    dbsnp="b151",
    acmg_sf="3.1",
    hpo="20230606",
    orphapacket="10.1",
    patho_mms="20220730",
    mehari_tx="0.3.0",
    clinvar_release=CLINVAR_RELEASE,
    clinvar_version=CLINVAR_VERSION,
    tracks="0",
    refseq_fe_37="105.20201022",
    refseq_fe_38="110",
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
    _, version = tmp.strip().split(" ", 1)
    return version


def downloader_version() -> str:
    """Return the downloader version."""
    if RUNS_IN_CI:
        return "0.0.0"
    else:
        subprocess.check_output(["git", "describe", "--tags"], text=True).strip()[1:]


#: The package versions from environment.
PACKAGE_VERSIONS = PackageVersions(
    downloader=downloader_version(),
    annonars=get_version("annonars"),
    viguno=get_version("viguno"),
    worker=get_version("varfish-server-worker"),
)
