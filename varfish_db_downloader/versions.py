"""Declaration of data versions."""

import os
import subprocess
from datetime import datetime

import attrs

#: Value to use for "today".
TODAY = os.environ.get("TODAY", datetime.today().strftime("%Y-%m-%d"))


@attrs.frozen()
class DataVersions:
    """Container with data versions."""

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


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    ensembl_37="87",
    ensembl_38="109",
    ensembl="109",
    today=TODAY,
    dbnsfp="4.4",
    dbscsnv="1.1",
    cadd="1.6",
    gnomad_constraints="2.1.1",
    gnomad_mtdna="3.1",
    gnomad_v2="2.1.1",
    gnomad_v3="3.1.2",
    gnomad_sv="2.1.1",
    dbvar="2023-05-16",
    dgv="2020-02-25",
    dgv_gs="2016-05-15",
    exac_cnv="0.3.1",
    g1k_svs="phase3-v2",
    helixmtdb="20200327",
    ucsc_cons_37="2016-10-07",
    ucsc_cons_38="2019-09-06",
    ucsc_rmsk_37="2020-03-22",
    ucsc_rmsk_38="2022-10-18",
    ucsc_genomic_super_dups_37="2011-10-25",
    ucsc_genomic_super_dups_38="2014-10-19",
    ucsc_alt_seq_liftover_37="2020-03-22",
    ucsc_alt_seq_liftover_38="2022-11-03",
    ucsc_fix_seq_liftover_37="2020-05-24",
    ucsc_fix_seq_liftover_38="2022-11-03",
    refseq_37="105",
    refseq_38="GCF_000001405.40-RS_2023_03",
    dbsnp="b151",
)


@attrs.frozen()
class PackageVersions:
    """Container with package versions."""

    #: Version of ``annona-rs`` executable.
    annonars: str
    #: Version of ``varfish-server-worker`` executable.
    worker: str


def get_version(executable: str) -> str:
    """Return version of ``executable``."""
    tmp: str = subprocess.check_output([executable, "--version"], text=True)
    _, version = tmp.strip().split(" ", 1)
    return version


#: The package versions from environment.
PACKAGE_VERSIONS = PackageVersions(
    annonars=get_version("annonars"),
    worker=get_version("varfish-server-worker"),
)
