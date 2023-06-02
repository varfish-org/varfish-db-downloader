"""Declaration of data versions."""

import attrs


@attrs.frozen()
class DataVersions:
    #: String to use for GRCh37 ENSEMBL version.
    ensembl_37: str
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
    #: Version of dbVar.
    dbvar: str


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    ensembl_37="87",
    ensembl="109",
    today="2023-06-01",
    dbnsfp="4.4",
    dbscsnv="1.1",
    cadd="1.6",
    gnomad_constraints="2.1.1",
    gnomad_mtdna="3.1",
    gnomad_v2="2.1.1",
    gnomad_v3="3.1.2",
    dbvar="2023-05-16",
)
