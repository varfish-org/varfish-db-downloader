"""Declaration of data versions."""

import attrs


@attrs.frozen()
class DataVersions:
    #: Version of dbNSFP.
    dbnsfp: str
    #: Version of dbscSNV.
    dbscsnv: str
    #: Version of CADD.
    cadd: str
    #: Version of gnomAD v2.
    gnomad_v2: str
    #: Version of gnomAD v3.
    gnomad_v3: str


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    dbnsfp="4.4",
    dbscsnv="1.1",
    cadd="1.6",
    gnomad_v2="2.1.1",
    gnomad_v3="3.1.2",
)
