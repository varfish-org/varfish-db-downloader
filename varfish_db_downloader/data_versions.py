"""Declaration of data versions."""

import attrs


@attrs.frozen()
class DataVersions:
    #: Version of dbNSFP.
    dbnsfp: str


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    dbnsfp="4.4",
)
