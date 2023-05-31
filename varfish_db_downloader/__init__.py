"""VarFish DB Downloader Python Utility Code"""

import re
import warnings

import pkg_resources

_is_released_version = False
try:
    __version__ = pkg_resources.get_distribution("varfish_db_downloader").version
    if re.match(r"^\d+\.\d+\.\d+$", __version__) is not None:
        _is_released_version = True
except pkg_resources.DistributionNotFound:
    warnings.warn("can't get __version__ because %s package isn't installed" % __package__, Warning)
    __version__ = None
