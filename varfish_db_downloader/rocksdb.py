"""Inspection of built RocksDB databases."""

import pathlib
import re

#: Matches the table format version in a RocksDB OPTIONS file.
RE_FORMAT_VERSION = re.compile(r"^\s*format_version=(\d+)\s*$", re.MULTILINE)


def max_table_format_version(rocksdb_dir: pathlib.Path) -> int:
    """Return the highest table ``format_version`` in the OPTIONS files of ``rocksdb_dir``."""
    versions = [
        int(version)
        for options in rocksdb_dir.glob("OPTIONS-*")
        for version in RE_FORMAT_VERSION.findall(options.read_text())
    ]
    if not versions:
        raise ValueError(f"no format_version in OPTIONS files of {rocksdb_dir}")
    return max(versions)
