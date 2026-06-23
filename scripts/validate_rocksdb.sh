#!/usr/bin/env bash

# Validate a finished RocksDB directory.
# - Final DB must contain at least one SST file.
# - Final DB must not contain WAL (*.log) files.
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <rocksdb-dir>" >&2
    exit 2
fi

rocksdb_dir="$1"

if [[ ! -d "$rocksdb_dir" ]]; then
    echo "RocksDB directory does not exist: $rocksdb_dir" >&2
    exit 1
fi

if ! find "$rocksdb_dir" -maxdepth 1 -type f -name '*.sst' -print -quit | grep -q .; then
    echo "RocksDB validation failed: no .sst files in $rocksdb_dir" >&2
    exit 1
fi

if find "$rocksdb_dir" -maxdepth 1 -type f -name '*.log' -print -quit | grep -q .; then
    echo "RocksDB validation failed: found WAL .log files in $rocksdb_dir" >&2
    exit 1
fi
