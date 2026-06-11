#!/usr/bin/env bash

cleanup_partial_rocksdb() {
    local output_rocksdb="$1"
    local rc=$?

    if [[ -d "$output_rocksdb" ]]; then
        echo "Cleaning incomplete RocksDB at $output_rocksdb" >&2
        rm -rf "$output_rocksdb"
    fi

    return "$rc"
}