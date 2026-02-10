#!/usr/bin/env bash
#
# Create a dummy .bin.zst file for testing/CI purposes.
# Creates a small binary file and compresses it with zstd.
#
# Usage: ./create_dummy_bin_zst.sh <output_file.bin.zst>
#

set -eu

if [ $# -ne 1 ]; then
    echo "Usage: $0 <output_file.bin.zst>"
    echo "Example: $0 dummy.bin.zst"
    exit 1
fi

OUTPUT_FILE="$1"

if [[ ! "$OUTPUT_FILE" =~ \.bin\.zst$ ]]; then
    echo "Error: Output file must have .bin.zst extension"
    exit 1
fi

# Create a minimal binary file with some dummy data
TMPFILE=$(mktemp)
trap 'rm -f "$TMPFILE"' EXIT

# Write a simple binary structure (e.g., a header with magic bytes and some data)
# This creates a 1KB file with dummy content
printf "DUMMY_BIN_FILE\x00\x00\x00" > "$TMPFILE"
dd if=/dev/zero bs=1024 count=1 2>/dev/null >> "$TMPFILE"

# Compress with zstd
zstd -q -f "$TMPFILE" -o "$OUTPUT_FILE"

echo "Created dummy file: $OUTPUT_FILE"
echo "File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
