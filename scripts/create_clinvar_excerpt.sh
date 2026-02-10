#!/usr/bin/env bash
#
# Extract ClinVar data archive, take first 1000 lines of each file, and repack.
#
# Usage: ./create_clinvar_excerpt.sh <input.tar.gz> <output_dir>
#

# no pipefail as head is used which always triggers pipefail
set -eu

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input.tar.gz> <output_dir>"
    echo "Example: $0 clinvar-data-extract-vars-20260104+0.18.5.tar.gz excerpt-data/"
    exit 1
fi

INPUT_ARCHIVE="$1"
OUTPUT_DIR="$2"

if [ ! -f "$INPUT_ARCHIVE" ]; then
    echo "Error: Input file '$INPUT_ARCHIVE' not found"
    exit 1
fi

# Create temporary directory
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo "Extracting archive..."
tar xzf "$INPUT_ARCHIVE" -C "$TMPDIR"

# Get the directory name from the archive
ARCHIVE_DIR=$(echo $INPUT_ARCHIVE | sed 's/\.tar\.gz$//')
WORK_DIR="$TMPDIR/$ARCHIVE_DIR"

echo "Processing files in $ARCHIVE_DIR..."

# Process each .jsonl.gz file
for file in "$WORK_DIR"/*.jsonl.gz; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        echo "  Processing $filename (taking first 1000 lines)..."
        
        # Decompress, take first 1000 lines, and recompress
        zcat "$file" | head -n 1000 | gzip > "$WORK_DIR/${filename}.tmp"
        mv "$WORK_DIR/${filename}.tmp" "$file"
    fi
done

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Create the output archive
OUTPUT_ARCHIVE=$(basename "$INPUT_ARCHIVE")
echo "Creating output archive: $OUTPUT_ARCHIVE"

cd "$OUTPUT_DIR"
tar -C "$TMPDIR" -czf "$OUTPUT_ARCHIVE" "$ARCHIVE_DIR"

echo "Done! Output: $OUTPUT_DIR/$OUTPUT_ARCHIVE"
echo "Archive contains:"
tar tzf "$OUTPUT_ARCHIVE"
