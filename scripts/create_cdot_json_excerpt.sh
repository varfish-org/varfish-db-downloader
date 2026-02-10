#!/usr/bin/env bash
#
# Create a small excerpt of cdot JSON data for testing.
# Takes the first N genes and transcripts while maintaining valid JSON structure.
#
# Usage: ./create_cdot_json_excerpt.sh <input.json.gz> <output.json.gz> [num_entries]
#

set -eu

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.json.gz> <output.json.gz> [num_entries]"
    echo "Example: $0 cdot-0.2.32.Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz excerpt.json.gz 50"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"
NUM_ENTRIES="${3:-50}"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

echo "Creating excerpt with first $NUM_ENTRIES genes and transcripts..."

# Use jq to:
# 1. Take only the first N entries from .genes (as object)
# 2. Take only the first N entries from .transcripts (as object)
# 3. Keep all other metadata intact
zcat "$INPUT_FILE" | jq --arg num "$NUM_ENTRIES" '
  .genes = (.genes | to_entries | .[:($num | tonumber)] | from_entries) |
  .transcripts = (.transcripts | to_entries | .[:($num | tonumber)] | from_entries)
' | gzip > "$OUTPUT_FILE"

echo "Done! Created $OUTPUT_FILE"
echo "File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
echo "Number of genes: $(zcat "$OUTPUT_FILE" | jq '.genes | length')"
echo "Number of transcripts: $(zcat "$OUTPUT_FILE" | jq '.transcripts | length')"
