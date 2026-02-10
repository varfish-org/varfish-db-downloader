#!/usr/bin/env bash
#
# Create a small excerpt of HGNC JSON data for testing.
# Takes the first N entries and maintains valid JSON structure.
#
# Usage: ./create_hgnc_json_excerpt.sh <input.json> <output.json> [num_entries]
#

set -eu

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.json> <output.json> [num_entries]"
    echo "Example: $0 hgnc_complete_set_2026-01-06.json excerpt-data/hgnc_excerpt.json 100"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"
NUM_ENTRIES="${3:-100}"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

echo "Creating excerpt with first $NUM_ENTRIES entries..."

# Use jq to:
# 1. Take only the first N docs from .response.docs
# 2. Update numFound and numFoundExact to reflect the new count
# 3. Keep the responseHeader intact
jq --arg num "$NUM_ENTRIES" '
  .response.docs |= .[:($num | tonumber)] |
  .response.numFound = ($num | tonumber) |
  .response.numFoundExact = true
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Done! Created $OUTPUT_FILE"
echo "File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
echo "Number of docs: $(jq '.response.numFound' "$OUTPUT_FILE")"
