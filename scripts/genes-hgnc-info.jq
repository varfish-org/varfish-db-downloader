# Read in large HGNC full dataset file and create JSONL file from it.

# Extract /response/docs and print each entry to a single line.
#
# Within each line, sort the keys alphabetically.
(
    .response.docs[] |
    to_entries |
    sort_by(.key) |
    from_entries
)
