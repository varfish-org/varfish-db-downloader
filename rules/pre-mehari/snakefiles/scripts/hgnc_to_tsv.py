# scripts/hgnc_to_tsv.py

import csv
import sys

# --- Configuration ---

# The exact header required in the output file, in order.
FINAL_OUTPUT_HEADER = [
    "hgnc_id", "symbol", "name", "locus_group", "locus_type", "status",
    "location", "location_sortable", "alias_symbol", "alias_name",
    "prev_symbol", "prev_name", "gene_family", "gene_family_id",
    "date_approved_reserved", "date_symbol_changed", "date_name_changed",
    "date_modified", "entrez_id", "ensembl_gene_id", "vega_id", "ucsc_id",
    "ena", "refseq_accession", "ccds_id", "uniprot_ids", "pubmed_id",
    "mgd_id", "rgd_id", "lsdb", "cosmic", "omim_id", "mirbase", "homeodb",
    "snornabase", "bioparadigms_slc", "orphanet", "pseudogene_org",
    "horde_id", "merops", "imgt", "iuphar", "kznf_gene_catalog",
    "mamit_trnadb", "cd", "lncrnadb", "enzyme_id", "intermediate_filament_db",
    "rna_central_ids", "lncipedia", "gtrnadb", "agr", "mane_select",
    "ucsc_id_novers"
]

# Dictionary to map old header names from the input file to the new names.
# Only fields that need renaming are listed here.
HEADER_MAP = {
    "gene_group": "gene_family",
    "gene_group_id": "gene_family_id",
    "pseudogene.org": "pseudogene_org",
    "mamit-trnadb": "mamit_trnadb",
    "rna_central_id": "rna_central_ids",
}

# --- Main Logic ---

try:
    # Accessing Snakemake objects
    input_tsv_path = snakemake.input.tsv
    output_tsv_path = snakemake.output.tsv
    release_info_path = snakemake.output.release_info
    version = snakemake.params.version
    genomebuild = snakemake.wildcards.genomebuild

    # --- Process the main TSV file ---
    with open(input_tsv_path, 'r', encoding='utf-8') as infile, \
         open(output_tsv_path, 'w', newline='', encoding='utf-8') as outfile:

        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')

        # Read the original header from the input file
        input_header = next(reader)

        # Create a dictionary from each row for robust column access by name
        # We use a DictReader-like approach manually for more control
        def row_to_dict(header, row):
            return {header[i]: val for i, val in enumerate(row)}

        # Write the new, correct header to the output file
        writer.writerow(FINAL_OUTPUT_HEADER)

        # Process each data row
        for row in reader:
            # Create a dictionary from the row data using the original header
            row_data = row_to_dict(input_header, row)

            # Create a new dictionary for the output data, renaming keys as we go
            processed_data = {HEADER_MAP.get(k, k): v for k, v in row_data.items()}

            # Create the new 'ucsc_id_novers' field
            ucsc_id = processed_data.get("ucsc_id", "")
            if ucsc_id and "." in ucsc_id:
                # This is more robust than removing the last 2 chars
                processed_data["ucsc_id_novers"] = ucsc_id.split('.', 1)[0]
            else:
                # Handle cases with no ID or no version suffix
                processed_data["ucsc_id_novers"] = ucsc_id

            # Build the final output row in the correct order specified by FINAL_OUTPUT_HEADER.
            # Use .get(col, "") to gracefully handle any missing columns (like 'gencc' which is dropped)
            # and to ensure empty fields are written as empty strings.
            output_row = [processed_data.get(col, "") for col in FINAL_OUTPUT_HEADER]
            writer.writerow(output_row)

    # --- Create the release_info file ---
    with open(release_info_path, 'w', encoding='utf-8') as f:
        f.write("table\tversion\tgenomebuild\tnull_value\n")
        f.write(f"Hgnc\t{version}\t{genomebuild}\t\n")

    print(f"Successfully processed {input_tsv_path} into {output_tsv_path}")

except Exception as e:
    # Log any errors for easier debugging
    print(f"Error processing HGNC file: {e}", file=sys.stderr)
    sys.exit(1)