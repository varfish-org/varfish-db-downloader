import json
import csv
import gzip

input_header = snakemake.input.header
input_json = snakemake.input.json
output_tsv = snakemake.output.tsv
output_release_info = snakemake.output.release_info

try:
    # Read JSON file
    with gzip.open(input_json, 'r') as json_file:
        data = json.load(json_file)
    
    # Check if 'genes' key exists
    if 'genes' not in data:
        raise ValueError("JSON must contain a 'genes' object")
    
    genes = data['genes']

    with open(input_header, 'r', encoding='utf-8') as header_file:
        header = header_file.read().strip().splitlines()
    
    # Write TSV file
    with open(output_tsv, 'w', newline='', encoding='utf-8') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        
        # Write header
        writer.writerow(header)
        
        # Write data rows
        for entrez_id, gene_info in genes.items():
            hgnc_id = gene_info.get('hgnc', '')
            if hgnc_id:
                writer.writerow([entrez_id, f"HGNC:{hgnc_id}"])   
    
    # Write the content to the output file
    with open(output_release_info, 'w') as f:
        f.write("\t".join(["table", "version", "genomebuild", "null_value"]) + "\n")
        f.write("\t".join(["RefseqToHgnc", snakemake.params.version, snakemake.wildcards.genomebuild, ""]) + "\n")

except json.JSONDecodeError as e:
    print(f"Error: Invalid JSON format - {e}")

except Exception as e:
    print(f"Error: {e}")
