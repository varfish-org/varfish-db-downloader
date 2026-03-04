#!/usr/bin/env python3
"""Generate manifest JSON file for pre-mehari output."""

import json
import os
import sys
from pathlib import Path


def count_lines(tsv_file):
    """Count lines in a TSV file, excluding the header."""
    try:
        with open(tsv_file, 'r') as f:
            # Skip header and count remaining lines
            lines = sum(1 for _ in f) - 1
            return max(0, lines)  # Ensure non-negative
    except Exception as e:
        print(f"Warning: Could not count lines in {tsv_file}: {e}", file=sys.stderr)
        return 0


def parse_release_info(release_info_file):
    """Parse a release_info file to extract table, version, and genomebuild."""
    try:
        with open(release_info_file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                return None
            # Skip header line, read data line
            data_line = lines[1].strip()
            parts = data_line.split('\t')
            if len(parts) >= 3:
                table = parts[0]
                version = parts[1]
                genomebuild = parts[2]
                return [table, version, genomebuild]
    except Exception as e:
        print(f"Warning: Could not parse {release_info_file}: {e}", file=sys.stderr)
    return None


def generate_table_key(table_name, db_group=None):
    """Generate the JSON key for a table.

    Maps table names to their PostgreSQL format:
    - HelixMtDb -> frequencies_helixmtdb
    - Clinvar -> clinvar_clinvar
    - Dbsnp -> dbsnp_dbsnp
    - etc.
    """
    table_lower = table_name.lower()

    # Define mappings for known tables
    mappings = {
        'helixmtdb': 'frequencies_helixmtdb',
        'mitomap': 'frequencies_mitomap',
        'mtdb': 'frequencies_mtdb',
        'clinvar': 'clinvar_clinvar',
        'dbsnp': 'dbsnp_dbsnp',
        'extraanno': 'extra_annos_extraanno',
        'extraannofield': 'extra_annos_extraannofield',
        'hgnc': 'geneinfo_hgnc',
        'knowngeneaa': 'conservation_knowngeneaa',
        'acmg': 'geneinfo_acmg',
        'hpo': 'geneinfo_hpo',
        'hponame': 'geneinfo_hponame',
        'mim2genemedgen': 'geneinfo_mim2genemedgen',
        'refseqtogenesymbol': 'geneinfo_refseqtogenesymbol',
        'ensembltogenesymbol': 'geneinfo_ensembltogenesymbol',
        'refseqtohgnc': 'geneinfo_refseqtohgnc',
    }

    return mappings.get(table_lower, f"unknown_{table_lower}")


def generate_manifest(release_dir, output_file):
    """Generate manifest JSON from a release directory."""
    release_path = Path(release_dir)

    if not release_path.exists():
        print(f"Error: Release directory {release_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    manifest = {}
    import_info = []

    # Find all release_info files
    release_info_files = list(release_path.glob('**/*.release_info'))

    for release_info_file in release_info_files:
        # Parse release_info
        info = parse_release_info(release_info_file)
        if info:
            import_info.append(info)

            # Find corresponding TSV file
            tsv_file = release_info_file.with_suffix('.tsv')
            if tsv_file.exists():
                table_name = info[0]
                line_count = count_lines(tsv_file)

                # Generate manifest key
                manifest_key = generate_table_key(table_name)

                # Handle dbsnp specially - sum all chromosome files
                if table_name.lower() == 'dbsnp':
                    if manifest_key not in manifest:
                        manifest[manifest_key] = 0
                    manifest[manifest_key] += line_count
                else:
                    manifest[manifest_key] = line_count

    # Sort import_info by table name
    import_info.sort(key=lambda x: x[0])

    # Remove duplicates from import_info (for multi-file tables like dbsnp)
    seen = set()
    unique_import_info = []
    for info in import_info:
        key = tuple(info)
        if key not in seen:
            seen.add(key)
            unique_import_info.append(info)

    # Add import_info to manifest
    manifest['importer_importinfo'] = unique_import_info

    # Write manifest to file
    with open(output_file, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"Manifest generated: {output_file}")
    print(f"  Tables: {len(manifest) - 1}")
    print(f"  Import info entries: {len(unique_import_info)}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: generate_manifest.py <release_dir> <output_json>", file=sys.stderr)
        sys.exit(1)

    release_dir = sys.argv[1]
    output_file = sys.argv[2]

    generate_manifest(release_dir, output_file)
