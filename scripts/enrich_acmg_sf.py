#!/usr/bin/env python3
"""
Enrich ACMG SF raw data with gene identifiers from HGNC REST API.

This script takes the raw ACMG SF data and enriches it with:
- HGNC ID
- ENSEMBL gene ID
- NCBI gene ID

The enrichment is done via the HGNC REST API.
"""

import argparse
import csv
import sys
import time
from typing import Dict, Optional

import requests


def fetch_gene_info_from_hgnc(gene_symbol: str) -> Optional[Dict[str, str]]:
    """
    Fetch gene information from HGNC REST API.

    Args:
        gene_symbol: The HGNC gene symbol

    Returns:
        Dictionary with hgnc_id, ensembl_gene_id, ncbi_gene_id, or None if not found
    """
    url = f"https://rest.genenames.org/fetch/symbol/{gene_symbol}"
    headers = {"Accept": "application/json"}

    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()

        data = response.json()

        if "response" in data and "docs" in data["response"] and len(data["response"]["docs"]) > 0:
            doc = data["response"]["docs"][0]

            return {
                "hgnc_id": doc.get("hgnc_id", ""),
                "ensembl_gene_id": doc.get("ensembl_gene_id", ""),
                "ncbi_gene_id": doc.get("entrez_id", ""),
                "gene_symbol": doc.get("symbol", gene_symbol),
            }
        else:
            print(f"Warning: No data found for gene symbol: {gene_symbol}", file=sys.stderr)
            return None

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for {gene_symbol}: {e}", file=sys.stderr)
        return None


def normalize_phenotype_category(category: str) -> str:
    """Normalize phenotype category names."""
    # Map "Other" to "Miscellaneous" and "Metabolic" variations
    category = category.strip()
    if category == "Other":
        return "Miscellaneous"
    # Handle special cases like "Cardiovascular  Metabolic" (double space)
    category = " ".join(category.split())
    return category


def normalize_sf_version(version: str) -> str:
    """Normalize SF list version to always include decimal point."""
    version = version.strip()
    # If the version is just an integer (like "1", "2", "3"), add ".0"
    if version and "." not in version:
        return f"{version}.0"
    return version


def process_acmg_sf_file(input_file: str, output_file: str, delay: float = 0.2):
    """
    Process the raw ACMG SF file and enrich it with gene identifiers.

    Args:
        input_file: Path to the raw TSV file
        output_file: Path to the output enriched TSV file
        delay: Delay in seconds between API calls (to be respectful to the API)
    """
    output_columns = [
        "hgnc_id",
        "ensembl_gene_id",
        "ncbi_gene_id",
        "gene_symbol",
        "mim_gene_id",
        "disease_phenotype",
        "disorder_mim",
        "phenotype_category",
        "inheritance",
        "sf_list_version",
        "variants_to_report",
    ]

    # Cache for gene information to avoid redundant API calls
    gene_cache = {}

    enriched_rows = []

    print(f"Reading input file: {input_file}")
    with open(input_file, "r", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            gene_symbol = row["Gene"].strip()

            # Fetch gene info from cache or API
            if gene_symbol not in gene_cache:
                print(f"Fetching data for: {gene_symbol}")
                gene_info = fetch_gene_info_from_hgnc(gene_symbol)
                if gene_info:
                    gene_cache[gene_symbol] = gene_info
                    time.sleep(delay)  # Be respectful to the API
                else:
                    # Create empty entry if not found
                    gene_cache[gene_symbol] = {
                        "hgnc_id": "",
                        "ensembl_gene_id": "",
                        "ncbi_gene_id": "",
                        "gene_symbol": gene_symbol,
                    }

            gene_info = gene_cache[gene_symbol]

            # Build enriched row
            enriched_row = {
                "hgnc_id": gene_info["hgnc_id"],
                "ensembl_gene_id": gene_info["ensembl_gene_id"],
                "ncbi_gene_id": gene_info["ncbi_gene_id"],
                "gene_symbol": gene_info["gene_symbol"],
                "mim_gene_id": row["Gene MIM"],
                "disease_phenotype": row.get("Disease/Phentyope", row.get("Disease/Phenotype", "")),
                "disorder_mim": row["Disorder MIM"],
                "phenotype_category": normalize_phenotype_category(row["Phenotype Category"]),
                "inheritance": row["Inheritance"],
                "sf_list_version": normalize_sf_version(row["SF List Version"]),
                "variants_to_report": row["Variants to report"],
            }

            enriched_rows.append(enriched_row)

    print(f"Writing output file: {output_file}")
    with open(output_file, "w", encoding="utf-8", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=output_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(enriched_rows)

    print(f"Done! Processed {len(enriched_rows)} rows.")
    print(f"Unique genes fetched: {len(gene_cache)}")


def main():
    parser = argparse.ArgumentParser(
        description="Enrich ACMG SF raw data with gene identifiers from HGNC"
    )
    parser.add_argument(
        "input",
        help="Path to the raw ACMG SF TSV file (e.g., data/acmg_sf/3.2/acmg_sf_raw.tsv)",
    )
    parser.add_argument(
        "output",
        help="Path to the output enriched TSV file (e.g., data/acmg_sf/3.2/acmg_sf.tsv)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.2,
        help="Delay in seconds between API calls (default: 0.2)",
    )

    args = parser.parse_args()

    process_acmg_sf_file(args.input, args.output, args.delay)


if __name__ == "__main__":
    main()
