#!/usr/bin/env python
"""Helper script to extract gene-disease association from orphapacket."""

import csv
import json
import pathlib
import sys


def main():
    symbol_to_hgnc = {}
    with open(sys.argv[1], "rt") as inputf:
        reader = csv.DictReader(inputf, delimiter="\t")
        for record in reader:
            symbol_to_hgnc[record["gene_symbol"]] = record["hgnc_id"]

    print(f"# xlink entries: {len(symbol_to_hgnc)}", file=sys.stderr)

    base_path = pathlib.Path(sys.argv[2])
    print("\t".join(["hgnc_id", "orpha_id", "disease_name"]))
    for json_path in sorted(base_path.glob("*.json")):
        with json_path.open("rt") as inputf:
            data = json.load(inputf)
            elem_top = data["Orphapacket"]
            if elem_top.get("DisorderType", {}).get("value") != "Disease":
                continue  # skip categories
            disease_name = elem_top["Label"]
            orpha_id = elem_top["PURL"].replace("http://www.orpha.net/ORDO/Orphanet_", "ORPHA:")
            elem_genes = elem_top.get("Genes", [])
            for elem_gene in elem_genes:
                gene_symbol = elem_gene["Gene"]["Symbol"]
                hgnc_id = symbol_to_hgnc.get(gene_symbol)
                if hgnc_id:  # skip if no HGNC ID exists, maybe withdrawn?
                    print("\t".join(map(str, [hgnc_id, orpha_id, disease_name])))


if __name__ == "__main__":
    main()
