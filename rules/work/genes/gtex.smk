import json
import csv
import gzip
import sys
import typing

import attrs
import cattrs


@attrs.frozen(auto_attribs=True)
class GtexTissueRecord:
    tissue: str
    tpms: typing.List[float] = attrs.field(factory=list)


@attrs.frozen(auto_attribs=True)
class GtexGeneRecord:
    hgnc_id: str
    ensembl_gene_id: str
    ensembl_gene_version: str
    records: typing.List[GtexTissueRecord]


rule genes_gtex_v8_download:  # -- download GTex v8 gene expression data
    output:
        attributes="work/download/genes/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        attributes_md5="work/download/genes/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.md5",
        genes_tpm="work/download/genes/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
        genes_tpm_md5="work/download/genes/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.attributes} \
            https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

        wget --no-check-certificate \
            -O {output.genes_tpm} \
            https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

        md5sum {output.attributes} > {output.attributes_md5}
        md5sum {output.genes_tpm} > {output.genes_tpm_md5}
        """


rule genes_gtex_v8_map:  # -- map GTex v8 gene files for annonars
    input:
        attributes="work/download/genes/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        genes_tpm="work/download/genes/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
        genes_xlink=f"output/full/mehari/genes-xlink-{DV.today}/genes-xlink.tsv",
    output:
        genes_tpm="work/genes/annonars/gtex_v8/genes_tpm.jsonl.gz",
    run:
        # Load mapping from sample ID to sample tissue details
        smtsd_count = {}
        with open(input.attributes, "rt") as inputf:
            reader = csv.DictReader(inputf, delimiter="\t")
            sampid_to_smtsd = {}
            for row in reader:
                sampid_to_smtsd[row["SAMPID"]] = row["SMTSD"]
                smtsd_count.setdefault(row["SMTSD"], 0)
                smtsd_count[row["SMTSD"]] += 1
        print("Sample counts per tissue:", file=sys.stderr)
        for smtsd, count in sorted(smtsd_count.items(), key=lambda x: x[1], reverse=True):
            print(f"{smtsd}: {count}", file=sys.stderr)
            # Load mapping from ENSEMBL to HGNC gene ID
        with open(input.genes_xlink, "rt") as inputf:
            reader = csv.DictReader(inputf, delimiter="\t")
            ensembl_to_hgnc = {row["ensembl_gene_id"]: row["hgnc_id"] for row in reader}

            # Map GTEx v8 gene expression data to counts JSONL data for annonars
        print("Transmogrifying expression data...", file=sys.stderr)
        with gzip.open(input.genes_tpm, "rt") as inputf, gzip.open(
            output.genes_tpm, "wt"
        ) as outputf:
            for _ in range(2):
                next(inputf)
            reader = csv.DictReader(inputf, delimiter="\t")
            for row in reader:
                ensembl_gene_id, ensembl_gene_version = row["Name"].split(".", 1)
                hgnc_id = ensembl_to_hgnc.get(ensembl_gene_id)
                if hgnc_id is None:
                    print(f"Skipping {ensembl_gene_id}.{ensembl_gene_version}", file=sys.stderr)
                    continue

                tissue_records = {}

                for sampid, tpm in row.items():
                    if not sampid.startswith("GTEX-"):
                        continue
                    smtds = sampid_to_smtsd[sampid]
                    if smtds not in tissue_records:
                        tissue_records[smtds] = GtexTissueRecord(tissue=smtds)
                    tissue_records[smtds].tpms.append(float(tpm))

                gene_record = GtexGeneRecord(
                    hgnc_id=hgnc_id,
                    ensembl_gene_id=ensembl_gene_id,
                    ensembl_gene_version=ensembl_gene_version,
                    records=list(sorted(tissue_records.values(), key=lambda r: r.tissue)),
                )
                for tissue_record in tissue_records.values():
                    tissue_record.tpms.sort()
                print(
                    json.dumps(cattrs.unstructure(gene_record)),
                    file=outputf,
                )
        print("... done transmogrifying GTex data", file=sys.stderr)
