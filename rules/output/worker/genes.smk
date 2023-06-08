## Rules to create build worker genes database..

import os


rule output_worker_genes:  # -- build genes protobuf file
    input:
        acmg_sf="data/acmg/{v_acmg_sf}/acmg.tsv",
        gnomad_constraints="work/genes/gnomad/{v_gnomad_constraints}/gnomad_constraints.tsv",
        dbnsfp="work/genes/dbnsfp/{v_dbnsfp}/genes.tsv.gz",
        hgnc="work/genes/hgnc/{date}/hgnc_info.jsonl",
        ncbi="work/genes/entrez/{date}/gene_info.jsonl",
    output:
        "output/worker/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{date}+{v_worker}/rocksdb/IDENTITY",
    wildcard_constraints:
        v_acmg_sf=RE_VERSION,
        v_gnomad_constraints=RE_VERSION,
        v_dbnsfp=RE_VERSION,
        date=RE_VERSION,
        v_worker=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db genes build \
            --path-out-rocksdb $(dirname {output}) \
            --path-in-acmg {input.acmg_sf} \
            --path-in-gnomad-constraints {input.gnomad_constraints} \
            --path-in-dbnsfp {input.dbnsfp} \
            --path-in-hgnc {input.hgnc} \
            --path-in-ncbi {input.ncbi}
        """
