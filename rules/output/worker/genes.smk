## Rules to create build worker genes database..


rule output_worker_genes:  # -- build genes protobuf file
    input:
        acmg_sf="data/acmg/{v_acmg_sf}/acmg.tsv",
        gnomad_constraints="work/genes/gnomad/{v_gnomad_constraints}/gnomad_constraints.tsv",
        dbnsfp="work/genes/dbnsfp/{v_dbnsfp}/genes.tsv.gz",
        hgnc="work/genes/hgnc/{date}/hgnc_info.jsonl",
        ncbi="work/genes/entrez/{date}/gene_info.jsonl",
    output:
        rocksdb_identity=(
            "output/worker/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{date}+{v_worker}/"
            "rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/worker/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{date}+{v_worker}/"
            "spec.yaml"
        )
    wildcard_constraints:
        v_acmg_sf=RE_VERSION,
        v_gnomad_constraints=RE_VERSION,
        v_dbnsfp=RE_VERSION,
        date=RE_VERSION,
        v_worker=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db genes build \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --path-in-acmg {input.acmg_sf} \
            --path-in-gnomad-constraints {input.gnomad_constraints} \
            --path-in-dbnsfp {input.dbnsfp} \
            --path-in-hgnc {input.hgnc} \
            --path-in-ncbi {input.ncbi}

        varfish-db-downloader tpl \
            --template rules/output/worker/genes.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.v_acmg_sf}+{wildcards.v_gnomad_constraints}+{wildcards.v_dbnsfp}+{wildcards.date}+{wildcards.v_worker} \
            --value v_acmg_sf={wildcards.v_acmg_sf} \
            --value v_gnomad_constraints={wildcards.v_gnomad_constraints} \
            --value v_dbnsfp={wildcards.v_dbnsfp} \
            --value date={wildcards.date} \
            \
            --value v_worker={wildcards.v_worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
