## Rules to create build annonars genes database..


rule output_annonars_genes:  # -- build annonars genes RocksDB file
    input:
        acmg_sf="data/acmg/{v_acmg_sf}/acmg.tsv",
        gnomad_constraints="work/genes/gnomad/{v_gnomad_constraints}/gnomad_constraints.tsv",
        dbnsfp="work/genes/dbnsfp/{v_dbnsfp}/genes.tsv.gz",
        hgnc="work/genes/hgnc/{date}/hgnc_info.jsonl",
        ncbi="work/genes/entrez/{date}/gene_info.jsonl",
    output:
        rocksdb_identity=(
            "output/annonars/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{date}+{v_annonars}/"
            "rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/annonars/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{date}+{v_annonars}/"
            "spec.yaml"
        ),
    wildcard_constraints:
        v_acmg_sf=RE_VERSION,
        v_gnomad_constraints=RE_VERSION,
        v_dbnsfp=RE_VERSION,
        date=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars gene import \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --path-in-acmg {input.acmg_sf} \
            --path-in-gnomad-constraints {input.gnomad_constraints} \
            --path-in-dbnsfp {input.dbnsfp} \
            --path-in-hgnc {input.hgnc} \
            --path-in-ncbi {input.ncbi}

        varfish-db-downloader tpl \
            --template rules/output/annonars/genes.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.v_acmg_sf}+{wildcards.v_gnomad_constraints}+{wildcards.v_dbnsfp}+{wildcards.date}+{wildcards.v_annonars} \
            --value v_acmg_sf={wildcards.v_acmg_sf} \
            --value v_gnomad_constraints={wildcards.v_gnomad_constraints} \
            --value v_dbnsfp={wildcards.v_dbnsfp} \
            --value date={wildcards.date} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """