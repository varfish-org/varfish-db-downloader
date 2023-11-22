## Rules to create build annonars genes database..


rule output_annonars_genes:  # -- build annonars genes RocksDB file
    input:
        acmg_sf="data/acmg_sf/{v_acmg_sf}/acmg_sf.tsv",
        clingen_37="work/genes/clingen/{date}/ClinGen_gene_curation_list_GRCh37.tsv",
        clingen_38="work/genes/clingen/{date}/ClinGen_gene_curation_list_GRCh38.tsv",
        gnomad_constraints="work/genes/gnomad/{v_gnomad_constraints}/gnomad_constraints.tsv",
        dbnsfp="work/genes/dbnsfp/{v_dbnsfp}/genes.tsv.gz",
        hgnc="work/genes/hgnc/{date}/hgnc_info.jsonl",
        ncbi="work/genes/entrez/{date}/gene_info.jsonl",
        omim="work/genes/omim/{v_hpo}+{date}/omim_diseases.tsv",
        orpha="work/genes/orphapacket/{v_orpha}+{date}/orpha_diseases.tsv",
        rcnv="work/genes/rcnv/2022/rcnv_collins_2022.tsv",
        shet="work/genes/shet/2019/shet_weghorn_2019.tsv",
        gtex="work/genes/annonars/gtex_v8/genes_tpm.jsonl.gz",
        domino="work/genes/domino/20190219/domino.tsv",
    output:
        rocksdb_identity=(
            "output/full/annonars/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{v_hpo}+{v_orpha}+{date}+{v_annonars}/"
            "rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/full/annonars/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{v_hpo}+{v_orpha}+{date}+{v_annonars}/"
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
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        annonars gene import \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --path-in-acmg {input.acmg_sf} \
            --path-in-clingen-37 {input.clingen_37} \
            --path-in-clingen-38 {input.clingen_38} \
            --path-in-gnomad-constraints {input.gnomad_constraints} \
            --path-in-dbnsfp {input.dbnsfp} \
            --path-in-hgnc {input.hgnc} \
            --path-in-omim {input.omim} \
            --path-in-orpha {input.orpha} \
            --path-in-ncbi {input.ncbi} \
            --path-in-rcnv {input.rcnv} \
            --path-in-shet {input.shet} \
            --path-in-gtex {input.gtex} \
            --path-in-domino {input.domino}

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
            --value v_orphapacket={wildcards.v_orpha} \
        > {output.spec_yaml}
        """
