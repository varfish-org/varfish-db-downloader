## Rules to create build annonars genes database..


rule output_annonars_genes:  # -- build annonars genes RocksDB file
    input:
        acmg_sf="data/acmg_sf/{v_acmg_sf}/acmg_sf.tsv",
        clingen_37="work/genes/clingen/{date}/ClinGen_gene_curation_list_GRCh37.tsv",
        clingen_38="work/genes/clingen/{date}/ClinGen_gene_curation_list_GRCh38.tsv",
        gnomad_constraints="work/genes/gnomad/{v_gnomad_constraints}/gnomad_constraints.tsv",
        dbnsfp="work/genes/dbnsfp/{v_dbnsfp}/genes.tsv.gz",
        hgnc="work/genes/hgnc/{hgnc_quarterly_date}/hgnc_info.jsonl",
        ncbi="work/genes/entrez/{date}/gene_info.jsonl",
        omim="work/genes/omim/{v_hpo}+{date}+{hgnc_quarterly_date}/omim_diseases.tsv",
        orpha="work/genes/orphadata/{date}/orpha_diseases.tsv",
        panelapp="work/download/genes/panelapp/{date}/panelapp.jsonl",
        conditions="work/genes/conditions/{v_hpo}+{date}+{hgnc_quarterly_date}/conditions.jsonl",
        rcnv="work/genes/rcnv/2022/rcnv_collins_2022.tsv",
        shet="work/genes/shet/2019/shet_weghorn_2019.tsv",
        gtex="work/genes/annonars/gtex_v8/genes_tpm.jsonl.gz",
        domino="work/genes/domino/20190219/domino.tsv",
        decipher_hi="work/genes/decipher/v3/decipher_hi_prediction.tsv",
    output:
        rocksdb_identity=(
            "output/full/annonars/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{v_hpo}+{date}+{hgnc_quarterly_date}+{v_annonars}/"
            "rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/full/annonars/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{v_hpo}+{date}+{hgnc_quarterly_date}+{v_annonars}/"
            "spec.yaml"
        ),
        manifest=(
            "output/full/annonars/genes-{v_acmg_sf}+{v_gnomad_constraints}+{v_dbnsfp}+{v_hpo}+{date}+{hgnc_quarterly_date}+{v_annonars}/"
            "MANIFEST.txt"
        ),
    wildcard_constraints:
        v_acmg_sf=RE_VERSION,
        v_gnomad_constraints=RE_VERSION,
        v_dbnsfp=RE_VERSION,
        date=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        annonars gene import \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --path-in-acmg {input.acmg_sf} \
            --path-in-clingen-37 {input.clingen_37} \
            --path-in-clingen-38 {input.clingen_38} \
            --path-in-conditions {input.conditions} \
            --path-in-gnomad-constraints {input.gnomad_constraints} \
            --path-in-dbnsfp {input.dbnsfp} \
            --path-in-hgnc {input.hgnc} \
            --path-in-omim {input.omim} \
            --path-in-orpha {input.orpha} \
            --path-in-panelapp {input.panelapp} \
            --path-in-ncbi {input.ncbi} \
            --path-in-rcnv {input.rcnv} \
            --path-in-shet {input.shet} \
            --path-in-gtex {input.gtex} \
            --path-in-domino {input.domino} \
            --path-in-decipher-hi {input.decipher_hi}

        varfish-db-downloader tpl \
            --template rules/output/annonars/genes.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.v_acmg_sf}+{wildcards.v_gnomad_constraints}+{wildcards.v_dbnsfp}+{wildcards.date}+{wildcards.hgnc_quarterly_date}+{wildcards.v_annonars} \
            --value v_acmg_sf={wildcards.v_acmg_sf} \
            --value v_gnomad_constraints={wildcards.v_gnomad_constraints} \
            --value v_dbnsfp={wildcards.v_dbnsfp} \
            --value date={wildcards.date} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """
