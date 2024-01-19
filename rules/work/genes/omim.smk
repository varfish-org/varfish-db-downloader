## Rules related to Omim disease to HGNC ID annotation.


rule genes_omim:  # -- prepare HGNC to OMIM disease mapping
    input:
        mim2gene="work/download/genes/ncbi/{date}/mim2gene_medgen",
        xlink="output/full/mehari/genes-xlink-{date}/genes-xlink.tsv",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
    output:
        tsv="work/genes/omim/{v_hpo}+{date}/omim_diseases.tsv",
    shell:
        """
        set -x

        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        head -n 1 {input.mim2gene} | sed -e 's/ /_/g' -e 's/#//g' \
        > $TMPDIR/mim2gene.tsv
        tail -n +2 {input.mim2gene} \
        | sed -e 's/^/OMIM:/g' \
        >> $TMPDIR/mim2gene.tsv

        grep -v ^# {input.hpoa} \
        > $TMPDIR/phenotype.hpoa

        qsv join -d '\t' \
            entrez_id {input.xlink} \
            GeneID $TMPDIR/mim2gene.tsv \
        | qsv select 'hgnc_id,MIM_number' \
        | qsv rename 'hgnc_id,omim_id' \
        | tr ',' '\t' \
        | qsv sort \
        > $TMPDIR/mim2gene_hgnc.tsv

        echo hgnc_id,omim_id,disease_name \
        > $TMPDIR/output.csv

        qsv join -d '\t' \
            omim_id $TMPDIR/mim2gene_hgnc.tsv \
            database_id $TMPDIR/phenotype.hpoa \
        | qsv select 'hgnc_id,omim_id,disease_name' \
        | tail -n +2 \
        | sort -t , -k1,2V -u \
        >> $TMPDIR/output.csv

        qsv fmt -t '\t' $TMPDIR/output.csv \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv}.md5
        """
