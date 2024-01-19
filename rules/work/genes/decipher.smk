## Rules related to DECIPHER gene information.


rule genes_decipher_hi_download:  # -- download DECIPHER HI predictions
    output:
        bed="work/download/genes/decipher/v3/HI_Predictions_Version3.bed.gz",
        bed_md5="work/download/genes/decipher/v3/HI_Predictions_Version3.bed.gz.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.bed} \
            https://www.deciphergenomics.org/files/downloads/HI_Predictions_Version3.bed.gz

        md5sum {output.bed} > {output.bed_md5}
        """


rule genes_decipher_hi_convert:  # -- convert DECIPHER HI predictions to TSV
    input:
        hgnc=f"output/full/mehari/genes-xlink-{DV.today}/genes-xlink.tsv",
        bed="work/download/genes/decipher/v3/HI_Predictions_Version3.bed.gz",
    output:
        tsv="work/genes/decipher/v3/decipher_hi_prediction.tsv.gz",
        tsv_md5="work/genes/decipher/v3/decipher_hi_prediction.tsv.gz.md5",
    shell:
        r"""
        set -x

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        echo -e "gene_symbol\tp_hi\thi_index" > $TMPDIR/tmp.tsv

        zcat {input.bed} \
        | tail -n +2 \
        | cut -f 4 \
        | tr '|' '\t' \
        | sed -e 's/%$//g' \
        >> $TMPDIR/tmp.tsv

        qsv select -d $'\t' '!omim_ids' {input.hgnc} \
        | tr ',' '\t' \
        > $TMPDIR/tmp3.tsv

        qsv join -d $'\t' \
            gene_symbol $TMPDIR/tmp3.tsv \
            gene_symbol $TMPDIR/tmp.tsv \
        > $TMPDIR/tmp2.csv

        qsv select \
            'hgnc_id,gene_symbol,p_hi,hi_index' \
            $TMPDIR/tmp2.csv \
        | qsv rename 'hgnc_id,hgnc_symbol,p_hi,hi_index' \
        | tr ',' '\t' \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv_md5}
        """
