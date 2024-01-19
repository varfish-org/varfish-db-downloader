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

        qsv join \
            gene_symbol {input.hgnc} \
            gene_symbol $TMPDIR/tmp.tsv \
        > $TMPDIR/tmp2.tsv

        ( \
            echo -e "hgnc_id\thgnc_symbol\tp_hi\thi_index"; \
            tail -n +2 $TMPDIR/tmp2.tsv \
            | tr ',' '\t' \
            | cut -f 1,6-8 \
            | LC_ALL=C sort \
        ) \
        | gzip -c \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv_md5}
        """
