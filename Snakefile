from snakemake import shell
shell.prefix("set -x; set -euo pipefail; ")

rule all:
    input:
        'GRCh37/ncbi_gene/ncbi_gene_info.tsv',
        'GRCh37/ncbi_gene/ncbi_gene_rif.tsv',

rule ncbi_gene_download:
    output:
        release_info_ncbi_gene_info='GRCh37/ncbi_gene/ncbi_gene_info.release_info',
        release_info_ncbi_gene_rif='GRCh37/ncbi_gene/ncbi_gene_rif.release_info',
        ncbi_gene_info_tsv='GRCh37/ncbi_gene/ncbi_gene_info.tsv',
        ncbi_gene_rif_tsv='GRCh37/ncbi_gene/ncbi_gene_rif.tsv',
    shell:
        r"""
        mkdir -p $(dirname {output.ncbi_gene_info_tsv})/download
        pushd $(dirname {output.ncbi_gene_info_tsv})/download

        if [[ ! -e Homo_sapiens.ags.gz.md5 ]]; then
            wget \
                -O Homo_sapiens.ags.gz \
                ftp://ftp.ncbi.nih.gov/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz
            md5sum Homo_sapiens.ags.gz >Homo_sapiens.ags.gz.md5
        fi

        if [[ ! -e linux64.gene2xml ]]; then
            wget \
                -O linux64.gene2xml.gz \
                ftp://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/gene2xml/linux64.gene2xml.gz
            gzip -d -c linux64.gene2xml.gz > linux64.gene2xml
            chmod +x linux64.gene2xml
        fi

        cd ..

        ./download/linux64.gene2xml -b T -c T -i download/Homo_sapiens.ags.gz \
        | tee >(
            python3 ../../tools/refseq_xml_to_tsv.py \
                --input /dev/stdin \
                --output /dev/stdout \
                --info-type generif \
            > ncbi_gene_rif.tsv
            ) \
        | python3 ../../tools/refseq_xml_to_tsv.py \
            --input /dev/stdin \
            --output /dev/stdout \
            --info-type summary \
        > ncbi_gene_info.tsv

        md5sum ncbi_gene_info.tsv >ncbi_gene_info.tsv.md5
        md5sum ncbi_gene_rif.tsv >ncbi_gene_rif.tsv.md5

        popd
        echo -e "table\tversion\nncbi_gene_info\t$(date +%Y-%m-%d)" > {output.release_info_ncbi_gene_info}
        echo -e "table\tversion\nncbi_gene_rif\t$(date +%Y-%m-%d)" > {output.release_info_ncbi_gene_rif}
        """
