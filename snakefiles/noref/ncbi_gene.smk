# Extract gene summaries and Gene RIF (gene Reference Into Function) from NCBI/Entrez.


rule noref_ncbi_gene_download:
    output:
        ags="noref/ncbi_gene/{download_date}/download/Homo_sapiens.ags.gz",
        ags_md5="noref/ncbi_gene/{download_date}/download/Homo_sapiens.ags.gz.md5",
        gene2xml="noref/ncbi_gene/{download_date}/download/linux64.gene2xml",
        gene2xml_md5="noref/ncbi_gene/{download_date}/download/linux64.gene2xml.md5",
    shell:
        r"""
        cd $(dirname {output.ags})

        if [[ ! -e Homo_sapiens.ags.gz.md5 ]]; then
            wget \
                -O Homo_sapiens.ags.gz \
                http://ftp.ncbi.nih.gov/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz
            md5sum Homo_sapiens.ags.gz >Homo_sapiens.ags.gz.md5
        fi

        if [[ ! -e linux64.gene2xml.md5 ]]; then
            wget \
                -O linux64.gene2xml.gz \
                http://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/gene2xml/linux64.gene2xml.gz
            gzip -d -c linux64.gene2xml.gz > linux64.gene2xml
            chmod +x linux64.gene2xml
            md5sum linux64.gene2xml > linux64.gene2xml.md5
        fi
        """


rule result_noref_ncbi_gene_process:
    input:
        ags="noref/ncbi_gene/{download_date}/download/Homo_sapiens.ags.gz",
        gene2xml="noref/ncbi_gene/{download_date}/download/linux64.gene2xml",
    output:
        release_info_ncbi_gene_info="noref/ncbi_gene/{download_date}/NcbiGeneInfo.release_info",
        release_info_ncbi_gene_rif="noref/ncbi_gene/{download_date}/NcbiGeneRif.release_info",
        ncbi_gene_info_tsv="noref/ncbi_gene/{download_date}/NcbiGeneInfo.tsv",
        ncbi_gene_rif_tsv="noref/ncbi_gene/{download_date}/NcbiGeneRif.tsv",
    shell:
        r"""
        ./{input.gene2xml} -b T -c T -i {input.ags} \
        | tee >(
            python3 tools/refseq_xml_to_tsv.py \
                --input /dev/stdin \
                --output /dev/stdout \
                --info-type generif \
            > {output.ncbi_gene_rif_tsv}
            ) \
        | python3 tools/refseq_xml_to_tsv.py \
            --input /dev/stdin \
            --output /dev/stdout \
            --info-type summary \
        > {output.ncbi_gene_info_tsv}

        pushd $(dirname {output.ncbi_gene_info_tsv})
        md5sum $(basename {output.ncbi_gene_info_tsv}) > $(basename {output.ncbi_gene_info_tsv}).md5
        md5sum $(basename {output.ncbi_gene_rif_tsv}) > $(basename {output.ncbi_gene_rif_tsv}).md5
        popd

        echo -e "table\tversion\tgenomebuild\tnull_value\nNcbiGeneInfo\t$(date +%Y-%m-%d)\t\t" > {output.release_info_ncbi_gene_info}
        echo -e "table\tversion\tgenomebuild\tnull_value\nNcbiGeneRif\t$(date +%Y-%m-%d)\t\t" > {output.release_info_ncbi_gene_rif}
        """
