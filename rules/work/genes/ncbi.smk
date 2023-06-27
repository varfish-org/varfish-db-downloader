## Rules related to gene data from NCBI.


rule genes_ncbi_download_mim2gene:  # -- download NCBI MedGen mim2gene
    output:
        download="work/download/genes/ncbi/{date}/mim2gene_medgen",
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        wget --no-check-certificate \
            -O {output.download} \
            https://ftp.ncbi.nih.gov/gene/DATA/mim2gene_medgen
        """


rule genes_ncbi_entrez_download:  # -- download NCBI Entrez files
    output:
        ags="work/download/genes/ncbi/{date}/Homo_sapiens.ags.gz",
        ags_md5="work/download/genes/ncbi/{date}/Homo_sapiens.ags.gz.md5",
        gene2xml="work/download/genes/ncbi/{date}/linux64.gene2xml",
        gene2xml_md5="work/download/genes/ncbi/{date}/linux64.gene2xml.md5",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        wget --no-check-certificate \
            -O $(dirname {output.ags})/Homo_sapiens.ags.gz \
            https://ftp.ncbi.nih.gov/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz

        wget --no-check-certificate \
            -O $TMPDIR/linux64.gene2xml.gz \
            https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/gene2xml/linux64.gene2xml.gz

        gzip -d -c $TMPDIR/linux64.gene2xml.gz \
        > {output.gene2xml}
        chmod u+x {output.gene2xml}

        md5sum {output.ags} >{output.ags_md5}
        md5sum {output.gene2xml} >{output.gene2xml_md5}
        """


rule genes_ncbi_entrez_process:  # -- process NCBI Entrez files
    input:
        ags="work/download/genes/ncbi/{date}/Homo_sapiens.ags.gz",
        gene2xml="work/download/genes/ncbi/{date}/linux64.gene2xml",
    output:
        jsonl="work/genes/entrez/{date}/gene_info.jsonl",
        jsonl_md5="work/genes/entrez/{date}/gene_info.jsonl.md5",
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        ./{input.gene2xml} -b T -c T -i {input.ags} \
        | python3 scripts/refseq_xml_to_json.py \
        > {output.jsonl}

        md5sum {output.jsonl} >{output.jsonl_md5}
        """
