## Rules related to ENSEMBL gene information.


# TODO fix this genombuild mess!
rule genes_ensembl_create_xlink:  # -- create ENSEMBL gene information xlink table
    output:
        tsv="work/genes/ensembl/{ensembl}/ensembl_xlink.tsv",
        tsv_md5="work/genes/ensembl/{ensembl}/ensembl_xlink.tsv.md5",
    params:
        ensembl_archive_url=DV.ensembl_38_archive_url,
    shell:
        r"""
        # Check wehther ensembl version is correct.
        export TMPDIR=$(mktemp -d)

        echo -e "ensembl_gene_id\tensembl_transcript_id\tentrez_id\tgene_symbol" \
        >{output.tsv}

        wget --no-check-certificate \
            -O $TMPDIR/tmp \
            '{params.ensembl_archive_url}/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "entrezgene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>' \

        sort -u $TMPDIR/tmp \
        >> {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """


rule genes_ensembl_download_maps_grch37:  # -- download files for ENST-ENSG mapping (GRCh37)
    output:
        download_txt="work/genes/ensembl/grch37/87/download/knowntoEnsembl.txt.gz",
        download_gtf="work/genes/ensembl/grch37/87/download/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_gtf} \
            'https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'
        wget --no-check-certificate \
            -O {output.download_txt} \
            'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownToEnsembl.txt.gz'
        """


rule genes_ensembl_process_maps_grch37:  # -- process ENST-ENSG mapping (GRCh37)
    input:
        download_txt="work/genes/ensembl/grch37/87/download/knowntoEnsembl.txt.gz",
        download_gtf="work/genes/ensembl/grch37/87/download/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
    output:
        tsv="work/genes/enst_ensg/grch37/87/enst_ensg.tsv",
        tsv_md5="work/genes/enst_ensg/grch37/87/enst_ensg.tsv.md5",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        awk \
            -F $'\t' \
            -f scripts/genes-enst-ensg.awk \
            <(zcat {input.download_gtf}) \
        | sort \
        > $TMPDIR/tmp1.txt

        zcat {input.download_txt} \
        | sed -e 's/\..//g' \
        | sort -k2,2 \
        >> $TMPDIR/tmp2.txt

        echo -e "real_enst\tenst\tensg" > {output.tsv}
        join -t $'\t' -1 2 -2 1 $TMPDIR/tmp2.txt $TMPDIR/tmp1.txt \
        >> {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """


rule genes_ensembl_download_maps_grch38:  # -- download files for ENST-ENSG mapping (GRCh38)
    output:
        download_gtf="work/genes/ensembl/grch38/{ensembl_version}/download/Homo_sapiens.GRCh38.{ensembl_version}.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_gtf} \
            https://ftp.ensembl.org/pub/release-{wildcards.ensembl_version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{wildcards.ensembl_version}.gtf.gz
        """


rule genes_ensembl_process_maps_grch38:  # -- process ENST-ENSG mapping (GRCh38)
    input:
        download_gtf="work/genes/ensembl/grch38/{ensembl_version}/download/Homo_sapiens.GRCh38.{ensembl_version}.gtf.gz",
    output:
        tsv="work/genes/enst_ensg/grch38/{ensembl_version}/enst_ensg.tsv",
        tsv_md5="work/genes/enst_ensg/grch38/{ensembl_version}/enst_ensg.tsv.md5",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        echo -e "enst\tensg" > {output.tsv}
        awk \
            -F $'\t' \
            -f scripts/genes-enst-ensg.awk \
            <(zcat {input.download_gtf}) \
        | sort \
        >> {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """
