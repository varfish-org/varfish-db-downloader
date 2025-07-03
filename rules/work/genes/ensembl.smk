## Rules related to ENSEMBL gene information.


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


rule genes_ensembl_download_maps_grch3X:  # -- download files for ENST-ENSG mapping
    output:
        download_gtf="work/genes/ensembl/{genomebuild}/{ensembl_version}/download/ensembl.{cdot_ensembl_gtf}.gz",
    params:
        cdot=DV.cdot
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_gtf} \
            https://github.com/SACGF/cdot/releases/download/data_v{params.cdot}/cdot-{params.cdot}.ensembl.{wildcards.cdot_ensembl_gtf}.gz
        """


rule genes_ensembl_download_maps_grch37_knowngene:  # -- download files for ENST-ENSG mapping (GRCh37)
    output:
        download_txt="work/genes/ensembl/grch37/{ensembl_version}/download/knowntoEnsembl.txt.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_txt} \
            'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownToEnsembl.txt.gz'
        # knownToEnsembl.txt.gz is frozen for hg19
        """


def input_genes_ensembl_process_maps_grch37(wildcards):
    return {
        "download_txt": f"work/genes/ensembl/grch37/{wildcards.ensembl_version}/download/knowntoEnsembl.txt.gz",
        "download_gtf": f"work/genes/ensembl/grch37/{wildcards.ensembl_version}/download/{DV.cdot_ensembl_gtf_37}.gz",
    }


rule genes_ensembl_process_maps_grch37:  # -- process ENST-ENSG mapping (GRCh37)
    input:
        unpack(input_genes_ensembl_process_maps_grch37)
    output:
        tsv="work/genes/enst_ensg/grch37/{ensembl_version}/enst_ensg.tsv",
        tsv_md5="work/genes/enst_ensg/grch37/{ensembl_version}/enst_ensg.tsv.md5",
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


def input_genes_ensembl_process_maps_grch38(wildcards):
    return {
        "download_gtf": f"work/genes/ensembl/grch38/{wildcards.ensembl_version}/download/{DV.cdot_ensembl_gtf_38}.gz",
    }


rule genes_ensembl_process_maps_grch38:  # -- process ENST-ENSG mapping (GRCh38)
    input:
        unpack(input_genes_ensembl_process_maps_grch38)
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
