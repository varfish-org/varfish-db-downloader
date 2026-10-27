## Rules related to ENSEMBL gene information.


rule genes_ensembl_download_entrez_grch38:  # -- download ENSEMBL to Entrez mapping (GRCh38)
    output:
        download_tsv="work/genes/ensembl/grch38/{ensembl_version}/download/Homo_sapiens.GRCh38.{ensembl_version}.entrez.tsv.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_tsv} \
            '{DV.ensembl_38_archive_ftp}/release-{wildcards.ensembl_version}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{wildcards.ensembl_version}.entrez.tsv.gz'
        """


rule genes_ensembl_create_xlink:  # -- create ENSEMBL gene information xlink table
    input:
        download_gtf="work/genes/ensembl/grch38/{ensembl}/download/Homo_sapiens.GRCh38.{ensembl}.chr_patch_hapl_scaff.gtf.gz",
        download_entrez="work/genes/ensembl/grch38/{ensembl}/download/Homo_sapiens.GRCh38.{ensembl}.entrez.tsv.gz",
    output:
        tsv="work/genes/ensembl/{ensembl}/ensembl_xlink.tsv",
        tsv_md5="work/genes/ensembl/{ensembl}/ensembl_xlink.tsv.md5",
    shell:
        r"""
        echo -e "ensembl_gene_id\tensembl_transcript_id\tentrez_id\tgene_symbol" \
        >{output.tsv}

        awk \
            -F $'\t' \
            -f scripts/genes-ensembl-xlink.awk \
            <(zcat {input.download_entrez}) \
            <(zcat {input.download_gtf}) \
        | sort -u \
        >> {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """


rule genes_ensembl_download_maps_grch3X:  # -- download files for ENST-ENSG mapping
    wildcard_constraints:
        ensembl_version=r"\d+",
    output:
        download_gtf="work/genes/ensembl/{genomebuild}/{ensembl_version}/download/Homo_sapiens.{genomebuild_cap}.{ensembl_version}.gtf.gz",
    params:
        url_base=lambda wildcards: DV.ensembl_38_archive_ftp
        if wildcards.genomebuild == "grch38"
        else DV.ensembl_37_archive_ftp,
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_gtf} \
            '{params.url_base}/release-{wildcards.ensembl_version}/gtf/homo_sapiens/Homo_sapiens.{wildcards.genomebuild_cap}.{wildcards.ensembl_version}.gtf.gz'
        """


# The GTF including patches and haplotypes has all genes that the BioMart tables had.
rule genes_ensembl_download_gtf_chr_patch_hapl_scaff_grch3X:  # -- download GTF incl. patches and haplotypes
    output:
        download_gtf="work/genes/ensembl/{genomebuild}/{ensembl_version}/download/Homo_sapiens.{genomebuild_cap}.{ensembl_version}.chr_patch_hapl_scaff.gtf.gz",
    params:
        url_base=lambda wildcards: DV.ensembl_38_archive_ftp
        if wildcards.genomebuild == "grch38"
        else DV.ensembl_37_archive_ftp,
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_gtf} \
            '{params.url_base}/release-{wildcards.ensembl_version}/gtf/homo_sapiens/Homo_sapiens.{wildcards.genomebuild_cap}.{wildcards.ensembl_version}.chr_patch_hapl_scaff.gtf.gz'
        """


def input_genes_ensembl_process_maps_grch37(wildcards):
    return {
        "download_txt": f"data/ucsc/hg19/knowntoEnsembl.txt.gz",
        "download_gtf": f"work/genes/ensembl/grch37/{wildcards.ensembl_version}/download/Homo_sapiens.GRCh37.{wildcards.ensembl_version}.gtf.gz",
    }


rule genes_ensembl_process_maps_grch37:  # -- process ENST-ENSG mapping (GRCh37)
    input:
        unpack(input_genes_ensembl_process_maps_grch37),
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
        "download_gtf": f"work/genes/ensembl/grch38/{wildcards.ensembl_version}/download/Homo_sapiens.GRCh38.{wildcards.ensembl_version}.gtf.gz",
    }


rule genes_ensembl_process_maps_grch38:  # -- process ENST-ENSG mapping (GRCh38)
    input:
        unpack(input_genes_ensembl_process_maps_grch38),
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
