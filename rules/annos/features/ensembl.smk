## Rules related to ENSEMBL features.


rule annos_features_ensembl_gene_regions_grch37_download:  # -- download ENSEMBL gene regions files
    output:
        gtf="work/download/annos/grch37/ensembl_genes/{version}/Homo_sapiens.GRCh37.{version}.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.gtf} \
            'https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.{version}.gtf.gz'
        """


rule annos_features_ensembl_gene_regions_grch37_process:  # -- process ENSEMBL gene regions files
    input:
        gtf="work/download/annos/grch37/ensembl_genes/{version}/Homo_sapiens.GRCh37.{version}.gtf.gz",
    output:
        tsv="work/annos/grch37/features/ensembl/{version}/ensembl_genes.bed.gz",
        tsv_md5="work/annos/grch37/features/ensembl/{version}/ensembl_genes.bed.gz.md5",
        tsv_tbi="work/annos/grch37/features/ensembl/{version}/ensembl_genes.bed.gz.tbi",
        tsv_tbi_md5="work/annos/grch37/features/ensembl/{version}/ensembl_genes.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/features-ensembl-gene-regions.awk \
            <(zcat {input.gtf}) \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | sort-bed - \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
