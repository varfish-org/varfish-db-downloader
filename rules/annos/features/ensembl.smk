## Rules related to ENSEMBL features.


rule annos_features_ensembl_gene_regions_grch37_download:  # -- download ENSEMBL gene regions files
    output:
        download_gtf="features/grch37/gene_regions/download/Homo_sapiens.GRCh37.87.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download_gtf} \
            'https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'
        """


rule annos_features_ensembl_gene_regions_grch37_process:  # -- process ENSEMBL gene regions files
    input:
        download_gtf="features/grch37/gene_regions/download/Homo_sapiens.GRCh37.87.gtf.gz",
    output:
        tsv="features/grch37/gene_regions/ensembl.bed.gz",
        tsv_md5="features/grch37/gene_regions/ensembl.bed.gz.md5",
        tsv_tbi="features/grch37/gene_regions/ensembl.bed.gz.tbi",
        tsv_tbi_md5="features/grch37/gene_regions/ensembl.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/features-ensembl-gene-regions.awk \
            <(zcat {input.download_gtf}) \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | sort-bed \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
