## Rules related to RefSeq features.


rule annos_features_ensembl_gene_regions_grch37_download:  # -- download ENSEMBL gene regions files
    output:
        download_acc="features/grch37/gene_regions/download/chr_accessions_GRCh37.p13",
        download_gtf="features/grch37/gene_regions/download/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {input.download_acc} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/Assembled_chromosomes/chr_accessions_GRCh37.p13'

        wget --no-check-certificate \
            -O {input.download_gtf} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz'
        """


rule annos_features_ensembl_gene_regions_grch37_process:  # -- process ENSEMBL gene regions files
    input:
        download_acc="features/grch37/gene_regions/download/chr_accessions_GRCh37.p13",
        download_gtf="features/grch37/gene_regions/download/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
    output:
        tsv="features/grch37/gene_regions/refseq.bed.gz",
        tsv_md5="features/grch37/gene_regions/refseq.bed.gz.md5",
        tsv_tbi="features/grch37/gene_regions/refseq.bed.gz.tbi",
        tsv_tbi_md5="features/grch37/gene_regions/refseq.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/features-refseq-gene-regions.awk \
            {output.download_acc} \
            <(zcat {output.download_gtf}) \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | sort-bed \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
