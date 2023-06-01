## Rules related to RefSeq features.


rule annos_features_refseq_gene_regions_grch37_download:  # -- download ENSEMBL gene regions files
    output:
        acc="work/download/annos/grch37/refseq/chr_accessions_GRCh37.p13",
        gtf="work/download/annos/grch37/refseq/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.acc} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/Assembled_chromosomes/chr_accessions_GRCh37.p13'

        wget --no-check-certificate \
            -O {output.gtf} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz'
        """


rule annos_features_refseq_gene_regions_grch37_process:  # -- process ENSEMBL gene regions files
    input:
        acc="work/download/annos/grch37/refseq/chr_accessions_GRCh37.p13",
        gtf="work/download/annos/grch37/refseq/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
    output:
        tsv="work/annos/grch37/features/refseq/refseq_genes.bed.gz",
        tsv_md5="work/annos/grch37/features/refseq/refseq_genes.bed.gz.md5",
        tsv_tbi="work/annos/grch37/features/refseq/refseq_genes.bed.gz.tbi",
        tsv_tbi_md5="work/annos/grch37/features/refseq/refseq_genes.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/features-refseq-gene-regions.awk \
            {input.acc} \
            <(zcat {input.gtf}) \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | sort-bed - \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
