## Rules related to RefSeq features.


rule annos_features_refseq_gene_regions_download_grch37:  # -- download ENSEMBL gene regions files (GRCh37)
    output:
        acc="work/download/annos/grch37/refseq/{version}/chr_accessions_GRCh37.p13",
        gtf="work/download/annos/grch37/refseq/{version}/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.acc} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.{wildcards.version}/Assembled_chromosomes/chr_accessions_GRCh37.p13'

        wget --no-check-certificate \
            -O {output.gtf} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/{wildcards.version}.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz'
        """


rule annos_features_refseq_gene_regions_download_grch38:  # -- download ENSEMBL gene regions files (GRCh38)
    output:
        report="work/download/annos/grch38/refseq/{version}/GCF_000001405.40_GRCh38.p14_assembly_report.txt",
        acc="work/download/annos/grch38/refseq/{version}/chr_accessions_GRCh38.p14",
        gtf="work/download/annos/grch38/refseq/{version}/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
    shell:
        r"""
        prefix=https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14

        wget --no-check-certificate \
            -O {output.report} \
            "$prefix/GCF_000001405.40_GRCh38.p14_assembly_report.txt"

        echo -e "#Chromosome\tRefSeq Accession.version\tRefSeq\tgi\tGenBank Accession.version\tGenBank gi" \
        > {output.acc}
        cat {output.report} \
        | tr -d '\r' \
        | awk -F $'\t' 'BEGIN {{ OFS=FS }} ($1 !~ /^#/) {{ print $10, $7, $9, $5, "." }}' \
        >> {output.acc}

        wget --no-check-certificate \
            -O {output.gtf} \
            "$prefix/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
        """


def input_annos_features_refseq_gene_regions_process(wildcards):
    """Input function for ``rule annos_features_refseq_gene_regions_process``."""
    if wildcards.genome_build == "grch37":
        return {
            "acc": f"work/download/annos/grch37/refseq/{wildcards.version}/chr_accessions_GRCh37.p13",
            "gtf": f"work/download/annos/grch37/refseq/{wildcards.version}/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
        }
    else:
        return {
            "acc": f"work/download/annos/grch38/refseq/{wildcards.version}/chr_accessions_GRCh38.p14",
            "gtf": f"work/download/annos/grch38/refseq/{wildcards.version}/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
        }


rule annos_features_refseq_gene_regions_process:  # -- process RefSeq gene regions files
    input:
        unpack(input_annos_features_refseq_gene_regions_process),
    output:
        tsv="work/annos/{genome_build}/features/refseq/{version}/refseq_genes.bed.gz",
        tsv_md5="work/annos/{genome_build}/features/refseq/{version}/refseq_genes.bed.gz.md5",
        tsv_tbi="work/annos/{genome_build}/features/refseq/{version}/refseq_genes.bed.gz.tbi",
        tsv_tbi_md5="work/annos/{genome_build}/features/refseq/{version}/refseq_genes.bed.gz.tbi.md5",
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
