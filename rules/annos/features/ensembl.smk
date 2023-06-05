## Rules related to ENSEMBL features.


rule annos_features_ensembl_gene_regions_grch37_download:  # -- download ENSEMBL gene regions files (GRCh37)
    output:
        gtf="work/download/annos/grch37/ensembl_genes/{version}/Homo_sapiens.GRCh37.{version}.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.gtf} \
            'https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.{wildcards.version}.gtf.gz'
        """


rule annos_features_ensembl_gene_regions_grch38_download:  # -- download ENSEMBL gene regions files (GRCh38)
    output:
        gtf="work/download/annos/grch38/ensembl_genes/{version}/Homo_sapiens.GRCh38.{ensembl_version}.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.gtf} \
            'https://ftp.ensembl.org/pub/release-{wildcards.ensembl_version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{wildcards.ensembl_version}.gtf.gz'
        """


def input_annos_features_ensembl_gene_regions_process(wildcards):
    """Input function for ``rule annos_features_ensembl_gene_regions_process``."""
    genome_release_nolower = {
        "grch37": "GRCh37",
        "grch38": "GRCh38",
    }[wildcards.genome_release]
    return {
        "gtf": (
            f"work/download/annos/{wildcards.genome_release}/ensembl_genes/{wildcards.version}/"
            f"Homo_sapiens.{genome_release_nolower}.{wildcards.version}.gtf.gz"
        )
    }


rule annos_features_ensembl_gene_regions_process:  # -- process ENSEMBL gene regions files
    input:
        unpack(input_annos_features_ensembl_gene_regions_process),
    output:
        tsv="work/annos/{genome_release}/features/ensembl/{version}/ensembl_genes.bed.gz",
        tsv_md5="work/annos/{genome_release}/features/ensembl/{version}/ensembl_genes.bed.gz.md5",
        tsv_tbi="work/annos/{genome_release}/features/ensembl/{version}/ensembl_genes.bed.gz.tbi",
        tsv_tbi_md5="work/annos/{genome_release}/features/ensembl/{version}/ensembl_genes.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/features-ensembl-gene-regions.awk \
            <(zcat {input.gtf}) \
        | (set +e; egrep '^#|^X|^Y|^M|^[1-9]|^chrX|^chrY|^chrM|^chr[1-9]'; set -e) \
        | sort-bed - \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
