## Rules related to ENSEMBL features.


def input_annos_features_ensembl_gene_regions_process(wildcards):
    """Input function for ``rule annos_features_ensembl_gene_regions_process``."""
    ensembl_version = {
        "grch37": DV.ensembl_37,
        "grch38": DV.ensembl_38,
    }[wildcards.genome_release]
    genomebuild = {
        "grch37": "GRCh37",
        "grch38": "GRCh38",
    }[wildcards.genome_release]
    return {
        "gtf": (
            f"work/genes/ensembl/{wildcards.genome_release}/{ensembl_version}/download/Homo_sapiens.{genomebuild}.{wildcards.version}.gtf.gz"
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
        | egrep -v '_gl|_alt|_random|Un|_fix' \
        | sort-bed - \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
