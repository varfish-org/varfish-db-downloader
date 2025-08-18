## Rules related to RefSeq features.


rule annos_features_refseq_gene_regions_download_grch37:
    output:
        acc="work/download/annos/grch37/refseq/{version}/chr_accessions_{assembly}",
        gtf="work/download/annos/grch37/refseq/{version}/{assembly}_genomic.gtf.gz",
    params:
        version_stripped=lambda wildcards: wildcards.version.partition(".")[0],
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.acc} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.{params.version_stripped}/Assembled_chromosomes/chr_accessions_GRCh37.p13'

        wget --no-check-certificate \
            -O {output.gtf} \
            '{DV.refseq_base_url}/{wildcards.version}/{wildcards.assembly}/{wildcards.assembly}_genomic.gtf.gz'
        """


rule annos_features_refseq_gene_regions_download_grch38:
    output:
        report="work/download/annos/grch38/refseq/{version}/{assembly}_assembly_report.txt",
        acc="work/download/annos/grch38/refseq/{version}/chr_accessions_{assembly}",
        gtf="work/download/annos/grch38/refseq/{version}/{assembly}_genomic.gtf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.report} \
            "{DV.refseq_base_url}/{wildcards.version}/{wildcards.assembly}/{wildcards.assembly}_assembly_report.txt"

        echo -e "#Chromosome\tRefSeq Accession.version\tRefSeq\tgi\tGenBank Accession.version\tGenBank gi" \
        > {output.acc}
        cat {output.report} \
        | tr -d '\r' \
        | awk -F $'\t' 'BEGIN {{ OFS=FS }} ($1 !~ /^#/) {{ print $10, $7, $9, $5, "." }}' \
        >> {output.acc}

        wget --no-check-certificate \
            -O {output.gtf} \
            '{DV.refseq_base_url}/{wildcards.version}/{wildcards.assembly}/{wildcards.assembly}_genomic.gtf.gz'
        """


def input_annos_features_refseq_gene_regions_process(wildcards):
    """Input function for ``rule annos_features_refseq_gene_regions_process``."""
    if wildcards.genomebuild == "grch37":
        assembly = DV.refseq_ref_37_assembly
    else:
        assembly = DV.refseq_ref_38_assembly

    return {
        "acc": f"work/download/annos/{wildcards.genomebuild}/refseq/{wildcards.version}/chr_accessions_{assembly}",
        "gtf": f"work/download/annos/{wildcards.genomebuild}/refseq/{wildcards.version}/{assembly}_genomic.gtf.gz",
    }


rule annos_features_refseq_gene_regions_process:  # -- process RefSeq gene regions files
    input:
        unpack(input_annos_features_refseq_gene_regions_process),
    output:
        tsv="work/annos/{genomebuild}/features/refseq/{version}/refseq_genes.bed.gz",
        tsv_md5="work/annos/{genomebuild}/features/refseq/{version}/refseq_genes.bed.gz.md5",
        tsv_tbi="work/annos/{genomebuild}/features/refseq/{version}/refseq_genes.bed.gz.tbi",
        tsv_tbi_md5="work/annos/{genomebuild}/features/refseq/{version}/refseq_genes.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/features-refseq-gene-regions.awk \
            {input.acc} \
            <(zcat {input.gtf}) \
        | (set +e; egrep '^#|^X|^Y|^M|^[1-9]|^chrX|^chrY|^chrM|^chr[1-9]'; set -e) \
        | egrep -v '_gl|_alt|_random|Un|fix' \
        | sort-bed - \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
