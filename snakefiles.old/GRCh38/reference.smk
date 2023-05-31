# Download GRCh38 reference sequence hs38.


rule grch38_download_hs38:
    output:
        fasta="GRCh38/reference/hs38/download/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz",
        md5="GRCh38/reference/hs38/download/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz.md5",
    log:
        "GRCh38/reference/hs38/download/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz.log",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.fasta} \
            -o {log} \
            http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

        pushd $(dirname {output.fasta})
        md5sum $(basename {output.fasta}) >$(basename {output.fasta}).md5
        """


rule grch38_get_hs38:
    input:
        fasta="GRCh38/reference/hs38/download/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz",
    output:
        fasta="GRCh38/reference/hs38/hs38.fa",
        md5="GRCh38/reference/hs38/hs38.fa.md5",
        fai="GRCh38/reference/hs38/hs38.fa.fai",
    shell:
        r"""
        GZIP="-q" zcat {input.fasta} > {output.fasta}

        samtools faidx {output.fasta}
        pushd $(dirname {output.fasta})
        md5sum $(basename {output.fasta}) >$(basename {output.fasta}).md5
        """
