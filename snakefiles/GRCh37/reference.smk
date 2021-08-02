# Download GRCh37 reference sequence hs37d5 and GRCh38 reference sequence hs38.


rule grch37_download_hs37d5:
    output:
        fasta="GRCh37/reference/hs37d5/hs37d5.fa.gz",
        md5="GRCh37/reference/hs37d5/hs37d5.fa.gz.md5",
    log:
        "GRCh37/reference/hs37d5/hs37d5.fa.gz.log",
    shell:
        r"""
        wget \
            -O {output.fasta} \
            -o {log} \
            http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

        pushd $(dirname {output.fasta})
        md5sum $(basename {output.fasta}) >$(basename {output.fasta}).md5
        """


rule grch37_get_hs37d5:
    input:
        "GRCh37/reference/hs37d5/hs37d5.fa.gz",
    output:
        fasta="GRCh37/reference/hs37d5/hs37d5.fa",
        fai="GRCh37/reference/hs37d5/hs37d5.fa.fai",
    log:
        fasta="GRCh37/reference/hs37d5/hs37d5.fa.log",
    shell:
        r"""
        GZIP="-q" zcat {input} \
        > {output.fasta} || true

        samtools faidx {output.fasta}
        pushd $(dirname {output.fasta})
        md5sum $(basename {output.fasta}) >$(basename {output.fasta}).md5
        """
