## Rules related to human reference genome sequence.

#: Download URLs
REFERENCE_URLS = {
    "grch37": (
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
        "phase2_reference_assembly_sequence/hs37d5.fa.gz"
    ),
    "grch38": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
        "seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
    ),
}


rule reference_download:  # -- download reference genome sequence
    output:
        download="work/download/reference/{genome_build}/reference.fa.gz",
    run:
        ref_url = REFERENCE_URLS[wildcards.genome_build]
        shell(
            r"""
            aria2c \
                --check-certificate=false \
                --file-allocation=trunc \
                --out={output.download} \
                --split=8 \
                --max-concurrent-downloads=8 \
                --max-connection-per-server=8 \
                {ref_url}
            """
        )


rule reference_process:  # -- post-process reference sequence after download
    input:
        download="work/download/reference/{genome_build}/reference.fa.gz",
    output:
        fasta="work/reference/{genome_build}/reference.fa",
        fasta_md5="work/reference/{genome_build}/reference.fa.md5",
        fasta_fai="work/reference/{genome_build}/reference.fa.fai",
        fasta_fai_md5="work/reference/{genome_build}/reference.fa.fai.md5",
    shell:
        r"""
        pigz -d -c {input.download} >{output.fasta}

        samtools faidx {output.fasta}

        md5sum {output.fasta} >{output.fasta_md5}
        md5sum {output.fasta_fai} >{output.fasta_fai_md5}
        """
