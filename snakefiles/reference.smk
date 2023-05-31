# Download of reference FASTA files.

REFERENCE_URLS = {
    "grch37": (
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
        "phase2_reference_assembly_sequence/hs37d5.fa.gz"
    ),
    "GRCh38": (
        "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
        "seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
    ),
}

rule reference_download:
    output:
        download="reference/{genome_build}/reference/download/reference.fa.gz",
        fasta="reference/{genome_build}/reference/reference.fa",
        fasta_fai="reference/{genome_build}/reference/reference.fa.fai",
    run:
        ref_url = REFERENCE_URLS[wildcards.genome_build]
        shell(
            r"""
            aria2c \
                --out={output.download} \
                --split=8 \
                --max-concurrent-downloads=8 \
                --max-connection-per-server=8 \
                {ref_url}

            pigz -d -c {output.download} >{output.fasta}
            samtools faidx {output.fasta}
            """
        )
