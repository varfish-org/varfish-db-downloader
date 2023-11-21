## Rules to create build annonars functional annotation database..


rule work_annonars_functional_download_37:  # -- download functional data for GRCh37
    output:
        "work/download/refseq/grch37/{version}/{assembly}_genomic.gff.gz",
    shell:
        r"""
        wget -O {output} \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/{wildcards.version}/{wildcards.assembly}/{wildcards.assembly}_genomic.gff.gz
        """


rule work_annonars_functional_download_38:  # -- download functional data for GRCh37
    output:
        "work/download/refseq/grch38/{version}/{assembly}_genomic.gff.gz",
    shell:
        r"""
        wget -O {output} \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/{wildcards.version}/{wildcards.assembly}/{wildcards.assembly}_genomic.gff.gz
        """


def output_annonars_functional_input(wildcards):
    if wildcards.genome_release == "grch37":
        return f"work/download/refseq/grch37/{DV.refseq_fe_37}/GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
    else:
        return f"work/download/refseq/grch38/{DV.refseq_fe_38}/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"


rule output_annonars_functional:  # -- build annonars functional RocksDB file
    input:
        output_annonars_functional_input,
    output:
        rocksdb_identity=(
            "output/full/annonars/functional-{genome_release}-{v_refseq}+{v_annonars}/"
            "rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/full/annonars/functional-{genome_release}-{v_refseq}+{v_annonars}/spec.yaml"
        ),
    wildcard_constraints:
        v_refseq=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        zgrep '^#\|RefSeqFE' {input} > $TMPDIR/tmp.gff

        annonars functional import -vvv \
            --genome-release {wildcards.genome_release} \
            --path-in-gff $TMPDIR/tmp.gff \
            --path-out-rocksdb $(dirname {output.rocksdb_identity})

        varfish-db-downloader tpl \
            --template rules/output/annonars/functional.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.v_refseq}+{wildcards.v_annonars} \
            --value v_refseq={wildcards.v_refseq} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
