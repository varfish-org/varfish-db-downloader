rule tmp_grch38_refseq_exons:
    output:
        bed="tmp/GRCh38/refseq_exons.bed",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        wget --no-check-certificate -O - 'http://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc' \
        | awk 'BEGIN {{ OFS="\t" }} !/^#/ {{ print $2, $1 }}' \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/names

        wget --no-check-certificate -O - 'http://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz' \
        | zgrep -v '^#' \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/genes

        LC_ALL=C join --check-order -t $'\t' -j 1 $TMPDIR/names $TMPDIR/genes \
        | cut -f 2- \
        | sed -e 's/[=;:,]/ /g' \
        | awk 'BEGIN {{ OFS="\t" }} ($3 == "exon") {{ print $1, $4, $5 }}' \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | sed -e 's/[";]//g' \
        | grep '^[1-9XYM]' \
        > {output.bed}
        """


rule tmp_grch38_ensembl_exons:
    output:
        bed="tmp/GRCh38/ensembl_exons.bed",
    shell:
        r"""
        wget --no-check-certificate -O - 'http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz' \
        | zcat \
        | awk '
            BEGIN {{ OFS="\t" }}
            ($3 == "exon") {{ print $1, $4, $5 }}' \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | sed -e 's/[";]//g' \
        | grep '^[1-9XYM]' \
        > {output.bed}
        """
