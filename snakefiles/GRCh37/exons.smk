rule tmp_grch37_refseq_exons:
    output:
        bed="tmp/GRCh37/refseq_exons.bed",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        wget -O - 'http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/Assembled_chromosomes/chr_accessions_GRCh37.p13' \
        | awk 'BEGIN {{ OFS="\t" }} !/^#/ {{ print $2, $1 }}' \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/names

        wget -O - 'http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz' \
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

rule tmp_grch37_ensembl_exons:
    output:
        bed="tmp/GRCh37/ensembl_exons.bed",
    shell:
        r"""
        wget -O - 'http://ftp.ensembl.org/pub/grch37/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz' \
        | zcat \
        | awk '
            BEGIN {{ OFS="\t" }}
            ($3 == "exon") {{ print $1, $4, $5 }}' \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | sed -e 's/[";]//g' \
        | grep '^[1-9XYM]' \
        > {output.bed}
        """
