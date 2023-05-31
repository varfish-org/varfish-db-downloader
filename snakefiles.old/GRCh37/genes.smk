rule GRCh37_refseq_genes_download:
    output:
        tsv="GRCh37/refseq_genes/r105/download/GeneInterval.tsv",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        wget --no-check-certificate -O - 'http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/Assembled_chromosomes/chr_accessions_GRCh37.p13' \
        | awk 'BEGIN {{ OFS="\t" }} !/^#/ {{ print $2, $1 }}' \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/names

        wget --no-check-certificate -O - 'http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz' \
        | zgrep -v '^#' \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/genes

        LC_ALL=C join --check-order -t $'\t' -j 1 $TMPDIR/names $TMPDIR/genes \
        | cut -f 2- \
        | sed -e 's/[=;:,]/ /g' \
        | awk 'BEGIN {{ OFS="\t" }} ($3 == "gene") {{ print $1, $4, $5, $15; }}' \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | sed -e 's/[";]//g' \
        | grep '^[1-9XYM]' \
        > {output.tsv}
        """


rule result_GRCh37_refseq_genes_tsv:
    input:
        tsv="GRCh37/refseq_genes/r105/download/GeneInterval.tsv",
        header="header/genes.txt",
    output:
        tsv="GRCh37/refseq_genes/r105/GeneInterval.tsv",
        release_info="GRCh37/refseq_genes/r105/GeneInterval.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="GRCh37\t"$1; $3=$3"\t\trefseq"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nRefSeqGenes\t105\tGRCh37\t" > {output.release_info}
        """


rule GRCh37_ensembl_genes_download:
    output:
        tsv="GRCh37/ensembl_genes/r104/download/GeneInterval.tsv",
    shell:
        r"""
        wget --no-check-certificate -O - 'http://ftp.ensembl.org/pub/grch37/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz' \
        | zcat \
        | awk '
            BEGIN {{ OFS="\t" }}
            ($0 ~ /protein_coding/ && $3 == "gene") {{
                print $1, $4, $5, $10;
            }}' \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | sed -e 's/[";]//g' \
        | grep '^[1-9XYM]' \
        > {output.tsv}
        """


rule result_GRCh37_ensembl_genes_tsv:
    input:
        tsv="GRCh37/ensembl_genes/r104/download/GeneInterval.tsv",
        header="header/genes.txt",
    output:
        tsv="GRCh37/ensembl_genes/r104/GeneInterval.tsv",
        release_info="GRCh37/ensembl_genes/r104/GeneInterval.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="GRCh37\t"$1; $3=$3"\t\tensembl"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblGenes\t104\tGRCh37\t" > {output.release_info}
        """
