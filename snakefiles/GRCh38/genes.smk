rule GRCh38_refseq_genes_download:
    output:
        tsv="GRCh38/refseq_genes/r39/download/GeneInterval.tsv",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        wget -O - 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc' \
        | awk 'BEGIN {{ OFS="\t" }} !/^#/ {{ print $2, $1 }}' \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/names

        wget -O - 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz' \
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


rule result_GRCh38_refseq_genes_tsv:
    input:
        tsv="GRCh38/refseq_genes/r39/download/GeneInterval.tsv",
        header="header/genes.txt",
    output:
        tsv="GRCh38/refseq_genes/r39/GeneInterval.tsv",
        release_info="GRCh38/refseq_genes/r39/GeneInterval.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="GRCh38\t"$1; $3=$3"\t\trefseq"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nRefSeqGenes\t39\tGRCh38\t" > {output.release_info}
        """


rule GRCh38_ensembl_genes_download:
    output:
        tsv="GRCh38/ensembl_genes/r104/download/GeneInterval.tsv",
    shell:
        r"""
        wget -O - 'ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz' \
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


rule result_GRCh38_ensembl_genes_tsv:
    input:
        tsv="GRCh38/ensembl_genes/r104/download/GeneInterval.tsv",
        header="header/genes.txt",
    output:
        tsv="GRCh38/ensembl_genes/r104/GeneInterval.tsv",
        release_info="GRCh38/ensembl_genes/r104/GeneInterval.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="GRCh38\t"$1; $3=$3"\t\tensembl"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblGenes\t104\tGRCh38\t" > {output.release_info}
        """
