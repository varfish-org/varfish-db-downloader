# Obtain current dump of HGNC gene information info.


rule grch3x_refseq_to_hgnc_download:
    output:
        jsongz="work/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{gff}/download/cdot-{cdot_release}.{gff}.json.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.jsongz} \
            https://github.com/SACGF/cdot/releases/download/data_v{wildcards.cdot_release}/cdot-{wildcards.cdot_release}.{wildcards.gff}.json.gz
        """


rule result_grch3x_hgnc_to_tsv:
    input:
        header="rules/pre-mehari/header/hgnc.txt",
        tsv="work/download/annos/hgnc/{quarterly_release_date}/hgnc_complete_set.tsv",
    output:
        tsv="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{gff}/Hgnc.tsv",
        release_info="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{gff}/Hgnc.release_info",
    params:
        version=lambda wildcards: wildcards.quarterly_release_date,
    shell:
        r"""
        # Check if the 22nd column header is "ucsc_id"
        header_check=$(head -n 1 {input.tsv} | awk -F $'\t' '{{print $22}}')
        if [ "$header_check" != "ucsc_id" ]; then
            echo "Error: Column 22 is not named 'ucsc_id'. Exiting."
            exit 1
        fi
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F $'\t' '
                BEGIN {{
                    OFS = FS
                }}
                {{
                    $NF = $NF "\t" substr($22, 1, length($22) - 2)
                    print
                }}'
        ) \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nHgnc\t{params.version}\t{wildcards.genomebuild}\t" > {output.release_info}
        """


rule result_grch3x_refseq_to_hgnc_to_tsv:
    input:
        header="rules/pre-mehari/header/refseqtohgnc.txt",
        json="work/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{gff}/download/cdot-{cdot_release}.{gff}.json.gz",
    output:
        tsv="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{gff}/RefseqToHgnc.tsv",
        release_info="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{gff}/RefseqToHgnc.release_info",
    params:
        version=lambda wildcards: wildcards.gff,
    script:
        "scripts/cdot-to-tsv.py"