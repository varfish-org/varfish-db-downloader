# Obtain current dump of HGNC gene information info.


rule grch3x_refseq_to_hgnc_download:
    output:
        jsongz="work/download/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{refseq_version}/{filename}",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.jsongz} \
            https://github.com/SACGF/cdot/releases/download/data_v{wildcards.cdot_release}/{wildcards.filename}
        """


rule result_grch3x_hgnc_to_tsv:
    input:
        header="rules/pre-mehari/header/hgnc.txt",
        tsv="work/download/annos/hgnc/{quarterly_release_date}/hgnc_complete_set.tsv",
    output:
        tsv="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{refseq_version}/Hgnc.tsv",
        release_info="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{refseq_version}/Hgnc.release_info",
    params:
        version=lambda wc: wc.quarterly_release_date,
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


def input_result_grch3x_refseq_to_hgnc_to_tsv(wildcards):
    hgnc_gffs = {
        "GRCh37": DV.cdot_refseq_gff_json_37,
        "GRCh38": DV.cdot_refseq_gff_json_38,
    }
    return {
        "header": "rules/pre-mehari/header/refseqtohgnc.txt",
        "json": (
            f"work/download/pre-mehari/{wildcards.genomebuild}/hgnc/"
            f"{wildcards.quarterly_release_date}+{wildcards.cdot_release}+{wildcards.refseq_version}"
            f"/cdot-{wildcards.cdot_release}.{hgnc_gffs[wildcards.genomebuild]}.json.gz"
        ),
    }


rule result_grch3x_refseq_to_hgnc_to_tsv:
    input:
        unpack(input_result_grch3x_refseq_to_hgnc_to_tsv)
    output:
        tsv="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{refseq_version}/RefseqToHgnc.tsv",
        release_info="output/pre-mehari/{genomebuild}/hgnc/{quarterly_release_date}+{cdot_release}+{refseq_version}/RefseqToHgnc.release_info",
    script:
        "scripts/cdot-to-tsv.py"