# Obtain current dump of mim2gene_medgen file from NCBI.


rule result_grchXX_mim2gene_medgen_to_tsv:
    input:
        header="rules/pre-mehari/header/mim2gene.txt",
        txt="work/download/genes/ncbi/{date}/mim2gene_medgen",
    output:
        tsv="output/pre-mehari/noref/mim2gene/{date}/Mim2geneMedgen.tsv",
        release_info="output/pre-mehari/noref/mim2gene/{date}/Mim2geneMedgen.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.txt}
        ) \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nMim2geneMedgen\t{wildcards.date}\t\t-" > {output.release_info}
        """
