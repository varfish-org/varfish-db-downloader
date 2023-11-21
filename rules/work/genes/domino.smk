## Rules related to DOMINO gene annotation.


rule genes_domino:  # -- postprocess file for import
    input:
        tsv="bundled-data/domino/score_all_final_19.02.19.txt",
    output:
        tsv="work/genes/domino/20190219/domino.tsv",
        tsv_md5="work/genes/domino/20190219/domino.tsv.md5",
    shell:
        """
        cut -f 1-2 {input.tsv} \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv}.md5
        """
