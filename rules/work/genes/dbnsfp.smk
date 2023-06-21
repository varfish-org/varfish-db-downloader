## Rules related to dbNSFP gene information.


rule genes_dbnsfp_genes_copy:  # -- fixup dbNSFP gene file inconsistency (`s/NA/./g`)
    input:
        tsv="work/download/annos/grch37/seqvars/dbnsfp/{version}a/dbNSFP{version}_gene.complete.gz",
    output:
        tsv="work/genes/dbnsfp/{version}/genes.tsv.gz",
        tsv_md5="work/genes/dbnsfp/{version}/genes.tsv.gz.md5",
    shell:
        r"""
        zcat {input.tsv} \
        | tr -d '\r' \
        | sed \
            -e ':repeat; s/\t\NA\t/\t.\t/g; t repeat' \
            -e 's/\tNA$/\t./g' \
        | pigz -c \
        > {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """
