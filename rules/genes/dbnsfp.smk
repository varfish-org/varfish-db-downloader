## Rules related to dbNSFP gene information.


rule genes_dbnsfp_genes_copy:  # -- copy over dbNSFP genes file
    input:
        tsv="work/download/grch37/seqvars/dbnsfp/{version}a/dbNSFP{version}_gene.complete.gz",
    output:
        tsv="work/genes/dbnsfp/{version}/genes.tsv.gz",
        tsv_md5="work/genes/dbnsfp/{version}/genes.tsv.gz.md5",
    shell:
        r"""
        zcat {input.tsv} \
        | pigz -c \
        > {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """
