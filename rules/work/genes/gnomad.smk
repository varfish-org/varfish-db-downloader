## Rules related to gnomAD gene constraints.


rule genes_gnomad_download:  # -- download gnomAD gene constraints
    output:
        bgz="work/download/genes/gnomad/2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        bgz_md5="work/download/genes/gnomad/2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.bgz} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

        md5sum {output.bgz} >{output.bgz_md5}
        """


def run_genes_gnomad_constraints_v2_1_1_to_tsv(input, output, wildcards):
    """Extra function because of snakefmt issues."""
    columns_src = [
        "transcript",
        "exp_lof",
        "exp_mis",
        "exp_syn",
        "mis_z",
        "obs_lof",
        "obs_mis",
        "obs_syn",
        "oe_lof",
        "oe_lof_lower",
        "oe_lof_upper",
        "oe_mis",
        "oe_mis_lower",
        "oe_mis_upper",
        "oe_syn",
        "oe_syn_lower",
        "oe_syn_upper",
        "pLI",
        "syn_z",
        "exac_pLI",
        "exac_obs_lof",
        "exac_exp_lof",
        "exac_oe_lof",
    ]
    columns_src_str = ",".join(columns_src)
    columns_tmp = ["ensembl_transcript_id"] + columns_src[1:]
    columns_tmp_str = ",".join(columns_tmp)
    columns_dst = ["ensembl_gene_id", "entrez_id", "gene_symbol"] + columns_src[1:]
    columns_dst_str = ",".join(columns_dst)
    shell(
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        zcat {input.bgz} \
        | tr '\t' ',' \
        > $TMPDIR/tmp.txt

        qsv select {columns_src_str} $TMPDIR/tmp.txt \
        | qsv rename {columns_tmp_str} \
        | qsv sort -u \
        | tr ',' '\t' \
        > $TMPDIR/tmp.tsv

        qsv join -d '\t' \
            ensembl_transcript_id $TMPDIR/tmp.tsv \
            ensembl_transcript_id {input.xlink_ensembl} \
        | qsv select {columns_dst_str} \
        | tr ',' '\t' \
        > {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """
    )


rule genes_gnomad_convert:  # -- create gnomAD gene constraints TSV
    input:
        bgz="work/download/genes/gnomad/2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        xlink_ensembl=f"work/genes/ensembl/{DV.ensembl}/ensembl_xlink.tsv",
    output:
        tsv="work/genes/gnomad/2.1.1/gnomad_constraints.tsv",
        tsv_md5="work/genes/gnomad/2.1.1/gnomad_constraints.tsv.md5",
    run:
        run_genes_gnomad_constraints_v2_1_1_to_tsv(input, output, wildcards)
