## Rules related to gnomAD gene constraints.


rule genes_gnomad_download_v2_1_1:  # -- download gnomAD gene constraints v2.1.1
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


rule genes_gnomad_convert_v2_1_1:  # -- create gnomAD gene constraints TSV (v2.1.1)
    input:
        bgz="work/download/genes/gnomad/2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        xlink_ensembl=f"work/genes/ensembl/{DV.ensembl}/ensembl_xlink.tsv",
    output:
        tsv="work/genes/gnomad/2.1.1/gnomad_constraints.tsv",
        tsv_md5="work/genes/gnomad/2.1.1/gnomad_constraints.tsv.md5",
    run:
        run_genes_gnomad_constraints_v2_1_1_to_tsv(input, output, wildcards)


rule genes_gnomad_download_v4_0:  # -- download gnomAD gene constraints v4.0
    output:
        tsv="work/download/genes/gnomad/4.0/gnomad.v4.0.constraint_metrics.tsv",
        tsv_md5="work/download/genes/gnomad/4.0/gnomad.v4.0.constraint_metrics.tsv.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.tsv} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv

        md5sum {output.tsv} >{output.tsv_md5}
        """


def run_genes_gnomad_constraints_v4_0_to_tsv(input, output, wildcards):
    """Extra function because of snakefmt issues.

    Note that the names in the output file are taken from the v2.1.1 file.
    """
    columns_src_dst = {
        "transcript": "ensembl_transcript_id",
        "lof_hc_lc.exp": "exp_lof",
        "mis.exp": "exp_mis",
        "syn.exp": "exp_syn",
        "mis.z_score": "mis_z",
        "lof_hc_lc.obs": "obs_lof",
        "mis.obs": "obs_mis",
        "syn.obs": "obs_syn",
        "lof.oe": "oe_lof",
        "lof.oe_ci.lower": "oe_lof_lower",
        "lof.oe_ci.upper": "oe_lof_upper",
        "mis.oe": "oe_mis",
        "mis.oe_ci.lower": "oe_mis_lower",
        "mis.oe_ci.upper": "oe_mis_upper",
        "syn.oe": "oe_syn",
        "syn.oe_ci.lower": "oe_syn_lower",
        "syn.oe_ci.upper": "oe_syn_upper",
        "lof.pLI": "pLI",
        "syn.z_score": "syn_z",
    }
    columns_src_str = ",".join(columns_src_dst.keys())
    columns_tmp_str = ",".join(columns_src_dst.values())
    columns_dst = ["ensembl_gene_id", "entrez_id", "gene_symbol"] + list(columns_src_dst.values())[
        1:
    ]
    columns_dst_str = ",".join(columns_dst)
    shell(
        r"""
        export TMPDIR=/tmp/x
        mkdir -p $TMPDIR
        # trap "rm -rf $TMPDIR" EXIT

        # Extract transcripts that are MANE Select for genes that have one.
        head -n 1 {input.tsv} \
        > $TMPDIR/tmp.tsv
        tail -n +2 {input.tsv} \
        | sort -k3,3r \
        | sort -k1,1 -u \
        | sed -e 's/","/"__COMMA__"/g' \
        >> $TMPDIR/tmp.tsv

        # Convert to CSV so `qsv` has an easier time.
        cat $TMPDIR/tmp.tsv \
        | tr '\t' ',' \
        > $TMPDIR/tmp.txt

        qsv select {columns_src_str} $TMPDIR/tmp.txt \
        | qsv rename {columns_tmp_str} \
        | tr ',' '\t' \
        | sed -e 's/__COMMA__/,/g' \
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


rule genes_gnomad_convert_v4_0:  # -- create gnomAD gene constraints TSV (v4.0)
    input:
        tsv="work/download/genes/gnomad/4.0/gnomad.v4.0.constraint_metrics.tsv",
        xlink_ensembl=f"work/genes/ensembl/{DV.ensembl}/ensembl_xlink.tsv",
    output:
        tsv="work/genes/gnomad/4.0/gnomad_constraints.tsv",
        tsv_md5="work/genes/gnomad/4.0/gnomad_constraints.tsv.md5",
    run:
        run_genes_gnomad_constraints_v4_0_to_tsv(input, output, wildcards)
