import polars as pl

def main(
    constraint_tsv_path, ensembl_xlink_tsv_path, output_tsv_path, output_tsv_md5_path
):
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

    additional_columns = ["exac_pLI", "exac_obs_lof", "exac_exp_lof", "exac_oe_lof"]

    columns_dst = (
        ["ensembl_gene_id", "entrez_id", "gene_symbol"]
        + list(columns_src_dst.values())[1:]
        + additional_columns
    )

    df = pl.read_csv(constraint_tsv_path, separator="\t", null_values="NA")
    ensembl_xlink = pl.read_csv(
        ensembl_xlink_tsv_path, separator="\t", null_values="NA"
    )

    df = (
        df.sort(by=["gene", "mane_select", "canonical", "cds_length", "transcript"])
        .filter(pl.col("gene_id").str.starts_with("ENSG"))
        .group_by(["gene"])
        .last()
        .sort(["gene", "mane_select", "canonical"])
        .with_columns(mane_or_canonical=pl.col("canonical") | pl.col("mane_select"))
    )
    discard = df.filter(pl.col("mane_or_canonical") == False).select(
        ["gene", "gene_id"]
    )
    print("No mane_select or canonical transcript found for ", discard.select("gene"))
    df = (
        df.filter(pl.col("mane_or_canonical") == True)
        .rename(columns_src_dst)
        .with_columns(pl.lit(None).alias(col) for col in additional_columns)
        .join(ensembl_xlink, on=["ensembl_transcript_id"])
        .select(columns_dst)
        .fill_null("NA")
    )
    df.write_csv(output_tsv_path, separator="\t")
    shell(f"md5sum {output_tsv_path} > {output_tsv_md5_path}")


if __name__ == "__main__":
    main(
        snakemake.input.tsv,
        snakemake.input.xlink_ensembl,
        snakemake.output.tsv,
        snakemake.output.tsv_md5,
    )
