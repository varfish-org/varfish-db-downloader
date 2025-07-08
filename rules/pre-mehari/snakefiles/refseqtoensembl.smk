rule result_noref_refseq_to_ensembl_tsv:
    input:
        tsv="work/genes/ensembl/{ensembl}/ensembl_xlink.tsv"
    output:
        tsv="output/pre-mehari/noref/refseqtoensembl/{ensembl}/RefseqToEnsembl.tsv",
        release_info="output/pre-mehari/noref/refseqtoensembl/{ensembl}/RefseqToEnsembl.release_info",
    shell:
        r"""
        awk -F $"\t" 'BEGIN{{OFS=FS}}{{print $3,$1,$2}}' {input.tsv} > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToEnsembl\t{wildcards.ensembl}\t\t-" > {output.release_info}
        """
