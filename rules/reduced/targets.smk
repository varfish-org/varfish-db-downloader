## Rules to create targets BED file for reduced output directory.

rule reduced_reduced_refseq_exons:  # -- create reduced data exons BED file
    input:
        unpack(input_annos_features_refseq_gene_regions_process),
    output:
        tsv="reduced-{set_name}/targets/{genome_build}/refseq/{version}/refseq_target_exoms.bed.gz",
        tsv_md5="reduced-{set_name}/targets/{genome_build}/refseq/{version}/refseq_target_exoms.bed.gz.md5",
        tsv_tbi="reduced-{set_name}/targets/{genome_build}/refseq/{version}/refseq_target_exoms.bed.gz.tbi",
        tsv_tbi_md5="reduced-{set_name}/targets/{genome_build}/refseq/{version}/refseq_target_exoms.bed.gz.tbi.md5",
    shell:
        r"""
        if [ "{wildcards.set_name}" == dev ]; then
            gene_symbol_re="^({DEV_GENE_SYMBOLS})$"
        else
            gene_symbol_re=".*"
        fi

        awk \
            -F $'\t' \
            -v "gene_symbol_re=$gene_symbol_re" \
            -f scripts/reduced-refseq-exons.awk \
            {input.acc} \
            <(zcat {input.gtf}) \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | sort-bed - \
        | bedops --range {EXON_PADDING} --merge - \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """
