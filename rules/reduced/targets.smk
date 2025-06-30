## Rules to create targets BED file for reduced output directory.


rule reduced_refseq_exons:  # -- create reduced data exons BED file
    input:
        unpack(input_annos_features_refseq_gene_regions_process),
    output:
        bed="output/reduced-{set_name}/targets/{genomebuild}/refseq/{version}/refseq_target_exons.bed",
        bed_md5="output/reduced-{set_name}/targets/{genomebuild}/refseq/{version}/refseq_target_exons.bed.md5",
        bed_gz="output/reduced-{set_name}/targets/{genomebuild}/refseq/{version}/refseq_target_exons.bed.gz",
        bed_gz_md5="output/reduced-{set_name}/targets/{genomebuild}/refseq/{version}/refseq_target_exons.bed.gz.md5",
        bed_gz_tbi="output/reduced-{set_name}/targets/{genomebuild}/refseq/{version}/refseq_target_exons.bed.gz.tbi",
        bed_gz_tbi_md5="output/reduced-{set_name}/targets/{genomebuild}/refseq/{version}/refseq_target_exons.bed.gz.tbi.md5",
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
        | egrep '^#|^X|^Y|^M|^[1-9]|^chrX|^chrY|^chrM|^chr[1-9]' \
        | egrep -v '_fix\s|_random\s|_alt\s' \
        | sort-bed - \
        | bedops --range {EXON_PADDING} --merge - \
        > {output.bed}

        bgzip -c {output.bed} \
        > {output.bed_gz}

        tabix -f {output.bed_gz}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_gz} >{output.bed_gz_md5}
        md5sum {output.bed_gz_tbi} >{output.bed_gz_tbi_md5}
        """
