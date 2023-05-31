rule tracks_grch37_ucsc_genomic_super_dups:
    output:
        bed="tracks/grch37/ucsc_genomicSuperDups.bed.gz",
        bed_md5="tracks/grch37/ucsc_genomicSuperDups.bed.gz.md5",
        bed_tbi="tracks/grch37/ucsc_genomicSuperDups.bed.gz.tbi",
        bed_tbi_md5="tracks/grch37/ucsc_genomicSuperDups.bed.gz.tbi.md5",
        txt=temp("tracks/grch37/download/genomicSuperDups.txt.gz"),
    shell:
        r"""
        set -x

        mkdir -p $(dirname {output.txt})
        wget -O {output.txt} https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/genomicSuperDups.txt.gz

        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {output.txt} \
            | cut -f 2,3,4,5 \
            | sed -e 's/^chr//g' \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule tracks_grch37_ucsc_rmsk:
    output:
        bed="tracks/grch37/ucsc_rmsk.bed.gz",
        bed_md5="tracks/grch37/ucsc_rmsk.bed.gz.md5",
        bed_tbi="tracks/grch37/ucsc_rmsk.bed.gz.tbi",
        bed_tbi_md5="tracks/grch37/ucsc_rmsk.bed.gz.tbi.md5",
        txt=temp("tracks/grch37/download/rmsk.txt.gz"),
    shell:
        r"""
        set -x

        mkdir -p $(dirname {output.txt})
        wget -O {output.txt} https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/rmsk.txt.gz

        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {output.txt} \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ if ($12 == $13) {{ label = $13 "/" $11 }} else {{ label = $12 "/" $13 "/" $11 }} print $6, $7, $8, label }}' \
            | sed -e 's/^chr//g' \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule tracks_grch37_ucsc_x_seq_lift_over_psl:
    output:
        bed="tracks/grch37/ucsc_{prefix}SeqLiftOverPsl.bed.gz",
        bed_md5="tracks/grch37/ucsc_{prefix}SeqLiftOverPsl.bed.gz.md5",
        bed_tbi="tracks/grch37/ucsc_{prefix}SeqLiftOverPsl.bed.gz.tbi",
        bed_tbi_md5="tracks/grch37/ucsc_{prefix}SeqLiftOverPsl.bed.gz.tbi.md5",
        txt=temp("tracks/grch37/download/{prefix}SeqLiftOverPsl.txt.gz"),
    shell:
        r"""
        set -x

        mkdir -p $(dirname {output.txt})
        wget -O {output.txt} https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/{wildcards.prefix}SeqLiftOverPsl.txt.gz

        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {output.txt} \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ if ($11 ~ /{wildcards.prefix}/ && $15 !~ /random/ && $15 !~ /hap/) {{ print $15, $17, $18, $11 }} }}' \
            | awk -F $'\t' -f scripts/sort-bed.awk \
            | sed -e 's/^chr//g' \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
