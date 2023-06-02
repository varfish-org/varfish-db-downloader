## Rules related to features from UCSC Genome Browser.


rule features_ucsc_grch37_download:  # -- download of UCSC hg19 tracks
    output:
        txt="work/download/annos/grch37/features/ucsc/{version}/{filename}",
    shell:
        r"""
        # Check version.
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        wget -O $TMPDIR/listing.html https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database
        version=$(grep {wildcards.filename} $TMPDIR/listing.html | awk '{{ print $$2 }}')
        if [[ "$version" != "{wildcards.version}" ]]; then
            >&2 echo "Version mismatch for {wildcards.filename}: expected {version}, got $version"
            exit 1
        fi

        # Actually perform the download.
        wget -O {output.txt} https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/{wildcards.filename}
        """


rule features_ucsc_genomic_super_dups_grch37_process:  # -- processing of UCSC hg19 genomicSuperDups
    input:
        txt="work/download/annos/grch37/features/ucsc/{version}/genomicSuperDups.txt.gz",
    output:
        bed="work/annos/grch37/features/ucsc/{version}/genomicSuperDups.bed.gz",
        bed_md5="work/annos/grch37/features/ucsc/{version}/genomicSuperDups.bed.gz.md5",
        bed_tbi="work/annos/grch37/features/ucsc/{version}/genomicSuperDups.bed.gz.tbi",
        bed_tbi_md5="work/annos/grch37/features/ucsc/{version}/genomicSuperDups.bed.gz.tbi.md5",
    shell:
        r"""
        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {input.txt} \
            | cut -f 2,3,4,5 \
            | sed -e 's/^chr//g' \
            | (set +e; egrep '^#|^X|^Y|^M|^[1-9]'; set -e) \
            | (set +e; egrep -v '^Un|_random|_fix|_alt|_hap'; set -e) \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule features_ucsc_rmsk_grch37_process:  # -- processing of UCSC hg19 rmsk
    input:
        txt="work/download/annos/grch37/features/ucsc/{version}/rmsk.txt.gz",
    output:
        bed="work/annos/grch37/features/ucsc/{version}/rmsk.bed.gz",
        bed_md5="work/annos/grch37/features/ucsc/{version}/rmsk.bed.gz.md5",
        bed_tbi="work/annos/grch37/features/ucsc/{version}/rmsk.bed.gz.tbi",
        bed_tbi_md5="work/annos/grch37/features/ucsc/{version}/rmsk.bed.gz.tbi.md5",
    shell:
        r"""
        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {input.txt} \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ if ($12 == $13) {{ label = $13 "/" $11 }} else {{ label = $12 "/" $13 "/" $11 }} print $6, $7, $8, label }}' \
            | sed -e 's/^chr//g' \
            | (set +e; egrep '^#|^X|^Y|^M|^[1-9]'; set -e) \
            | (set +e; egrep -v '^Un|_random|_fix|_alt|_hap'; set -e) \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule features_ucsc_liftover_grch37_process:  # -- process of UCSC hg19 *SeqLiftOverPsl
    input:
        txt="work/download/annos/grch37/features/ucsc/{version}/{prefix}SeqLiftOverPsl.bed.gz",
    output:
        bed="work/annos/grch37/features/ucsc/{version}/{prefix}SeqLiftOverPsl.bed.gz",
        bed_md5="work/annos/grch37/features/ucsc/{version}/{prefix}SeqLiftOverPsl.bed.gz.md5",
        bed_tbi="work/annos/grch37/features/ucsc/{version}/{prefix}SeqLiftOverPsl.bed.gz.tbi",
        bed_tbi_md5="work/annos/grch37/features/ucsc/{version}/{prefix}SeqLiftOverPsl.bed.gz.tbi.md5",
    shell:
        r"""
        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {input.txt} \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ if ($11 ~ /{wildcards.prefix}/ && $15 !~ /random/ && $15 !~ /hap/) {{ print $15, $17, $18, $11 }} }}' \
            | sort-bed - \
            | sed -e 's/^chr//g' \
            | (set +e; egrep '^#|^X|^Y|^M|^[1-9]'; set -e) \
            | (set +e; egrep -v '^Un|_random|_fix|_alt|_hap'; set -e) \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
