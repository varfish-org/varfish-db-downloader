## Rules related to features from UCSC Genome Browser.


rule features_ucsc_download:  # -- download of UCSC tracks
    output:
        txt="work/download/annos/{genome_release}/features/ucsc/{version}/{filename}",
    shell:
        r"""
        # Check version.
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        if [[ "{wildcards.genome_release}" == "grch37" ]]; then
            ucsc_name=hg19
        else
            ucsc_name=hg38
        fi

        wget -O $TMPDIR/listing.html https://hgdownload.cse.ucsc.edu/goldenpath/$ucsc_name/database
        version=$(grep {wildcards.filename} $TMPDIR/listing.html| awk '{{ gsub(/-/, "", $3); print $3 }}')
        if [[ "$version" != "{wildcards.version}" ]]; then
            >&2 echo "Version mismatch for {wildcards.filename}: expected {version}, got $version"
            exit 1
        fi

        # Actually perform the download.
        wget -O {output.txt} https://hgdownload.cse.ucsc.edu/goldenpath/$ucsc_name/database/{wildcards.filename}
        """


rule features_ucsc_genomic_super_dups_process:  # -- processing of UCSC genomicSuperDups
    input:
        txt="work/download/annos/{genome_release}/features/ucsc/{version}/genomicSuperDups.txt.gz",
    output:
        bed="output/full/worker/track-features-ucsc-genomicsuperdups-{genome_release}-{version}/genomicSuperDups.bed.gz",
        bed_md5="output/full/worker/track-features-ucsc-genomicsuperdups-{genome_release}-{version}/genomicSuperDups.bed.gz.md5",
        bed_tbi="output/full/worker/track-features-ucsc-genomicsuperdups-{genome_release}-{version}/genomicSuperDups.bed.gz.tbi",
        bed_tbi_md5="output/full/worker/track-features-ucsc-genomicsuperdups-{genome_release}-{version}/genomicSuperDups.bed.gz.tbi.md5",
    shell:
        r"""
        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {input.txt} \
            | cut -f 2,3,4,5 \
            | {{ if [[ "{wildcards.genome_release}" == "grch37" ]]; then sed -e 's/^chr//g'; else cat; fi }} \
            | (set +e; egrep '^#|^X|^Y|^M|^[1-9]|^chrX|^chrY|^chrM|^chr[1-9]'; set -e) \
            | (set +e; egrep -v '^Un|^chrU|_random|_fix|_alt|_hap'; set -e) \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule features_ucsc_rmsk_process:  # -- processing of UCSC rmsk
    input:
        txt="work/download/annos/{genome_release}/features/ucsc/{version}/rmsk.txt.gz",
    output:
        bed="output/full/worker/track-features-ucsc-rmsk-{genome_release}-{version}/rmsk.bed.gz",
        bed_md5="output/full/worker/track-features-ucsc-rmsk-{genome_release}-{version}/rmsk.bed.gz.md5",
        bed_tbi="output/full/worker/track-features-ucsc-rmsk-{genome_release}-{version}/rmsk.bed.gz.tbi",
        bed_tbi_md5="output/full/worker/track-features-ucsc-rmsk-{genome_release}-{version}/rmsk.bed.gz.tbi.md5",
    shell:
        r"""
        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {input.txt} \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ if ($12 == $13) {{ label = $13 "/" $11 }} else {{ label = $12 "/" $13 "/" $11 }} print $6, $7, $8, label }}' \
            | {{ if [[ "{wildcards.genome_release}" == "grch37" ]]; then sed -e 's/^chr//g'; else cat; fi }} \
            | (set +e; egrep '^#|^X|^Y|^M|^[1-9]|^chrX|^chrY|^chrM|^chr[1-9]'; set -e) \
            | (set +e; egrep -v '^Un|^chrU|_random|_fix|_alt|_hap'; set -e) \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule features_ucsc_liftover_process:  # -- process of UCSC *SeqLiftOverPsl
    input:
        txt="work/download/annos/{genome_release}/features/ucsc/{version}/{prefix}SeqLiftOverPsl.txt.gz",
    output:
        bed="output/full/worker/track-features-ucsc-{prefix}seqliftoverpsl-{genome_release}-{version}/{prefix}SeqLiftOverPsl.bed.gz",
        bed_md5="output/full/worker/track-features-ucsc-{prefix}seqliftoverpsl-{genome_release}-{version}/{prefix}SeqLiftOverPsl.bed.gz.md5",
        bed_tbi="output/full/worker/track-features-ucsc-{prefix}seqliftoverpsl-{genome_release}-{version}/{prefix}SeqLiftOverPsl.bed.gz.tbi",
        bed_tbi_md5="output/full/worker/track-features-ucsc-{prefix}seqliftoverpsl-{genome_release}-{version}/{prefix}SeqLiftOverPsl.bed.gz.tbi.md5",
    shell:
        r"""
        (
            echo -e "#chrom\tbegin\tend\tlabel"
            zcat {input.txt} \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ if ($11 ~ /{wildcards.prefix}/ && $15 !~ /random/ && $15 !~ /hap/) {{ print $15, $17, $18, $11 }} }}' \
            | sort-bed - \
            | {{ if [[ "{wildcards.genome_release}" == "grch37" ]]; then sed -e 's/^chr//g'; else cat; fi }} \
            | (set +e; egrep '^#|^X|^Y|^M|^[1-9]|^chrX|^chrY|^chrM|^chr[1-9]'; set -e) \
            | (set +e; egrep -v '^Un|^chrU|_random|_fix|_alt|_hap'; set -e) \
        ) \
        | bgzip -c \
        > {output.bed}

        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
