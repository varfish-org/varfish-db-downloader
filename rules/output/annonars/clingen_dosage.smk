## Output rules related to ClinGen dosage sensitivity regions.


rule annos_features_clingen_dosage_download_to_bed:  # -- convert ClinGen dosage sensitivity to BEd
    input:
        tsv="work/download/annos/{genome_release}/features/clingen_dosage/{date}/clingen_region_curation_list.tsv",
    output:
        bed="output/full/annonars/clingen-dosage-{genome_release}/{date}/clingen_region_curation_list.bed.gz",
        bed_md5="output/full/annonars/clingen-dosage-{genome_release}/{date}/clingen_region_curation_list.bed.gz.md5",
        bed_tbi="output/full/annonars/clingen-dosage-{genome_release}/{date}/clingen_region_curation_list.bed.gz.tbi",
        bed_tbi_md5="output/full/annonars/clingen-dosage-{genome_release}/{date}/clingen_region_curation_list.bed.gz.tbi.md5",
        spec_yaml="output/full/annonars/clingen-dosage-{genome_release}/{date}/clingen_region_curation_list.spec.yaml",
    shell:
        r"""
        if [[ "{wildcards.genome_release}" == "grch37" ]]; then
            chr_prefix=
        else
            chr_prefix=chr
        fi

        tail -n +8 {input.tsv} \
        | awk -v chr_prefix=$chr_prefix -F $'\t' 'BEGIN {{ OFS=FS }}
        {{
            if ($4 == "tbd") {{
                next;  /* skip, unmatched region */
            }}

            region=$4;
            split($4, a, /[:-]/);
            sub(/^chr/, "", a[1]);
            print chr_prefix a[1], a[2] - 1, a[3], $0;
        }}' \
        | LC_ALL=C sort -k1,1V -k2,2n \
        | bgzip -c \
        > {output.bed}
        tabix -f {output.bed}

        md5sum {output.bed} > {output.bed_md5}
        md5sum {output.bed_tbi} > {output.bed_tbi_md5}

        varfish-db-downloader tpl \
            --template rules/output/annonars/clingen_dosage.spec.yaml \
            --value today={wildcards.date} \
            --value genome_release={wildcards.genome_release} \
            \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
