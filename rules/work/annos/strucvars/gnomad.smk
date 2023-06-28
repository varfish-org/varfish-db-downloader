## Rules related to gnomAD-SV.


rule annos_strucvars_gnomad_grch37_download:  # -- download gnomAD-SV files
    output:
        vcf="work/download/annos/grch37/strucvars/gnomad/2.1.1/gnomad_v2.1_sv.sites.vcf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz
        """


rule annos_strucvars_gnomad_grch37_process:  # -- process gnomAD-SV files
    input:
        vcf="work/download/annos/grch37/strucvars/gnomad/2.1.1/gnomad_v2.1_sv.sites.vcf.gz",
    output:
        bed="output/full/worker/track-strucvars-gnomad-grch37-2.1.1/gnomad.bed.gz",
        bed_md5="output/full/worker/track-strucvars-gnomad-grch37-2.1.1/gnomad.bed.gz.md5",
        bed_tbi="output/full/worker/track-strucvars-gnomad-grch37-2.1.1/gnomad.bed.gz.tbi",
        bed_tbi_md5="output/full/worker/track-strucvars-gnomad-grch37-2.1.1/gnomad.bed.gz.tbi.md5",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        echo -e "#chromosome\tbegin\tend\tsv_type\tn_homalt\tn_het" \
        > $TMPDIR/tmp.bed

        bcftools query \
            -e 'SVTYPE="MCNV"' \
            -f "%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\t%INFO/N_HOMALT\t%INFO/N_HET\n" \
            {input.vcf} \
        >> $TMPDIR/tmp.bed

        bgzip -c $TMPDIR/tmp.bed >{output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
