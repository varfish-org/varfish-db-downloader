## Rules related to gnomAD-SV and gnomAD-CNV v4.


rule annos_strucvars_gnomad_sv_4_grch38_download:  # -- download gnomAD-SV >=4.1 files
    output:
        vcf="work/download/annos/grch38/strucvars/gnomad_sv/{v_gnomad}/gnomad.v{v_gnomad}.sv.{token}.vcf.gz",
    wildcard_constraints:
        v_gnomad=RE_VERSION,
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.v_gnomad}/genome_sv/gnomad.v{wildcards.v_gnomad}.sv.{wildcards.token}.vcf.gz
        """


rule annos_strucvars_gnomad_sv_4_grch38_process:  # -- process gnomAD-SV v4 files
    input:
        vcf="work/download/annos/grch38/strucvars/gnomad_sv/{version}/gnomad.v{version}.sv.sites.vcf.gz",
    output:
        bed=f"output/full/tracks/track-strucvars-gnomad-sv-grch38-{{version}}+{DV.tracks}/gnomad-sv.bed.gz",
        bed_md5=f"output/full/tracks/track-strucvars-gnomad-sv-grch38-{{version}}+{DV.tracks}/gnomad-sv.bed.gz.md5",
        bed_tbi=f"output/full/tracks/track-strucvars-gnomad-sv-grch38-{{version}}+{DV.tracks}/gnomad-sv.bed.gz.tbi",
        bed_tbi_md5=f"output/full/tracks/track-strucvars-gnomad-sv-grch38-{{version}}+{DV.tracks}/gnomad-sv.bed.gz.tbi.md5",
    wildcard_constraints:
        version=RE_VERSION,
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        echo -e "#chromosome\tbegin\tend\tsv_type\tmale_n_homref\tmale_n_het\tmale_n_homalt\tmale_n_hemiref\tmale_n_hemialt\tfemale_n_homref\tfemale_n_het\tfemale_n_homalt\tcnv_n_total\tcnv_n_var" \
        > $TMPDIR/tmp.bed

        bcftools query \
            -f "%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\t%INFO/N_HOMREF_XY\t%INFO/N_HET_XY\t%INFO/N_HOMALT_XY\t%INFO/N_HEMIREF_XY\t%INFO/N_HEMIALT_XY\t%INFO/N_HOMREF_XY\t%INFO/N_HET_XY\t%INFO/N_HOMALT_XY\t%CN_NUMBER\t%CN_COUNT\n" \
            {input.vcf} \
        | awk -v OFS='\t' '{{
            if ($5 == ".") {{ $5 = 0; }}
            if ($6 == ".") {{ $6 = 0; }}
            if ($7 == ".") {{ $7 = 0; }}
            if ($8 == ".") {{ $8 = 0; }}
            if ($9 == ".") {{ $9 = 0; }}
            if ($10 == ".") {{ $10 = 0; }}
            if ($11 == ".") {{ $11 = 0; }}
            if ($12 == ".") {{ $12 = 0; }}
            if ($13 == ".") {{ $13 = 0; }}
            if ($14 == ".") {{
                $14 = 0
            }} else {{
                sum = 0
                split($14, a, ",")
                for (x in a) {{
                    sum += x
                }}
                $14 = sum
            }}
            print $0
        }}' \
        | sed -e 's/CPX/BND/g' -e 's/CTX/BND/g' \
        >> $TMPDIR/tmp.bed

        bgzip -c $TMPDIR/tmp.bed >{output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule annos_strucvars_gnomad_cnv_4_grch38_download:  # -- download gnomAD-CNV >=4.1 files
    output:
        vcf="work/download/annos/grch38/strucvars/gnomad_cnv/{v_gnomad}/gnomad.v{v_gnomad}.cnv.{token}.vcf.gz",
    wildcard_constraints:
        v_gnomad=RE_VERSION,
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.v_gnomad}/exome_cnv/gnomad.v{wildcards.v_gnomad}.cnv.{wildcards.token}.vcf.gz
        """


rule annos_strucvars_gnomad_cnv_4_grch38_process:  # -- process gnomAD-CNV >=4.1 files
    input:
        vcf="work/download/annos/grch38/strucvars/gnomad_cnv/{version}/gnomad.v{version}.cnv.all.vcf.gz",
    output:
        bed=f"output/full/tracks/track-strucvars-gnomad-cnv-grch38-{{version}}+{DV.tracks}/gnomad-cnv.bed.gz",
        bed_md5=f"output/full/tracks/track-strucvars-gnomad-cnv-grch38-{{version}}+{DV.tracks}/gnomad-cnv.bed.gz.md5",
        bed_tbi=f"output/full/tracks/track-strucvars-gnomad-cnv-grch38-{{version}}+{DV.tracks}/gnomad-cnv.bed.gz.tbi",
        bed_tbi_md5=f"output/full/tracks/track-strucvars-gnomad-cnv-grch38-{{version}}+{DV.tracks}/gnomad-cnv.bed.gz.tbi.md5",
    wildcard_constraints:
        version=RE_VERSION,
    shell:
        r"""
        set -euo pipefail

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        echo -e "#chromosome\tbegin\tend\tsv_type\tn_total\tn_var" \
        > $TMPDIR/tmp.bed

        bcftools query \
            -f "%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\t%INFO/SN\t%INFO/SC\n" \
            {input.vcf} \
        >> $TMPDIR/tmp.bed

        bgzip -c $TMPDIR/tmp.bed >{output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
