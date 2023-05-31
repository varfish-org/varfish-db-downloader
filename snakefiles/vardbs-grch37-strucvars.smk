rule vardbs_grch37_struc_vars_clinvar:
    output:
        bed="vardbs/grch37/strucvar/clinvar.bed.gz",
        bed_md5="vardbs/grch37/strucvar/clinvar.bed.gz.md5",
        bed_tbi="vardbs/grch37/strucvar/clinvar.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvar/clinvar.bed.gz.tbi.md5",
    shell:
        r"""
        set -x
        export LC_ALL=C

        awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-clinvar.awk \
            vardbs/grch37/strucvar/clinvar.tsv \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule vardbs_grch37_struc_vars_dbvar:
    output:
        download=expand(
            "vardbs/grch37/strucvar/download/GRCh37.nr_{type}.tsv.gz",
            type=["deletions", "duplications", "insertions"],
        ),
        bed="vardbs/grch37/strucvar/dbvar.bed.gz",
        bed_md5="vardbs/grch37/strucvar/dbvar.bed.gz.md5",
        bed_tbi="vardbs/grch37/strucvar/dbvar.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvar/dbvar.bed.gz.tbi.md5",
    shell:
        r"""
        set -x
        export LC_ALL=C

        for dst in {output.download}; do
            type=$(basename $dst | cut -d _ -f 2 | cut -d . -f 1)
            wget --no-check-certificate \
                -O $dst \
                https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/$type/GRCh37.nr_$type.tsv.gz
        done

        awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-dbvar.awk \
            <(zcat {output.download}) \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule vardbs_grch37_struc_vars_dgv:
    output:
        download_txt="vardbs/grch37/strucvar/download/GRCh37_hg19_variants_2020-02-25.txt",
        bed="vardbs/grch37/strucvar/dgv.bed.gz",
        bed_md5="vardbs/grch37/strucvar/dgv.bed.gz.md5",
        bed_tbi="vardbs/grch37/strucvar/dgv.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvar/dgv.bed.gz.tbi.md5",
    shell:
        r"""
        set -x
        export LC_ALL=C

        wget --no-check-certificate \
            -O {output.download_txt} \
            http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt

        awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-dgv.awk \
            {output.download_txt} \
        | grep -v _gl \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule vardbs_grch37_struc_vars_dgv_gs:
    output:
        download_gff3="vardbs/grch37/strucvar/download/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
        bed="vardbs/grch37/strucvar/dgv_gs.bed.gz",
        bed_md5="vardbs/grch37/strucvar/dgv_gs.bed.gz.md5",
        bed_tbi="vardbs/grch37/strucvar/dgv_gs.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvar/dgv_gs.bed.gz.tbi.md5",
    shell:
        r"""
        set -x
        export LC_ALL=C

        wget --no-check-certificate \
            -O {output.download_gff3} \
            http://dgv.tcag.ca/dgv/docs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3

        awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-dgv_gs.awk \
            {output.download_gff3} \
        | grep -v _gl \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule vardbs_grch37_struc_vars_exac:
    output:
        download_bed="vardbs/grch37/strucvar/download/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
        bed="vardbs/grch37/strucvar/exac.bed.gz",
        bed_md5="vardbs/grch37/strucvar/exac.bed.gz.md5",
        bed_tbi="vardbs/grch37/strucvar/exac.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvar/exac.bed.gz.tbi.md5",
    shell:
        r"""
        set -x
        export LC_ALL=C

        wget --no-check-certificate \
            -O {output.download_bed} \
            ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/cnv/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed

        awk \
            -f scripts/vardbs-grch37-strucvar-exac.awk \
            {output.download_bed} \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule vardbs_grch37_struc_vars_g1k:
    output:
        vcf="vardbs/grch37/strucvar/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
        bed="vardbs/grch37/strucvar/g1k.bed.gz",
        bed_md5="vardbs/grch37/strucvar/g1k.bed.gz.md5",
        bed_tbi="vardbs/grch37/strucvar/g1k.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvar/g1k.bed.gz.tbi.md5",
    shell:
        r"""
        set -x

        wget --no-check-certificate \
            -O {output.vcf} \
            https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/integrated_sv_map/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz

        zcat {output.vcf} \
        | awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-g1k.awk \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule vardbs_grch37_struc_vars_gnomad_sv:
    output:
        vcf="vardbs/grch37/strucvar/download/gnomad_v2.1_sv.sites.vcf.gz",
        bed="vardbs/grch37/strucvar/gnomad_sv.bed.gz",
        bed_md5="vardbs/grch37/strucvar/gnomad_sv.bed.gz.md5",
        bed_tbi="vardbs/grch37/strucvar/gnomad_sv.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvar/gnomad_sv.bed.gz.tbi.md5",
    shell:
        r"""
        set -x

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        wget --no-check-certificate \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz

        echo -e "#chromosome\tbegin\tend\tsv_type\tn_homalt\tn_het" \
        > $TMPDIR/tmp.bed

        bcftools query \
            -e 'SVTYPE="MCNV"' \
            -f "%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\t%INFO/N_HOMALT\t%INFO/N_HET\n" \
            {output.vcf} \
        >> $TMPDIR/tmp.bed
        bgzip -c $TMPDIR/tmp.bed >{output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
