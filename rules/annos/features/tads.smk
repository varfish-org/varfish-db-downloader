## Rules related to TAD.


rule annos_features_tads_grch37_download:  # -- download GRCh37 TAD files
    output:
        imr90="work/download/annos/grch37/tads/dixon2015/IMR90_domains_hg19.bed",
        hesc="work/download/annos/grch37/tads/dixon2015/hESC_domains_hg19.bed",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.imr90} \
            https://compbio.med.harvard.edu/modencode/webpage/hic/IMR90_domains_hg19.bed

        wget --no-check-certificate \
            -O {output.hesc} \
            https://compbio.med.harvard.edu/modencode/webpage/hic/hESC_domains_hg19.bed
        """


rule annos_features_tads_grch37_process:  # -- process GRCh37 TAD files
    input:
        imr90="work/download/annos/grch37/tads/dixon2015/IMR90_domains_hg19.bed",
        hesc="work/download/annos/grch37/tads/dixon2015/hESC_domains_hg19.bed",
    output:
        bed_imr90="work/annos/grch37/features/tads/dixon2015/imr90.bed",
        bed_imr90_md5="work/annos/grch37/features/tads/dixon2015/imr90.bed.md5",
        bed_hesc="work/annos/grch37/features/tads/dixon2015/hesc.bed",
        bed_hesc_md5="work/annos/grch37/features/tads/dixon2015/hesc.bed.md5",
    shell:
        r"""
        echo -e "#chrom\tbegin\tend" >{output.bed_imr90}
        sed -e 's/^chr//' {input.imr90} \
        >>{output.bed_imr90}

        echo -e "#chrom\tbegin\tend" >{output.bed_hesc}
        sed -e 's/^chr//' {input.hesc} \
        >>{output.bed_hesc}

        md5sum {output.bed_imr90} >{output.bed_imr90_md5}
        md5sum {output.bed_hesc} >{output.bed_hesc_md5}
        """
