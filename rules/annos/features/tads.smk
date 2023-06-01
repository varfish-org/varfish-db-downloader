## Rules related to TAD.


rule annos_features_tads_grch37_download:  # -- download GRCh37 TAD files
    output:
        download_imr90="features/grch37/tads/download/IMR90_domains_hg19.bed",
        download_hesc="features/grch37/tads/download/hESC_domains_hg19.bed",
    shell:
        r"""
        wget --no-check-certificate \
            -O {input.download_imr90} \
            https://compbio.med.harvard.edu/modencode/webpage/hic/IMR90_domains_hg19.bed

        wget --no-check-certificate \
            -O {input.download_hesc} \
            https://compbio.med.harvard.edu/modencode/webpage/hic/hESC_domains_hg19.bed
        """


rule annos_features_tads_grch37_process:  # -- process GRCh37 TAD files
    input:
        download_imr90="features/grch37/tads/download/IMR90_domains_hg19.bed",
        download_hesc="features/grch37/tads/download/hESC_domains_hg19.bed",
    output:
        bed_imr90="features/grch37/tads/imr90.bed",
        bed_imr90_md5="features/grch37/tads/imr90.bed.md5",
        bed_hesc="features/grch37/tads/hesc.bed",
        bed_hesc_md5="features/grch37/tads/hesc.bed.md5",
    shell:
        r"""
        echo -e "#chrom\tbegin\tend" >{output.bed_imr90}
        sed -e 's/^chr//' {output.download_imr90} \
        >>{output.bed_imr90}

        echo -e "#chrom\tbegin\tend" >{output.bed_hesc}
        sed -e 's/^chr//' {output.download_hesc} \
        >>{output.bed_hesc}

        md5sum {output.bed_imr90} >{output.bed_imr90_md5}
        md5sum {output.bed_hesc} >{output.bed_hesc_md5}
        """
