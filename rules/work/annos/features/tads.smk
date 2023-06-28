## Rules related to TAD.


rule annos_features_tads_download:  # -- download TAD ZIP files from 3dgenome.org
    output:
        zip="work/download/annos/{genome_release}/tads/{genome_release}.TADs.zip",
    shell:
        r"""
        if [[ "{wildcards.genome_release}" == grch37 ]]; then
            name=hg19
        else
            name=hg38
        fi

        wget --no-check-certificate \
            -O {output.zip} \
            http://3dgenome.fsm.northwestern.edu/downloads/$name.TADs.zip
        """


rule annos_features_tads_unzip:  # -- unzip TAD ZIP files
    input:
        zip="work/download/annos/{genome_release}/tads/{genome_release}.TADs.zip",
    output:
        hesc="work/download/annos/{genome_release}/tads/H1-ESC_Dixon2015-raw_TADs.txt",
    shell:
        r"""
        unzip -o -d $(dirname {input.zip}) -j {input.zip}

        cd $(dirname {input.zip})

        if $(set +e; ls | grep Dixon_2015 >/dev/null); then
            for f in *Dixon_2015*; do
                mv $f ${{f/Dixon_/Dixon}}
            done
        fi
        """


rule annos_features_tads_process:  # -- process TAD files
    input:
        hesc="work/download/annos/{genome_release}/tads/H1-ESC_Dixon2015-raw_TADs.txt",
    output:
        bed_hesc="output/full/tracks/track-tads-{genome_release}-dixon2015+{DV.tracks}/hesc.bed",
        bed_hesc_md5="output/full/tracks/track-tads-{genome_release}-dixon2015+{DV.tracks}/hesc.bed.md5",
    shell:
        r"""
        echo -e "#chrom\tbegin\tend" >{output.bed_hesc}
        cat {input.hesc} \
        | {{ if [[ "{wildcards.genome_release}" == "grch37" ]]; then sed -e 's/^chr//g'; else cat; fi }} \
        >>{output.bed_hesc}

        md5sum {output.bed_hesc} >{output.bed_hesc_md5}
        """
