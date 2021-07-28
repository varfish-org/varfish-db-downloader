rule GRCh37_tads_hesc_download:
    output:
        tsv="GRCh37/tads_hesc/dixon2012/download/hESC_domains_hg19.bed",
    shell:
        r"""
        base_url="http://compbio.med.harvard.edu/modencode/webpage/hic"

        wget -O - $base_url/hESC_domains_hg19.bed \
        | cut -b 4- \
        | sed 's/\r$//' \
        > {output.tsv}
        """


rule GRCh37_tads_hesc_tsv:
    input:
        header="header/tads.txt",
        tsv="GRCh37/tads_hesc/dixon2012/download/hESC_domains_hg19.bed",
    output:
        tsv="GRCh37/tads_hesc/dixon2012/TadInterval.tsv",
        release_info="GRCh37/tads_hesc/dixon2012/TadInterval.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="GRCh37\t"$1; $3=$3"\t"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nTadInterval:hesc\tDixon2012\tGRCh37\t" > {output.release_info}
        """


rule GRCh37_tads_imr90_download:
    output:
        tsv="GRCh37/tads_imr90/dixon2012/download/hESC_domains_hg19.bed",
    shell:
        r"""
        base_url="http://compbio.med.harvard.edu/modencode/webpage/hic"

        wget -O - $base_url/hESC_domains_hg19.bed \
        | cut -b 4- \
        | sed 's/\r$//' \
        > {output.tsv}
        """


rule GRCh37_tads_imr90_tsv:
    input:
        header="header/tads.txt",
        tsv="GRCh37/tads_imr90/dixon2012/download/hESC_domains_hg19.bed",
    output:
        tsv="GRCh37/tads_imr90/dixon2012/TadInterval.tsv",
        release_info="GRCh37/tads_imr90/dixon2012/TadInterval.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="GRCh37\t"$1; $3=$3"\t"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nTadInterval:imr90\tDixon2012\tGRCh37\t" > {output.release_info}
        """
