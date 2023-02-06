rule GRCh37_vista_download:
    output:
        tsv="GRCh37/vista/{download_date}/download/VistaEnhancer.tsv",
    shell:
        r"""
        2>/dev/null wget --no-check-certificate -O - 'https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page_size=100;show=1;search.result=yes;search.form=no;form=search;action=search;search.org=Both;search.sequence=1' \
        | sed -e 's/<pre>//g' \
        | sed -e 's/<\/pre>//g' \
        | grep '^>' \
        | sed -e 's/element /hs/g' \
        | grep Human \
        | sed -e 's/[:\|-]/ /g' \
        | awk 'BEGIN {{ OFS = "\t" }} {{ print $2, $3, $4, $5, $6 }}' \
        | cut -b 4- \
        | LC_ALL=C sort -k1,1g -k2,2n -k3,3n \
        > {output.tsv}
        """


rule result_GRCh37_vista_tsv:
    input:
        tsv="GRCh37/vista/{download_date}/download/VistaEnhancer.tsv",
        header="header/vista.txt",
    output:
        tsv="GRCh37/vista/{download_date}/VistaEnhancer.tsv",
        release_info="GRCh37/vista/{download_date}/VistaEnhancer.release_info",
    wildcard_constraints:
        download_date="[^/]+",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="GRCh37\t"$1; $3=$3"\t"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nVistaEnhancer\t{wildcards.download_date}\tGRCh37\t" > {output.release_info}
        """
