rule grchxx_ensembl_regulatory_download:
    output:
        tsv="{genome_build}/ensembl_regulatory/{download_date}/download/hsapiens_regulatory_feature__regulatory_feature__main.txt.gz",
    shell:
        r"""
        if [[ {wildcards.genome_build} == GRCh37 ]]; then
            infix=/grch37
        else
            infix=
        fi

        wget -O {output.tsv} http://ftp.ensembl.org/pub${{infix}}/release-104/mysql/regulation_mart_104/hsapiens_regulatory_feature__regulatory_feature__main.txt.gz
        """


rule grchxx_ensembl_regulatory_merge:
    input:
        tsv="{genome_build}/ensembl_regulatory/{download_date}/download/hsapiens_regulatory_feature__regulatory_feature__main.txt.gz",
    output:
        tsv="{genome_build}/ensembl_regulatory/{download_date}/download/EnsemblRegulatoryFeature.tsv",
    shell:
        r"""
        if [[ {wildcards.genome_build} == GRCh37 ]]; then
            chr=
        else
            chr=chr
        fi

        echo -e "Chromosome/scaffold name\tStart (bp)\tEnd (bp)\tRegulatory stable ID\tFeature type\tFeature type description\tSO term accession\tSO term name" \
        > {output.tsv}

        zcat {input} \
        | awk -F $'\t' 'BEGIN {{ OFS = FS }} {{ print $6, $10, $11, $7, $3, $1, $5, $8 }}' \
        | grep '^[0-9XYM]' \
        | LC_ALL=C sort -k1,1g -k2,2n -k3,3n \
        | sed -e "s/^/$chr/" \
        >> {output.tsv}
        """


rule result_grchxx_ensembl_regulatory_tsv:
    input:
        tsv="{genome_build}/ensembl_regulatory/{download_date}/download/EnsemblRegulatoryFeature.tsv",
        header="header/ensembl_regulatory.txt",
    output:
        tsv="{genome_build}/ensembl_regulatory/{download_date}/EnsemblRegulatoryFeature.tsv",
        release_info="{genome_build}/ensembl_regulatory/{download_date}/EnsemblRegulatoryFeature.release_info",
    wildcard_constraints:
        download_date="[^/]+",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.tsv} \
            | awk -F$"\t" 'BEGIN{{OFS=FS}}{{$1="{wildcards.genome_build}\t"$1; $3=$3"\t"; print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblRegulatoryFeature\t{wildcards.download_date}\t{wildcards.genome_build}\t" > {output.release_info}
        """
