## Download of clinvar files for worker.


rule annos_strucvars_clinvar_download:  # -- download/extract ClinVar files
    output:
        tar=f"work/download/annos/{{genome_release}}/strucvars/clinvar/{{clinvar_version}}/clinvar-strucvar-{{genome_release}}-{DV.clinvar_release}.tar.gz",
        tsv="work/download/annos/{genome_release}/strucvars/clinvar/{clinvar_version}/clinvar_strucvar.tsv.gz",
    wildcard_constraints:
        genome_release=RE_GENOME,
        clinvar_version=RE_VERSION,
    shell:
        r"""
        clinvar_version=$(echo "{wildcards.clinvar_version}" | sed -e 's/-//g' | cut -d '+' -f 1)

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        wget --no-check-certificate \
            -O {output.tar} \
            https://github.com/bihealth/annonars-data-clinvar/releases/download/clinvar-weekly-$clinvar_version/$(basename {output.tar})

        if [[ {wildcards.genome_release} == grch37 ]]; then
            release=GRCh37
        else
            release=GRCh38
        fi

        tar -C $TMPDIR -xvf $(readlink -f {output.tar})
        cut -f 2- $TMPDIR/clinvar-strucvar-{wildcards.genome_release}*/output.tsv \
        | awk \
            -F $'\t' \
            -v release=$release \
            -v clinvar_version={wildcards.clinvar_version} \
            \
            '
            BEGIN {{ OFS=FS }}
            (NR == 1) {{ print "#" $0; }}
            (NR > 1) {{ $8 = clinvar_version; print $0; }}
            ' \
        | gzip -c \
        > {output.tsv}
        """
