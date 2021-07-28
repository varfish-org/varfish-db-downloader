# Build ClinVar TSV files.


# Run clinvar_tsv tool to build for {clinvar_release} clinvar.
rule grch37_clinvar_run_clinvar_tsv:
    input:
        hs37d5="GRCh37/reference/hs37d5/hs37d5.fa",
        hs38="GRCh38/reference/hs38/hs38.fa",
    output:
        b38="GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/output/clinvar.b38.tsv.gz",
        b38_tbi="GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/output/clinvar.b38.tsv.gz.tbi",
        b37="GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/output/clinvar.b37.tsv.gz",
        b37_tbi="GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/output/clinvar.b37.tsv.gz.tbi",
        done=[
            touch("GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/.done"),
            touch("GRCh38/clinvar/{clinvar_release}/clinvar_tsv_main/.done"),
        ],
    log:
        "GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/clinvar_tsv.log",
    shell:
        r"""
        set -x

        clinvar_tsv \
            main \
            --b37-path $(readlink -f {input.hs37d5}) \
            --b38-path $(readlink -f {input.hs38}) \
            --work-dir $(dirname {output.done[0]}) \
        2>&1 | tee {log}
        """


# Copy TSV out to appropriate position, create release_info files.
rule grchXX_clinvar_postprocess:
    input:
        tsv="GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/output/clinvar.b{release}.tsv.gz",
        done="GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/.done",
    output:
        tsv="GRCh{release}/clinvar/{clinvar_release}/Clinvar.tsv",
        md5="GRCh{release}/clinvar/{clinvar_release}/Clinvar.tsv.md5",
        release_info="GRCh{release}/clinvar/{clinvar_release}/Clinvar.release_info",
    shell:
        r"""
        set -x

        zcat {input.tsv} >{output.tsv}

        pushd $(dirname {output.tsv})
        md5sum $(basename {output.tsv}) >$(basename {output.tsv}).md5
        popd

        echo -e "table\tversion\tgenomebuild\tnull_value\nClinvar\t$(date +%Y/%m/%d)\tGRCh{wildcards.release}\t" > {output.release_info}
        """
