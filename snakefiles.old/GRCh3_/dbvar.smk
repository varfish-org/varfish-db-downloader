rule grchXX_dbvar_latest_download:
    output:
        expand(
            "{{genome_build}}/dbVar/{{download_date}}/download/{{type}}/{{genome_build}}.nr_{{type}}.{ending}",
            ending=["bed.gz", "bedpe.gz", "tsv.gz", "acmg_genes.tsv.gz"],
        ),
    shell:
        r"""
        cd $(dirname $(dirname {output[0]}))

        # mirror server
        echo 'mirror -P 8 --include="{wildcards.genome_build}"' | lftp http://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/

        # set rights as the image has no execute rights to access the folders.
        chmod ug+wX -R .

        # check md5 sums
        find -iname '*.md5' | xargs -I+ bash -c 'cd $(dirname +) && md5sum -c $(basename +)'
        """


rule grchXX_dbvar_latest_tsv:
    input:
        header="header/dbvarsv.txt",
        tsv="{genome_build}/dbVar/{download_date}/download/{type}/{genome_build}.nr_{type}.tsv.gz",
    output:
        tsv=temp("{genome_build}/dbVar/{download_date}/DbVarSv:{type}.tsv"),
        release_info=temp("{genome_build}/dbVar/{download_date}/DbVarSv:{type}.release_info"),
    run:
        to_tsv(input.tsv, output.tsv, output.release_info, input.header)


rule result_grchXX_dbvar_merge_tsv:
    input:
        tsv_del="{genome_build}/dbVar/{download_date}/DbVarSv:deletions.tsv",
        tsv_dup="{genome_build}/dbVar/{download_date}/DbVarSv:duplications.tsv",
        tsv_ins="{genome_build}/dbVar/{download_date}/DbVarSv:insertions.tsv",
        release_info_del="{genome_build}/dbVar/{download_date}/DbVarSv:deletions.release_info",
        release_info_dup="{genome_build}/dbVar/{download_date}/DbVarSv:duplications.release_info",
        release_info_ins="{genome_build}/dbVar/{download_date}/DbVarSv:insertions.release_info",
    output:
        tsv="{genome_build}/dbVar/{download_date}/DbVarSv.tsv",
        release_info="{genome_build}/dbVar/{download_date}/DbVarSv.release_info",
    shell:
        r"""
        (
            cat {input.tsv_del}
            tail -n +2 {input.tsv_dup}
            tail -n +2 {input.tsv_ins}
        ) > {output.tsv}

        sed 's/:deletions//' {input.release_info_del} > {output.release_info}
        """
