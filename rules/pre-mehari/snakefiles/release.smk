NOREFS = [
    f"output/pre-mehari/noref/hpo/{DV.hpo}/Hpo.release_info",
    f"output/pre-mehari/noref/hpo/{DV.hpo}/HpoName.release_info",
    f"output/pre-mehari/noref/acmg/{DV.acmg_sf}/Acmg.release_info",
    f"output/pre-mehari/noref/mim2gene/{DV.today}/Mim2geneMedgen.release_info",
    f"output/pre-mehari/noref/refseqtogenesymbol/{DV.today}/RefseqToGeneSymbol.release_info",
]


def input_result_grch3x_release_server_db(wildcards):
    genomebuild = wildcards.genomebuild
    dbsnp = expand(
        "output/pre-mehari/{genomebuild}/dbSNP/{dbsnp}/Dbsnp.{chrom}.release_info",
        genomebuild=genomebuild_cap[genomebuild],
        dbsnp=DV.dbsnp,
        chrom=chroms[genomebuild],
    )
    return (
        NOREFS
        + dbsnp
        + [
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/hgnc/{DV.hgnc_quarterly}+{DV.cdot}+{refseq_versions[genomebuild]}/Hgnc.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/hgnc/{DV.hgnc_quarterly}+{DV.cdot}+{refseq_versions[genomebuild]}/RefseqToHgnc.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/clinvar/{DV.hgnc_quarterly}+{DV.clinvar_release}/Clinvar.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/HelixMTdb/{DV.helixmtdb}/HelixMtDb.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/ensembltogenesymbol/{ensembl_versions[genomebuild]}/EnsemblToGeneSymbol.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/MITOMAP/{DV.today}/Mitomap.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/mtDB/{DV.mtdb}/MtDb.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/extra_annos/{DV.cadd}/ExtraAnno.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/extra_annos/{DV.cadd}/ExtraAnnoField.release_info",
            f"output/pre-mehari/{genomebuild_cap[genomebuild]}/knowngeneaa/{cons_versions[genomebuild]}/KnowngeneAA.release_info",
            # f"output/pre-mehari/{genomebuild_cap[genomebuild]}/gnomAD_constraints/v{gnomad_versions[genomebuild]}/GnomadConstraints.release_info",
        ]
    )


rule result_grch3x_release_server_db:
    input:
        input_result_grch3x_release_server_db,
    output:
        tmp_import_versions=temp("output/pre-mehari/releases/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}/.import_versions.tsv"),
        import_versions="output/pre-mehari/releases/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}/import_versions.tsv",
    shell:
        r"""
        out_dir=$(dirname {output.import_versions})
        mkdir -p $out_dir

        for path in {input}; do
            genome=$(echo $path | cut -d / -f 3)
            db=$(echo $path | cut -d / -f 4)
            version=$(echo $path | cut -d / -f 5)
            abs=$(readlink -f $path)
            dir=$genome/$db/$version
            file=$(basename $path)
            dirfile=$dir/$file

            echo -e "${{genome}}\t${{db}}\t${{version}}" >> {output.tmp_import_versions}

            ( \
                cd $out_dir; \
                mkdir -p $dir; \
                test -e $dirfile || ln -rs $abs $dirfile; \
                test -e ${{dirfile%.release_info}}.tsv || ln -rs ${{abs%.release_info}}.tsv ${{dirfile%.release_info}}.tsv \
            )
        done

        (
            echo -e "build\ttable_group\tversion"
            sort -u {output.tmp_import_versions}
        ) > {output.import_versions}
        """


rule result_grch3x_release_manifest:
    input:
        done="output/pre-mehari/releases/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}/import_versions.tsv",
    output:
        manifest="output/full/pre-mehari/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}/manifest-postgres.json",
    shell:
        r"""
        release_dir=$(dirname {input.done})

        mkdir -p $(dirname {output.manifest})

        python rules/pre-mehari/snakefiles/scripts/generate_manifest.py \
            $release_dir \
            {output.manifest}
        """


rule result_grch3x_release_server_db_tar:
    input:
        "output/pre-mehari/releases/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}/import_versions.tsv",
        "output/full/pre-mehari/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}/manifest-postgres.json",
    output:
        tar="output/pre-mehari/releases/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}.tar.gz",
        sha256="output/pre-mehari/releases/{release_name}/varfish-postgres-db-{release_name}-{genomebuild}.tar.gz.sha256",
    shell:
        r"""
        in_dir=$(dirname {input})

        tar \
            --owner=0 \
            --group=0 \
            -C $in_dir/.. \
            -vhcf - \
            varfish-postgres-db-{wildcards.release_name}-{wildcards.genomebuild} \
        | pigz \
        > $(readlink -f {output.tar})
        pushd $(dirname {output.tar})
        sha256sum $(basename {output.tar}) >$(basename {output.tar}).sha256
        """
