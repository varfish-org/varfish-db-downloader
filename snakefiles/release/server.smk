def input_result_grch37_release_server_db(wildcards):
    dbsnp = expand("GRCh37/dbSNP/b155/Dbsnp.{chrom}.release_info", chrom=CHROMS)
    gnomad_exomes = expand(
        "GRCh37/gnomAD_exomes/r2.1.1/GnomadExomes.{chrom}.release_info", chrom=CHROMS
    )
    gnomad_genomes = expand(
        "GRCh37/gnomAD_genomes/r2.1.1/GnomadGenomes.{chrom}.release_info", chrom=CHROMS_NO_Y
    )
    other = [
        "GRCh37/clinvar/{download_date}/Clinvar.release_info",
        "GRCh37/dbVar/{download_date}/DbVarSv.release_info",
        "GRCh37/DGV/2016/DgvGoldStandardSvs.release_info",
        "GRCh37/DGV/2020/DgvSvs.release_info",
        "GRCh37/ensembl_genes/r104/GeneInterval.release_info",
        "GRCh37/ensembl_regulatory/{download_date}/EnsemblRegulatoryFeature.release_info",
        "GRCh37/ensembltogenesymbol/{download_date}/EnsemblToGeneSymbol.release_info",
        "GRCh37/ensembltorefseq/{download_date}/EnsemblToRefseq.release_info",
        "GRCh37/ExAC_constraints/r0.3.1/ExacConstraints.release_info",
        "GRCh37/ExAC/r1/ExacCnv.release_info",
        "GRCh37/ExAC/r1/Exac.release_info",
        "GRCh37/extra_annos/{download_date}/ExtraAnnoField.release_info",
        "GRCh37/extra_annos/{download_date}/ExtraAnno.release_info",
        "GRCh37/gnomAD_constraints/v2.1.1/GnomadConstraints.release_info",
        "GRCh37/gnomAD_SV/v2.1/GnomAdSv.release_info",
        "GRCh37/HelixMTdb/20200327/HelixMtDb.release_info",
        "GRCh37/hgmd_public/ensembl_r104/HgmdPublicLocus.release_info",
        "GRCh37/hgnc/{download_date}/Hgnc.release_info",
        "GRCh37/hgnc/{download_date}/RefseqToHgnc.release_info",
        "GRCh37/knowngeneaa/{download_date}/KnowngeneAA.release_info",
        "GRCh37/MITOMAP/{download_date}/Mitomap.release_info",
        "GRCh37/mtDB/{download_date}/MtDb.release_info",
        "GRCh37/refseq_genes/r105/GeneInterval.release_info",
        "GRCh37/tads_hesc/dixon2012/TadInterval.release_info",
        "GRCh37/tads_imr90/dixon2012/TadInterval.release_info",
        "GRCh37/thousand_genomes/phase3/ThousandGenomes.release_info",
        "GRCh37/thousand_genomes/phase3/ThousandGenomesSv.release_info",
        "GRCh37/vista/{download_date}/VistaEnhancer.release_info",
    ]
    all = dbsnp + gnomad_exomes + gnomad_genomes + other
    return [s.format(**config) for s in all]


rule result_grch37_release_server_db:
    input:
        input_result_grch37_release_server_db,
    output:
        tar="releases/{release_name}/varfish-server-background-db-{release_name}-grch37.tar.gz".format(
            **config
        ),
        sha256="releases/{release_name}/varfish-server-background-db-{release_name}-grch37.tar.gz.sha256".format(
            **config
        ),
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        mkdir -p $TMPDIR/varfish-server-background-db-{config[release_name]}-grch37

        import_versions=$TMPDIR/varfish-server-background-db-{config[release_name]}-grch37/import_versions.tsv

        echo -e "build\ttable_group\tversion" > $import_versions

        for path in {input}; do
            genome=$(echo $path | cut -d / -f 1)
            db=$(echo $path | cut -d / -f 2)
            version=$(echo $path | cut -d / -f 3)
            genome=$(echo $path | cut -d / -f 1)
            abs=$(readlink -f $path)

            echo -e "${{genome}}\t${{db}}\t${{version}}" >> $import_versions

            ( \
                cd $TMPDIR/varfish-server-background-db-{config[release_name]}-grch37; \
                mkdir -p $genome/$db/$version; \
                ln -s $abs $path; \
                ln -s ${{abs%.release_info}}.tsv ${{path%.release_info}}.tsv \
            )
        done

        tar \
            -C $TMPDIR \
            -vhczf $(readlink -f {output.tar}) \
            varfish-server-background-db-{config[release_name]}-grch37
        pushd $(dirname {output.tar})
        sha256sum $(basename {output.tar}) >$(basename {output.tar}).sha256
        """


def input_result_grch38_release_server_db(wildcards):
    dbsnp = expand("GRCh38/dbSNP/b155/Dbsnp.{chrom}.release_info", chrom=CHROMS)
    gnomad_exomes = expand(
        "GRCh38/gnomAD_exomes/r2.1.1/GnomadExomes.{chrom}.release_info", chrom=CHROMS
    )
    gnomad_genomes = expand(
        "GRCh38/gnomAD_genomes/r3.1.1/GnomadGenomes.{chrom}.release_info", chrom=CHROMS_NO_Y
    )

    other = [
        "GRCh38/clinvar/{download_date}/Clinvar.release_info",
        "GRCh38/dbVar/{download_date}/DbVarSv.release_info",
        "GRCh38/DGV/2020/DgvSvs.release_info",
        "GRCh38/DGV/2020/DgvGoldStandardSvs.release_info",
        "GRCh38/ensembl_genes/r104/GeneInterval.release_info",
        "GRCh38/ensembl_regulatory/{download_date}/EnsemblRegulatoryFeature.release_info",
        "GRCh38/ensembltogenesymbol/{download_date}/EnsemblToGeneSymbol.release_info",
        "GRCh38/ensembltorefseq/{download_date}/EnsemblToRefseq.release_info",
        "GRCh38/extra_annos/{download_date}/ExtraAnnoField.release_info",
        "GRCh38/extra_annos/{download_date}/ExtraAnno.release_info",
        "GRCh38/gnomAD_constraints/v2.1.1/GnomadConstraints.release_info",
        "GRCh38/HelixMTdb/20200327/HelixMtDb.release_info",
        "GRCh38/hgmd_public/ensembl_r104/HgmdPublicLocus.release_info",
        "GRCh38/hgnc/{download_date}/Hgnc.release_info",
        "GRCh38/hgnc/{download_date}/RefseqToHgnc.release_info",
        "GRCh38/knowngeneaa/{download_date}/KnowngeneAA.release_info",
        "GRCh38/MITOMAP/{download_date}/Mitomap.release_info",
        "GRCh38/mtDB/{download_date}/MtDb.release_info",
        "GRCh38/refseq_genes/r39/GeneInterval.release_info",
    ]
    all = dbsnp + gnomad_exomes + gnomad_genomes + other
    return [s.format(**config) for s in all]


rule result_grch38_release_server_db:
    input:
        input_result_grch38_release_server_db,
    output:
        tar="releases/{release_name}/varfish-server-background-db-{release_name}-grch38.tar.gz".format(
            **config
        ),
        sha256="releases/{release_name}/varfish-server-background-db-{release_name}-grch38.tar.gz.sha256".format(
            **config
        ),
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        mkdir -p $TMPDIR/varfish-server-background-db-{config[release_name]}-grch38

        import_versions=$TMPDIR/varfish-server-background-db-{config[release_name]}-grch38/import_versions.tsv

        echo -e "build\ttable_group\tversion" > $import_versions

        for path in {input}; do
            genome=$(echo $path | cut -d / -f 1)
            db=$(echo $path | cut -d / -f 2)
            version=$(echo $path | cut -d / -f 3)
            genome=$(echo $path | cut -d / -f 1)
            abs=$(readlink -f $path)

            echo -e "${{genome}}\t${{db}}\t${{version}}" >> $import_versions

            ( \
                cd $TMPDIR/varfish-server-background-db-{config[release_name]}-grch38; \
                mkdir -p $genome/$db/$version; \
                ln -s $abs $path; \
                ln -s ${{abs%.release_info}}.tsv ${{path%.release_info}}.tsv \
            )
        done

        tar \
            -C $TMPDIR \
            -vhczf $(readlink -f {output.tar}) \
            varfish-server-background-db-{config[release_name]}-grch38
        pushd $(dirname {output.tar})
        sha256sum $(basename {output.tar}) >$(basename {output.tar}).sha256
        """
