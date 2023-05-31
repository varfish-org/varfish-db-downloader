rule result_grch37_varfish_annotator_db:
    input:
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        clinvar="GRCh37/clinvar/{download_date}/clinvar_tsv_main/output/clinvar.b37.tsv.gz".format(
            **config
        ),
        gnomad_exomes=expand(
            "GRCh37/gnomAD_exomes/r2.1.1/download/gnomad.exomes.r2.1.1.sites.chr{chrom}.stripped.vcf.bgz",
            chrom=CHROMS,
        ),
        gnomad_genomes=expand(
            "GRCh37/gnomAD_genomes/r2.1.1/download/gnomad.genomes.r2.1.1.sites.chr{chrom}.stripped.vcf.bgz",
            chrom=CHROMS_NO_Y,
        ),
        exac="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz",
        thousand_genomes="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
        hgmd="GRCh37/hgmd_public/ensembl_r104/HgmdPublicLocus.tsv",
    output:
        "GRCh37/varfish-annotator-db/varfish-annotator-db-{release_name}-grch37.h2.db".format(
            **config
        ),
    shell:
        r"""
        _out={output}
        out=$PWD/${{_out%.h2.db}}

        varfish-annotator init-db \
            --release GRCh37 \
            --db-release-info "varfish-annotator:v{config[annotator_version]}" \
            --db-release-info "varfish-annotator-db:r{config[release_name]}" \
            --db-path $out \
            \
            --ref-path {input.reference} \
            \
            --db-release-info "exac:r1.0" \
            --exac-path {input.exac} \
            \
            --db-release-info "thousand_genomes:v3.20101123" \
            --thousand-genomes-path {input.thousand_genomes} \
            \
            --db-release-info "clinvar:{config[download_date]}" \
            --clinvar-path {input.clinvar} \
            \
            --db-release-info "gnomad_exomes:r2.1.1" \
            $(for path in {input.gnomad_exomes}; do \
                echo --gnomad-exomes-path $path; \
            done) \
            \
            --db-release-info "gnomad_genomes:r2.1.1" \
            $(for path in {input.gnomad_genomes}; do \
                echo --gnomad-genomes-path $path; \
            done) \
            \
            --db-release-info "hgmd_public:ensembl_r104" \
            --hgmd-public {input.hgmd}
        """


rule result_grch38_varfish_annotator_db:
    input:
        reference="GRCh38/reference/hs38/hs38.fa",
        clinvar="GRCh37/clinvar/{download_date}/clinvar_tsv_main/output/clinvar.b38.tsv.gz".format(
            **config
        ),
        gnomad_exomes=expand(
            "GRCh38/gnomAD_exomes/r2.1.1/download/gnomad.exomes.r2.1.1.sites.chr{chrom}.stripped.vcf.bgz",
            chrom=CHROMS,
        ),
        gnomad_genomes=expand(
            "GRCh38/gnomAD_genomes/r3.1.1/download/gnomad.genomes.r3.1.1.sites.chr{chrom}.stripped.vcf.bgz",
            chrom=CHROMS_NO_Y,
        ),
        hgmd="GRCh38/hgmd_public/ensembl_r104/HgmdPublicLocus.tsv",
    output:
        "GRCh38/varfish-annotator-db/varfish-annotator-db-{release_name}-grch38.h2.db".format(
            **config
        ),
    shell:
        r"""
        _out={output}
        out=$PWD/${{_out%.h2.db}}

        varfish-annotator init-db \
            --release GRCh38 \
            --db-release-info "varfish-annotator:v{config[annotator_version]}" \
            --db-release-info "varfish-annotator-db:r{config[release_name]}" \
            --db-path $out \
            \
            --ref-path {input.reference} \
            \
            --db-release-info "clinvar:{config[download_date]}" \
            --clinvar-path {input.clinvar} \
            \
            --db-release-info "gnomad_exomes:r2.1.1" \
            $(for path in {input.gnomad_exomes}; do \
                echo --gnomad-exomes-path $path; \
            done) \
            \
            --db-release-info "gnomad_genomes:r3.1.1" \
            $(for path in {input.gnomad_genomes}; do \
                echo --gnomad-genomes-path $path; \
            done) \
            \
            --db-release-info "hgmd_public:ensembl_r104" \
            --hgmd-public {input.hgmd}
        """
