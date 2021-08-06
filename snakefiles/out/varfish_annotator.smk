"""
    init-db      Initialize or update DB
      Usage: init-db [options]
        Options:
          --clinvar-path
            Path to Clinvar TSV file(s) to use for import, see documentation 
            for more information
        * --db-path
            Path to H2 file to initialize/update
          --db-release-info
            Provide database release information as "$db:$release" for storage 
            in DB
          --exac-path
            Path to ExAC VCF file to use for import, see documentation for 
            more information
          --gnomad-exomes-path
            Path to gnomAD exomes VCF file to use for import, see 
            documentation for more information
          --gnomad-genomes-path
            Path to gnomAD genomes VCF file to use for import, see 
            documentation for more information
          --help

          --hgmd-public
            Path to HTMD Public TSV file to use for import, see documentation 
            for more information
        * --ref-path
            Path to reference FASTA file, used for variant normalization
          --region
            Genomic region CHR:START-END (1-based) to import
          --thousand-genomes-path
            Path to 1000 genomes VCF file to use for import, see documentation 
            for more information
"""


rule result_grch37_varfish_annotator_db:
    input:
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        clinvar="GRCh37/clinvar/{download_date}/clinvar_tsv_main/output/clinvar.b37.tsv.gz".format(
            **config
        ),
        gnomad_exomes_chr1="GRCh37/gnomAD_exomes/r2.1.1/download/gnomad.exomes.r2.1.1.sites.chr1.normalized.vcf.bgz",
        gnomad_genomes_chr1="GRCh37/gnomAD_genomes/r2.1.1/download/gnomad.genomes.r2.1.1.sites.chr1.normalized.vcf.bgz",
        exac="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz",
        thousand_genomes="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
        hgmd="GRCh37/hgmd_public/ensembl_r104/HgmdPublicLocus.tsv",
    output:
        "releases/{release_name}/varfish-annotator-db-grch37-{release_name}.h2.db".format(**config),
    shell:
        r"""
        _out={output}
        out=$PWD/${{_out%.h2.db}}
        _gnomad_exomes_prefix={input.gnomad_exomes_chr1}
        gnomad_exomes_prefix=${{_gnomad_exomes_prefix%.chr1.normalized.vcf.bgz}}
        _gnomad_genomes_prefix={input.gnomad_genomes_chr1}
        gnomad_genomes_prefix=${{_gnomad_genomes_prefix%.chr1.normalized.vcf.bgz}}

        varfish-annotator init-db \
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
            $(for path in ${{gnomad_exomes_prefix}}.*.normalized.vcf.bgz; do \
                echo --gnomad-exomes-path $path; \
            done) \
            \
            --db-release-info "gnomad_genomes:r2.1.1" \
            $(for path in ${{gnomad_genomes_prefix}}.*.normalized.vcf.bgz; do \
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
        gnomad_exomes_chr1="GRCh38/gnomAD_exomes/r2.1.1/download/gnomad.exomes.r2.1.1.sites.chr1.normalized.vcf.bgz",
        gnomad_genomes_chr1="GRCh38/gnomAD_genomes/r3.1.1/download/gnomad.genomes.r3.1.1.sites.chr1.normalized.vcf.bgz",
        hgmd="GRCh38/hgmd_public/ensembl_r104/HgmdPublicLocus.tsv",
    output:
        "releases/{release_name}/varfish-annotator-db-grch38-{release_name}.h2.db".format(**config),
    shell:
        r"""
        _out={output}
        out=$PWD/${{_out%.h2.db}}
        _gnomad_exomes_prefix={input.gnomad_exomes_chr1}
        gnomad_exomes_prefix=${{_gnomad_exomes_prefix%.chr1.normalized.vcf.bgz}}
        _gnomad_genomes_prefix={input.gnomad_genomes_chr1}
        gnomad_genomes_prefix=${{_gnomad_genomes_prefix%.chr1.normalized.vcf.bgz}}

        varfish-annotator init-db \
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
            $(for path in ${{gnomad_exomes_prefix}}.*.normalized.vcf.bgz; do \
                echo --gnomad-exomes-path $path; \
            done) \
            \
            --db-release-info "gnomad_genomes:r3.1.1" \
            $(for path in ${{gnomad_genomes_prefix}}.*.normalized.vcf.bgz; do \
                echo --gnomad-genomes-path $path; \
            done) \
            \
            --db-release-info "hgmd_public:ensembl_r104" \
            --hgmd-public {input.hgmd}
        """
