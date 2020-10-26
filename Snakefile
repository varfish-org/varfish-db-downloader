import os

from snakemake import shell
from tools.sv_db_to_tsv import to_tsv
from datetime import date


shell.prefix("set -x; set -euo pipefail; ")


CHROMS = (list(range(1, 23)) + ['X'])

RELEASE_PATH = os.getenv("RELEASE_PATH", "releases")

DATA_RELEASE = os.getenv("DATA_RELEASE", date.today().strftime("%Y%m%d"))
CLINVAR_RELEASE = os.getenv("CLINVAR_RELEASE", date.today().strftime("%Y%m%d"))
JANNOVAR_RELEASE = os.getenv("JANNOVAR_RELEASE", date.today().strftime("%Y%m%d"))

VARFISH_SERVER_BACKGROUND_PATH = expand("{release_path}/varfish-server-background-db-{date}", release_path=RELEASE_PATH, date=DATA_RELEASE)
VARFISH_ANNOTATOR_DB_PATH = expand("{release_path}/varfish-annotator-db-{date}", release_path=RELEASE_PATH, date=DATA_RELEASE)
VARFISH_ANNOTATOR_PATH = expand("{release_path}/varfish-annotator-{date}", release_path=RELEASE_PATH, date=DATA_RELEASE)
JANNOVAR_DB_PATH = expand("{release_path}/jannovar-db-{date}", release_path=RELEASE_PATH, date=JANNOVAR_RELEASE)

VARFISH_ANNOTATOR_DB_FILE = expand("{varfish_annotator_db}.h2.db", varfish_annotator_db=VARFISH_ANNOTATOR_DB_PATH)
JANNOVAR_FILES = expand("{jannovar_db_path}/hg19_{database}.ser", jannovar_db_path=JANNOVAR_DB_PATH, database=["ensembl", "refseq_curated"])
VARFISH_ANNOTATOR_FILES = VARFISH_ANNOTATOR_DB_FILE + JANNOVAR_FILES + ["GRCh37/reference/hs37d5/hs37d5.fa", "GRCh37/reference/hs37d5/hs37d5.fa.fai"]
VARFISH_ANNOTATOR_FILES_OUT = [os.path.basename(file) for file in VARFISH_ANNOTATOR_FILES]
VARFISH_ANNOTATOR_H2_FILES = [
    *expand("GRCh37/clinvar/{clinvar_release}/clinvar_tsv_main/output/clinvar.b37.tsv.gz{index}", clinvar_release=CLINVAR_RELEASE, index=["", ".tbi"]),
    *expand("GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz{index}", index=["", ".tbi"]),
    *expand("GRCh37/gnomAD_exomes/r2.1/download/gnomad.exomes.r2.1.sites.chr{chrom}.normalized.vcf.bgz{index}", chrom=CHROMS, index=["", ".tbi"]),
    *expand("GRCh37/gnomAD_genomes/r2.1/download/gnomad.genomes.r2.1.sites.chr{chrom}.normalized.vcf.bgz{index}", chrom=CHROMS, index=["", ".tbi"]),
    *expand("GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz{index}", index=["", ".tbi"]),
    "GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.tsv",
    "GRCh37/reference/hs37d5/hs37d5.fa",
    "GRCh37/reference/hs37d5/hs37d5.fa.fai",
]
VARFISH_FILES = [
    # clinvar
    "GRCh37/clinvar/{}/Clinvar.tsv".format(CLINVAR_RELEASE),
    "GRCh37/clinvar/{}/Clinvar.release_info".format(CLINVAR_RELEASE),
    "GRCh38/clinvar/{}/Clinvar.tsv".format(CLINVAR_RELEASE),
    "GRCh38/clinvar/{}/Clinvar.release_info".format(CLINVAR_RELEASE),
    # dbsnp
    *expand("GRCh37/dbSNP/b151/Dbsnp.{chrom}.tsv", chrom=CHROMS+['Y', 'MT']),
    *expand("GRCh37/dbSNP/b151/Dbsnp.{chrom}.release_info", chrom=CHROMS+['Y', 'MT']),
    # exac
    "GRCh37/ExAC/r1/Exac.tsv",
    "GRCh37/ExAC/r1/Exac.release_info",
    # gnomad
    *expand("GRCh37/gnomAD_exomes/r2.1/GnomadExomes.{chrom}.tsv", chrom=CHROMS),
    *expand("GRCh37/gnomAD_exomes/r2.1/GnomadExomes.{chrom}.release_info", chrom=CHROMS),
    *expand("GRCh37/gnomAD_genomes/r2.1/GnomadGenomes.{chrom}.tsv", chrom=CHROMS),
    *expand("GRCh37/gnomAD_genomes/r2.1/GnomadGenomes.{chrom}.release_info", chrom=CHROMS),
    # hgmd
    "GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.tsv",
    "GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.release_info",
    # hgnc
    "GRCh37/hgnc/latest/Hgnc.tsv",
    "GRCh37/hgnc/latest/Hgnc.release_info",
    "GRCh37/hgnc/latest/RefseqToHgnc.tsv",
    "GRCh37/hgnc/latest/RefseqToHgnc.release_info",
    # knowngeneaa
    "GRCh37/knowngeneaa/latest/KnowngeneAA.tsv",
    "GRCh37/knowngeneaa/latest/KnowngeneAA.release_info",
    # ncbi gene
    "GRCh37/ncbi_gene/latest/NcbiGeneInfo.tsv",
    "GRCh37/ncbi_gene/latest/NcbiGeneInfo.release_info",
    "GRCh37/ncbi_gene/latest/NcbiGeneRif.tsv",
    "GRCh37/ncbi_gene/latest/NcbiGeneRif.release_info",
    # thousand genomes
    "GRCh37/thousand_genomes/phase3/ThousandGenomes.tsv",
    "GRCh37/thousand_genomes/phase3/ThousandGenomes.release_info",
    # kegg
    "GRCh37/kegg/april2011/KeggInfo.tsv",
    "GRCh37/kegg/april2011/KeggInfo.release_info",
    "GRCh37/kegg/april2011/EnsemblToKegg.tsv",
    "GRCh37/kegg/april2011/EnsemblToKegg.release_info",
    "GRCh37/kegg/april2011/RefseqToKegg.tsv",
    "GRCh37/kegg/april2011/RefseqToKegg.release_info",
    # mim2gene
    "GRCh37/mim2gene/latest/Mim2geneMedgen.tsv",
    "GRCh37/mim2gene/latest/Mim2geneMedgen.release_info",
    # hpo
    "GRCh37/hpo/latest/Hpo.tsv",
    "GRCh37/hpo/latest/Hpo.release_info",
    "GRCh37/hpo/latest/HpoName.tsv",
    "GRCh37/hpo/latest/HpoName.release_info",
    # acmg
    "GRCh37/acmg/v2.0/Acmg.tsv",
    "GRCh37/acmg/v2.0/Acmg.release_info",
    # VISTA enhancer validation results
    "GRCh37/vista/latest/VistaEnhancer.tsv",
    "GRCh37/vista/latest/VistaEnhancer.release_info",
    # ENSEMBL Regulatory Features
    "GRCh37/ensembl_regulatory/latest/EnsemblRegulatoryFeature.tsv",
    "GRCh37/ensembl_regulatory/latest/EnsemblRegulatoryFeature.release_info",
    # Dixon 2012 TAD BED files
    "GRCh37/tads_hesc/dixon2012/TadInterval.tsv",
    "GRCh37/tads_hesc/dixon2012/TadInterval.release_info",
    "GRCh37/tads_imr90/dixon2012/TadInterval.tsv",
    "GRCh37/tads_imr90/dixon2012/TadInterval.release_info",
    # Genes BED files
    "GRCh37/refseq_genes/r105/GeneInterval.tsv",
    "GRCh37/refseq_genes/r105/GeneInterval.release_info",
    "GRCh37/ensembl_genes/r96/GeneInterval.tsv",
    "GRCh37/ensembl_genes/r96/GeneInterval.release_info",
    # Gnomad Constraints
    "GRCh37/gnomAD_constraints/v2.1.1/GnomadConstraints.tsv",
    "GRCh37/gnomAD_constraints/v2.1.1/GnomadConstraints.release_info",
    # Exac Constraints
    "GRCh37/ExAC_constraints/r0.3.1/ExacConstraints.tsv",
    "GRCh37/ExAC_constraints/r0.3.1/ExacConstraints.release_info",
    # Ensembl to Refseq mapping
    "GRCh37/ensembltorefseq/latest/EnsemblToRefseq.tsv",
    "GRCh37/ensembltorefseq/latest/EnsemblToRefseq.release_info",
    # Refseq to Ensembl mapping
    "GRCh37/refseqtoensembl/latest/RefseqToEnsembl.tsv",
    "GRCh37/refseqtoensembl/latest/RefseqToEnsembl.release_info",
    # MGI mouse gene information
    "GRCh37/mgi/latest/MgiHomMouseHumanSequence.tsv",
    "GRCh37/mgi/latest/MgiHomMouseHumanSequence.release_info",
    # Gene symbol mapping
    "GRCh37/refseqtogenesymbol/latest/RefseqToGeneSymbol.tsv",
    "GRCh37/refseqtogenesymbol/latest/RefseqToGeneSymbol.release_info",
    "GRCh37/ensembltogenesymbol/latest/EnsemblToGeneSymbol.tsv",
    "GRCh37/ensembltogenesymbol/latest/EnsemblToGeneSymbol.release_info",
    # SVs -- dbVar
    "GRCh37/dbVar/latest/DbVarSv.tsv",
    "GRCh37/dbVar/latest/DbVarSv.release_info",
    "GRCh38/dbVar/latest/DbVarSv.tsv",
    "GRCh38/dbVar/latest/DbVarSv.release_info",
    # SVs -- DGV
    "GRCh37/DGV/2016/DgvGoldStandardSvs.tsv",
    "GRCh37/DGV/2016/DgvGoldStandardSvs.release_info",
    "GRCh37/DGV/2016/DgvSvs.tsv",
    "GRCh37/DGV/2016/DgvSvs.release_info",
    "GRCh38/DGV/2016/DgvSvs.tsv",
    "GRCh38/DGV/2016/DgvSvs.release_info",
    # SVs -- ExAC
    "GRCh37/ExAC/r1/ExacCnv.tsv",
    "GRCh37/ExAC/r1/ExacCnv.release_info",
    # SVs -- 1000G
    "GRCh37/thousand_genomes/phase3/ThousandGenomesSv.tsv",
    "GRCh37/thousand_genomes/phase3/ThousandGenomesSv.release_info",
    # SVs -- gnomAD sv
    "GRCh37/gnomAD_SV/v2/GnomAdSv.tsv",
    "GRCh37/gnomAD_SV/v2/GnomAdSv.release_info",
    # MT -- mtDB
    "GRCh37/mtDB/latest/MtDb.tsv",
    "GRCh37/mtDB/latest/MtDb.release_info",
    # MT -- HelixMTdb
    "GRCh37/HelixMTdb/20190926/HelixMtDb.tsv",
    "GRCh37/HelixMTdb/20190926/HelixMtDb.release_info",
    # MT -- MITOMAP
    "GRCh37/MITOMAP/20200116/Mitomap.tsv",
    "GRCh37/MITOMAP/20200116/Mitomap.release_info",
    # Extra annotations.
    "GRCh37/extra-annos/20200704/ExtraAnno.tsv",
    "GRCh37/extra-annos/20200704/ExtraAnno.release_info",
    "GRCh37/extra-annos/20200704/ExtraAnnoField.tsv",
    "GRCh37/extra-annos/20200704/ExtraAnnoField.release_info",
]
VARFISH_FILES_WITH_IMPORT = VARFISH_FILES + ["import_versions.tsv"]


rule all:
    input:
        expand("{linkout}/{file}", linkout=VARFISH_SERVER_BACKGROUND_PATH, file=VARFISH_FILES_WITH_IMPORT),
        expand("{linkout}/{file}", linkout=VARFISH_ANNOTATOR_DB_PATH, file=VARFISH_ANNOTATOR_H2_FILES),
        expand("{linkout}/{file}", linkout=VARFISH_ANNOTATOR_PATH, file=VARFISH_ANNOTATOR_FILES_OUT)


rule generate_import_versions:
    input:
        VARFISH_FILES
    output:
        "import_versions.tsv"
    shell:
        r"""
        (
            echo -e "build\ttable_group\tversion"
            for i in $(find GRCh37 GRCh38 -maxdepth 2 -mindepth 2 -type d -not -path '*/reference/*')
            do
                echo $i | sed 's/\//\t/g'
            done | sort
        ) > {output}
        """


rule data_freeze_server_background_db:
    input:
        VARFISH_FILES_WITH_IMPORT
    output:
        expand("{linkout}/{file}", linkout=VARFISH_SERVER_BACKGROUND_PATH, file=VARFISH_FILES_WITH_IMPORT)
    params:
        linkout=VARFISH_SERVER_BACKGROUND_PATH
    shell:
        r"""
        for file in {input}
        do
            ln -sfr $file {params.linkout}/$file
        done
        """


rule data_freeze_varfish_annotator_h2_db:
    input:
        VARFISH_ANNOTATOR_H2_FILES,
    output:
        expand("{linkout}/{file}", linkout=VARFISH_ANNOTATOR_DB_PATH, file=VARFISH_ANNOTATOR_H2_FILES)
    params:
        linkout=VARFISH_ANNOTATOR_DB_PATH
    shell:
        r"""
        for file in {input}
        do
            ln -sfr $file {params.linkout}/$file
        done
        """


rule data_freeze_build_varfish_annotator_h2_db:
    input:
        expand("{linkout}/{file}", linkout=VARFISH_ANNOTATOR_DB_PATH, file=VARFISH_ANNOTATOR_H2_FILES)
    output:
        VARFISH_ANNOTATOR_DB_FILE
    params:
        release_path=RELEASE_PATH,
        data_release=DATA_RELEASE,
        clinvar_release=CLINVAR_RELEASE
    shell:
        r"""
        ANNOTATOR_VERSION=$(conda list varfish-annotator-cli | tail -1 | awk '{{print $2}}')
        varfish-annotator init-db \
            --db-release-info "varfish-annotator:v$ANNOTATOR_VERSION" \
            --db-release-info "varfish-annotator-db:r{params.data_release}" \
            \
            --ref-path {params.release_path}/varfish-annotator-db-{params.data_release}/GRCh37/reference/hs37d5/hs37d5.fa \
            \
            --db-release-info "clinvar:{params.clinvar_release}" \
            --clinvar-path {params.release_path}/varfish-annotator-db-{params.data_release}/GRCh37/clinvar/{params.clinvar_release}/clinvar_tsv_main/output/clinvar.b37.tsv.gz \
            \
            --db-path ./{params.release_path}/varfish-annotator-db-{params.data_release} \
            \
            --db-release-info "exac:r1.0" \
            --exac-path {params.release_path}/varfish-annotator-db-{params.data_release}/GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz \
            \
            --db-release-info "gnomad_exomes:r2.1" \
            $(for path in {params.release_path}/varfish-annotator-db-{params.data_release}/GRCh37/gnomAD_exomes/r2.1/download/gnomad.exomes.r2.1.sites.chr*.normalized.vcf.bgz; do \
                echo --gnomad-exomes-path $path; \
            done) \
            \
            --db-release-info "gnomad_genomes:r2.1" \
            $(for path in {params.release_path}/varfish-annotator-db-{params.data_release}/GRCh37/gnomAD_genomes/r2.1/download/gnomad.genomes.r2.1.sites.chr*.normalized.vcf.bgz; do \
                echo --gnomad-genomes-path $path; \
            done) \
            \
            --db-release-info "thousand_genomes:v3.20101123" \
            --thousand-genomes-path {params.release_path}/varfish-annotator-db-{params.data_release}/GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz \
            \
            --db-release-info "hgmd_public:ensembl_r75" \
            --hgmd-public {params.release_path}/varfish-annotator-db-{params.data_release}/GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.tsv
        """


rule download_jannovar_dbs:
    output:
        JANNOVAR_FILES,
    shell:
        r"""
        output=$(dirname {output[0]})
        jannovar download \
            -d hg19/refseq_curated \
            --download-dir $output
        jannovar download \
            -d hg19/ensembl \
            --download-dir $output
        """


rule data_freeze_varfish_annotator:
    input:
        VARFISH_ANNOTATOR_FILES,
    output:
        expand("{linkout}/{file}", linkout=VARFISH_ANNOTATOR_PATH, file=VARFISH_ANNOTATOR_FILES_OUT)
    params:
        linkout=VARFISH_ANNOTATOR_PATH
    shell:
        r"""
        for file in {input}
        do
            ln -sfr $file {params.linkout}/$(basename $file)
        done
        """


include: "snakefiles/acmg.snake"
include: "snakefiles/clinvar.snake"
include: "snakefiles/dbsnp.snake"
include: "snakefiles/dbvar.snake"
include: "snakefiles/dgv.snake"
include: "snakefiles/ensembl_regulatory.snake"
include: "snakefiles/ensembltogenesymbol.snake"
include: "snakefiles/ensembltorefseq.snake"
include: "snakefiles/exac.snake"
include: "snakefiles/genes.snake"
include: "snakefiles/gnomad.snake"
include: "snakefiles/hgnc.snake"
include: "snakefiles/hgmd.snake"
include: "snakefiles/helixmtdb.snake"
include: "snakefiles/hpo.snake"
include: "snakefiles/kegg.snake"
include: "snakefiles/knowngeneaa.snake"
include: "snakefiles/mgi.snake"
include: "snakefiles/mim2gene.snake"
include: "snakefiles/mitomap.snake"
include: "snakefiles/mtdb.snake"
include: "snakefiles/ncbi_gene.snake"
include: "snakefiles/reference.snake"
include: "snakefiles/refseqtoensembl.snake"
include: "snakefiles/refseqtogenesymbol.snake"
include: "snakefiles/tads.snake"
include: "snakefiles/thousand_genomes.snake"
include: "snakefiles/vista.snake"
include: "snakefiles/extra_annos.snake"
