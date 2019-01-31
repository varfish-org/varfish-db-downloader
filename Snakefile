from snakemake import shell
shell.prefix("set -x; set -euo pipefail; ")

rule all:
    input:
        "GRCh37/reference/hs37d5/hs37d5.fa",
        "GRCh37/clinvar/latest/clinvar.tsv",
        "GRCh38/clinvar/latest/clinvar.tsv",
        "GRCh37/dbSNP/b151/All_20180423.vcf.gz",
        "GRCh37/gnomAD_exomes/r2.1/gnomad.exomes.r2.1.sites.vcf.bgz",
        "GRCh37/gnomAD_genomes/r2.1/gnomad.genomes.r2.1.sites.vcf.bgz",
        "GRCh37/hgmd_public/hgmd_public.tsv",
        "GRCh37/hgnc/hgnc.tsv",
        "GRCh37/knowngeneaa/knownGeneAA.tsv",
        "GRCh37/ncbi_gene/ncbi_gene_info.tsv",
        "GRCh37/ncbi_gene/ncbi_gene_rif.tsv",
        "GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz",
        "GRCh37/kegg/kegg.release_info",
        "GRCh37/mim2gene/mim2gene.release_info",
        "GRCh37/hpo/hpo.release_info",
        "GRCh38/reference/hs38/hs38.fa",

include: "snakefiles/clinvar.snake"
include: "snakefiles/dbsnp.snake"
include: "snakefiles/exac.snake"
include: "snakefiles/gnomad.snake"
include: "snakefiles/hgnc.snake"
include: "snakefiles/hgmd.snake"
include: "snakefiles/hpo.snake"
include: "snakefiles/kegg.snake"
include: "snakefiles/knowngeneaa.snake"
include: "snakefiles/ncbi_gene.snake"
include: "snakefiles/mim2gene.snake"
include: "snakefiles/reference.snake"
include: "snakefiles/thousand_genomes.snake"
