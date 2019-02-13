from snakemake import shell
shell.prefix("set -x; set -euo pipefail; ")


CHROMS = (list(range(1, 23)) + ['X'])


rule all:
    input:
        # reference files
        "GRCh37/reference/hs37d5/hs37d5.fa",
        "GRCh38/reference/hs38/hs38.fa",
        # clinvar
        "GRCh37/clinvar/latest/Clinvar.tsv",
        "GRCh37/clinvar/latest/Clinvar.release_info",
        "GRCh38/clinvar/latest/Clinvar.tsv",
        "GRCh38/clinvar/latest/Clinvar.release_info",
        # dbsnp
        "GRCh37/dbSNP/b151/Dbsnp.tsv",
        "GRCh37/dbSNP/b151/Dbsnp.release_info",
        # exac
        "GRCh37/ExAC/r1/Exac.tsv",
        "GRCh37/ExAC/r1/Exac.release_info",
        # gnomad
        expand("GRCh37/gnomAD_exomes/r2.1/GnomadExomes.{chrom}.tsv", chrom=CHROMS),
        expand("GRCh37/gnomAD_exomes/r2.1/GnomadExomes.{chrom}.release_info", chrom=CHROMS),
        expand("GRCh37/gnomAD_genomes/r2.1/GnomadGenomes.{chrom}.tsv", chrom=CHROMS),
        expand("GRCh37/gnomAD_genomes/r2.1/GnomadGenomes.{chrom}.release_info", chrom=CHROMS),
        # hgmd
        "GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.tsv",
        "GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.release_info",
        # hgnc
        "GRCh37/hgnc/latest/Hgnc.tsv",
        "GRCh37/hgnc/latest/Hgnc.release_info",
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
