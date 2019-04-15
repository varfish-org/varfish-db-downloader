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
        "GRCh37/vista/latest/Vista.tsv",
        "GRCh37/vista/latest/Vista.release_info.tsv",
        # ENSEMBL Regulatory Features
        "GRCh37/ensembl/latest/EnsemblRegulatory.tsv",
        "GRCh37/ensembl/latest/EnsemblRegulatory.release_info.tsv",
        # Dixon 2012 TAD BED files
        "GRCh37/tads/dixon2012/hESC_domains_hg19.bed",
        "GRCh37/tads/dixon2012/IMR90_domains_hg19.bed",
        "GRCh37/tads/dixon2012/Dixon2012Tads.release_info.bed",
        ########################################
        ### SVs files, need post processing step
        ########################################
        # SVs -- dbVar
        expand("{genomebuild}/dbVar/latest/download/{type}/{genomebuild}.nr_{type}.{ending}{md5}",
            genomebuild=['GRCh37', 'GRCh38'],
            type=['deletions', 'insertions', 'duplications'],
            ending=['bed.gz', 'bedpe.gz', 'tsv.gz', 'acmg_genes.tsv.gz'],
            md5=['', '.md5']
        ),
        # SVs -- DGV
        "GRCh37/DGV/2016/download/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
        "GRCh37/DGV/2016/download/GRCh37_hg19_variants_2016-05-15.txt",
        "GRCh37/DGV/2016/download/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3.md5",
        "GRCh37/DGV/2016/download/GRCh37_hg19_variants_2016-05-15.txt.md5",
        "GRCh38/DGV/2016/download/GRCh38_hg38_variants_2016-08-31.txt",
        "GRCh38/DGV/2016/download/GRCh38_hg38_variants_2016-08-31.txt.md5",
        # SVs -- ExAC
        "GRCh37/ExAC/r1/download/exac-final-cnv.gene.scores071316",
        "GRCh37/ExAC/r1/download/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
        "GRCh37/ExAC/r1/download/README.cnv_gene_scores",
        "GRCh37/ExAC/r1/download/README.cnv_bed",
        "GRCh37/ExAC/r1/download/md5sum.txt",
        # SVs -- 1000G
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz",
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz.tbi",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi",
        "GRCh37/thousand_genomes/phase3/download/integrated_call_samples_v3.20130502.ALL.panel",
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz.md5",
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz.tbi.md5",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.md5",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi.md5",
        "GRCh37/thousand_genomes/phase3/download/integrated_call_samples_v3.20130502.ALL.panel.md5",


include: "snakefiles/clinvar.snake"
include: "snakefiles/dbsnp.snake"
include: "snakefiles/dbvar.snake"
include: "snakefiles/dgv.snake"
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
include: "snakefiles/acmg.snake"
include: "snakefiles/vista.snake"
include: "snakefiles/ensembl_regulatory.snake"
include: "snakefiles/tads.snake"
