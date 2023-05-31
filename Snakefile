from tools.sv_db_to_tsv import to_tsv


rule default:
    input:
        "annos/grch37/cadd/.done",
        "annos/grch37/dbnsfp-4.4a/.done",
        "annos/grch37/dbnsfp-4.4c/.done",
        "annos/grch37/dbscsnv/.done",
        "annos/grch37/helixmtdb/helixmtdb.vcf",
        "annos/grch37/gnomad_mtdna/gnomad_mtdna.vcf.gz",
        "annos/grch37/ucsc_conservation/ucsc_conservation.tsv",
        "annos/grch37/dbsnp/dbsnp.vcf.gz",
        "annos/grch37/gnomad_exomes/.done",
        "annos/grch37/gnomad_genomes/.done",
        "annos/grch38/cadd/.done",
        "annos/grch38/dbnsfp-4.4a/.done",
        "annos/grch38/dbnsfp-4.4c/.done",
        "annos/grch38/gnomad_exomes/.done",
        "annos/grch38/gnomad_genomes/.done",
        "annos/grch38/gnomad_mtdna/gnomad_mtdna.vcf.gz",
        "annos/grch38/helixmtdb/helixmtdb.vcf",
        "features/grch37/tads/imr90.bed",
        "features/grch37/tads/hesc.bed",
        "features/grch37/gene_regions/refseq.bed.gz",
        "features/grch37/gene_regions/ensembl.bed.gz",
        "features/grch37/masked/repeat.bed.gz",
        "features/grch37/masked/segdup.bed.gz",
        "genes/hgnc/hgnc_info.jsonl",
        "genes/ncbi/gene_info.jsonl",
        "genes/dbnsfp/genes.tsv.gz",
        "genes/xlink/ensembl.tsv",
        "genes/xlink/hgnc.tsv",
        "genes/mim2gene/mim2gene.tsv",
        "tracks/grch37/ucsc_genomicSuperDups.bed.gz",
        "tracks/grch37/ucsc_rmsk.bed.gz",
        "tracks/grch37/ucsc_fixSeqLiftOverPsl.bed.gz",
        "tracks/grch37/ucsc_altSeqLiftOverPsl.bed.gz",
        "vardbs/grch37/strucvar/clinvar.bed.gz",
        "vardbs/grch37/strucvar/dbvar.bed.gz",
        "vardbs/grch37/strucvar/dgv.bed.gz",
        "vardbs/grch37/strucvar/dgv_gs.bed.gz",
        "vardbs/grch37/strucvar/g1k.bed.gz",
        "vardbs/grch37/strucvar/gnomad_sv.bed.gz",
        "vardbs/grch37/strucvar/exac.bed.gz",


include: "snakefiles/annos.smk"
include: "snakefiles/genes.smk"
include: "snakefiles/features.smk"
include: "snakefiles/vardbs-grch37-strucvars.smk"
include: "snakefiles/tracks-grch37.smk"
include: "snakefiles/reference.smk"
