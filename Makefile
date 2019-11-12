DATA_RELEASE = 20191107
ANNOTATOR_VERSION = 0.9


all: pack_server pack_annotator pack_jannovar


pack_server:
	tar chzvf \
		varfish-server-background-db-$(DATA_RELEASE).tar.gz \
		varfish-server-background-db-$(DATA_RELEASE)/
	sha256sum \
		varfish-server-background-db-$(DATA_RELEASE).tar.gz \
		> varfish-server-background-db-$(DATA_RELEASE).tar.gz.sha256


varfish-annotator-db-$(DATA_RELEASE).db.h2:
	varfish-annotator init-db \
		--db-release-info "varfish-annotator:v$(ANNOTATOR_VERSION)" \
		--db-release-info "varfish-annotator-db:r$(DATA_RELEASE)" \
		\
		--ref-path varfish-annotator-db-$(DATA_RELEASE)/GRCh37/reference/hs37d5/hs37d5.fa \
		\
		--db-release-info "clinvar:2019-06-22" \
		--clinvar-path varfish-annotator-db-$(DATA_RELEASE)/GRCh37/clinvar/latest/clinvar_tsv_main/output/clinvar_allele_trait_pairs.single.b37.tsv.gz \
		--clinvar-path varfish-annotator-db-$(DATA_RELEASE)/GRCh37/clinvar/latest/clinvar_tsv_main/output/clinvar_allele_trait_pairs.multi.b37.tsv.gz \
		\
		--db-path ./varfish-annotator-db-$(DATA_RELEASE) \
		\
		--db-release-info "exac:r1.0" \
		--exac-path varfish-annotator-db-$(DATA_RELEASE)/GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz \
		\
		--db-release-info "gnomad_exomes:r2.1" \
		$(shell for path in varfish-annotator-db-$(DATA_RELEASE)/GRCh37/gnomAD_exomes/r2.1/download/gnomad.exomes.r2.1.sites.chr*.normalized.vcf.bgz; do echo --gnomad-exomes-path $$path; done) \
		\
		--db-release-info "gnomad_genomes:r2.1" \
		$(shell for path in varfish-annotator-db-$(DATA_RELEASE)/GRCh37/gnomAD_genomes/r2.1/download/gnomad.genomes.r2.1.sites.chr*.normalized.vcf.bgz; do echo --gnomad-genomes-path $$path; done) \
		\
		--db-release-info "thousand_genomes:v3.20101123" \
		--thousand-genomes-path varfish-annotator-db-$(DATA_RELEASE)/GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz \
		\
		--db-release-info "hgmd_public:ensembl_r75" \
		--hgmd-public varfish-annotator-db-$(DATA_RELEASE)/GRCh37/hgmd_public/ensembl_r75/HgmdPublicLocus.tsv


pack_annotator: varfish-annotator-db-$(DATA_RELEASE).db.h2
	gzip -c \
		varfish-annotator-db-$(DATA_RELEASE).db.h2 \
		> varfish-annotator-db-$(DATA_RELEASE).db.h2.gz
	sha256sum \
		varfish-annotator-db-$(DATA_RELEASE).h2.db.gz \
		> varfish-annotator-db-$(DATA_RELEASE).h2.db.gz.sha256


download_jannovar:
	jannovar download \
		-d hg19/refseq_curated \
		--download-dir varfish-annotator-transcripts-$(DATA_RELEASE)
	jannovar download \
		-d hg19/ensembl \
		--download-dir varfish-annotator-transcripts-$(DATA_RELEASE)


pack_jannovar: download_jannovar
	tar czvf \
		varfish-annotator-transcripts-$(DATA_RELEASE).tar.gz \
		varfish-annotator-transcripts-$(DATA_RELEASE)/*.ser
	sha256sum \
		varfish-annotator-transcripts-$(DATA_RELEASE).tar.gz \
		> varfish-annotator-transcripts-$(DATA_RELEASE).tar.gz.sha256

