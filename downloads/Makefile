URL_CLINVAR_SINGLE = https://raw.githubusercontent.com/macarthur-lab/clinvar/master/output/b37/single/clinvar_allele_trait_pairs.single.b37.tsv.gz
URL_CLINVAR_MULTI = https://raw.githubusercontent.com/macarthur-lab/clinvar/master/output/b37/multi/clinvar_allele_trait_pairs.multi.b37.tsv.gz
URL_OMIM = ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen
URL_HGNC = ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
URL_HPO = http://compbio.charite.de/jenkins/job/hpo.annotations.2018/lastSuccessfulBuild/artifact/misc_2018/phenotype.hpoa
URL_EXAC = ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
URL_GNOMAD = https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz
URL_DBSNP = ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
URL_REF = ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

CLINVAR_SINGLE = ../databases/clinvar/download/$(notdir $(URL_CLINVAR_SINGLE))
CLINVAR_MULTI = ../databases/clinvar/download/$(notdir $(URL_CLINVAR_MULTI))
OMIM = ../databases/omim/download/$(notdir $(URL_OMIM))
HGNC = ../databases/hgnc/download/$(notdir $(URL_HGNC))
HPO = ../databases/hpo/download/$(notdir $(URL_HPO))
EXAC = ../databases/exac/download/$(notdir $(URL_EXAC))
GNOMAD = ../databases/gnomad/download/$(notdir $(URL_GNOMAD))
DBSNP = ../databases/dbsnp/download/$(notdir $(URL_DBSNP))
REF = $(patsubst %.gz,%,$(notdir $(URL_REF)))

clinvar: $(CLINVAR_SINGLE) $(CLINVAR_MULTI)
$(CLINVAR_SINGLE):
	wget -c $(URL_CLINVAR_SINGLE) -O $@
$(CLINVAR_MULTI):
	wget -c $(URL_CLINVAR_MULTI) -O $@

omim: $(OMIM)
$(OMIM):
	wget -c $(URL_OMIM) -O $@

hgnc: $(HGNC)
$(HGNC):
	wget -c $(URL_HGNC) -O $@

hpo: $(HPO)
$(HPO):
	wget -c $(URL_HPO) -O $@

exac: $(EXAC)
$(EXAC):
	wget -c $(URL_EXAC) -O $@
	wget -c $(URL_EXAC).tbi -O $@.tbi

gnomad: $(GNOMAD)
$(GNOMAD):
	wget -c $(URL_GNOMAD) -O $@
	wget -c $(URL_GNOMAD).tbi -O $@.tbi

dbsnp: $(DBSNP)
$(DBSNP):
	wget -c $(URL_DBSNP) -O $@
	wget -c $(URL_DBSNP).tbi -O $@.tbi

ref: $(REF)
$(REF):
	wget -c $(URL_REF) -O- \
	| gunzip -c > $@ || true

case:
	@echo -e
	@echo -e Place ped files \(yes, multiple are allowed\) in databases/case/download
	@echo -e

kegg:
	@echo -e
	@echo -e \#\# EnsemblToKegg
	@echo -e \* https://genome.ucsc.edu
	@echo -e \* Tools -\> Table Browser
	@echo -e \\t- group: All tables
	@echo -e \\t- assembly: GRCh37
	@echo -e \\t- table: keggPathway
	@echo -e \\t- output format: selected fields from primary and related tables
	@echo -e \\t- output file: genetokegg.tsv
	@echo -e \\t- get output
	@echo -e \\t\\t- Linked Tables -\> knownToEnsembl -\> allow selection from checked tables
	@echo -e \\t\\t- Linked Tables -\> ensGtp -\> allow selection from checked tables
	@echo -e \\t\\t- Select Fields from keggPathway -\> mapID
	@echo -e \\t\\t- ensGtp fields -\> gene
	@echo -e \\t\\t- get output
	@echo -e
	@echo -e \#\# RefseqToKegg
	@echo -e \* https://genome.ucsc.edu
	@echo -e \* Tools -\> Table Browser
	@echo -e \\t- group: All tables
	@echo -e \\t- assembly: GRCh37
	@echo -e \\t- table: keggPathway
	@echo -e \\t- output format: selected fields from primary and related tables
	@echo -e \\t- output file: genetokegg.tsv
	@echo -e \\t- get output
	@echo -e \\t\\t- Linked Tables -\> knownToKeggEntrez -\> allow selection from checked tables
	@echo -e \\t\\t- knownToKeggEntrez fields -\> keggEntrez
	@echo -e \\t\\t- get output
	@echo -e 
	@echo -e \#\# KeggInfo
	@echo -e \* https://genome.ucsc.edu
	@echo -e \* Tools -\> Table Browser
	@echo -e \\t- group: All tables
	@echo -e \\t- assembly: GRCh37
	@echo -e \\t- table: keggMapDesc
	@echo -e \\t- output format: selected fields from primary and related tables
	@echo -e \\t- output file: kegginfo.tsv
	@echo -e \\t- get output
	@echo -e \\t\\t- Select Fields from keggMapDesc -\> mapID
	@echo -e \\t\\t- Select Fields from keggMapDesc -\> description
	@echo -e \\t\\t- get output
	@echo -e
	@echo -e Save both files to databases/kegg/download
	@echo -e
