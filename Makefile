# scripts/ -- contains scripts for conversion. each $(NAME).sh gets a databases/$(NAME)/download folder and download rule
# helpers/scripts -- contains scripts for files needed for conversion. each $(NAME).sh gets a helpers/data/$(NAME)/download folder and download rule

# exclude rules from download.
EXTERNAL = annotation smallvariant coding_regions
DBS = $(patsubst scripts/%.sh,%,$(wildcard scripts/*.sh))
# the order of helpers is important
HELPERS = reference ensembl refseq coding_regions
DLS = $(patsubst %,dl_%,$(filter-out $(EXTERNAL),$(HELPERS) $(DBS)))
PATHS = $(patsubst %,databases/%/download,$(DBS)) $(patsubst %,helpers/data/%/download,$(HELPERS))
.PHONY: convert download init helpers help $(DBS) $(DLS) $(HELPERS) $(PATHS) $(EXTERNAL)

print:
	@echo $(DBS)
	@echo $(DLS)
	@echo $(PATHS)

convert: helpers $(DBS)
$(DBS):
	$(MAKE) -C scripts $@

helpers: download $(HELPERS)
$(HELPERS):
	$(MAKE) -C helpers/scripts $@

download: $(DLS)
$(DLS): init
	$(MAKE) -C downloads $(patsubst dl_%,%,$@)

init: $(PATHS)
$(PATHS):
	@mkdir -p $@

help:
	@echo init - create folder structure
	@echo download - download all \(CAUTION: takes time\), including reference
	@echo helpers - prepare helpers \(required for convert, downloads must be finished\)
	@echo convert - convert all \(downloads must be finished or make files in some other way available\)
	@echo 
	@echo Address single rules:
	@echo - $(DLS)
	@echo - $(HELPERS)
	@echo - $(DBS)
