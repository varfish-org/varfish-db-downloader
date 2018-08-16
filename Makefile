DBS = $(patsubst scripts/%.sh,%,$(wildcard scripts/*.sh))
DLS = $(patsubst %,dl_%,$(DBS))
PATHS = $(patsubst %,databases/%/download,$(DBS))
.PHONY: convert download dl_ref init help $(DBS) $(DLS)

convert: $(DBS)
$(DBS):
	$(MAKE) -C scripts $@

download: dl_ref $(DLS)
$(DLS): init
	$(MAKE) -C downloads $(patsubst dl_%,%,$@)

dl_ref:
	$(MAKE) -C downloads ref

init: $(PATHS)
$(PATHS):
	@mkdir -p $@

help:
	@echo init - create folder structure
	@echo convert - convert all \(downloads have to be finished or make files in some other way available\)
	@echo download - download all \(CAUTION: takes time\)
	@echo dl_ref - download reference
	@echo 
	@echo Address single rules:
	@echo - $(DLS)
	@echo - $(DBS)
	