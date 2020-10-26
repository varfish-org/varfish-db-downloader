DATA_RELEASE = 20201006
RELEASE_PATH = releases
JANNOVAR_RELEASE = $(DATA_RELEASE)


all: pack_server pack_annotator pack_jannovar


pack_server:
	cd $(RELEASE_PATH) && tar \
		--owner=0 \
		--group=0 \
		-chzvf varfish-server-background-db-$(DATA_RELEASE).tar.gz \
		varfish-server-background-db-$(DATA_RELEASE)/
	cd $(RELEASE_PATH) && sha256sum \
		varfish-server-background-db-$(DATA_RELEASE).tar.gz \
		> varfish-server-background-db-$(DATA_RELEASE).tar.gz.sha256


pack_annotator:
	cd $(RELEASE_PATH) && tar \
		--owner=0 \
		--group=0 \
		-chzvf varfish-annotator-db-$(DATA_RELEASE).tar.gz \
		varfish-annotator-db-$(DATA_RELEASE)/
	cd $(RELEASE_PATH) && sha256sum \
		varfish-annotator-db-$(DATA_RELEASE).tar.gz \
		> varfish-annotator-db-$(DATA_RELEASE).tar.gz.sha256


pack_jannovar:
	cd $(RELEASE_PATH) && tar \
		--owner=0 \
		--group=0 \
		-czvf jannovar-db-$(JANNOVAR_RELEASE).tar.gz \
		jannovar-db-$(JANNOVAR_RELEASE)/*.ser
	cd $(RELEASE_PATH) && sha256sum \
		jannovar-db-$(JANNOVAR_RELEASE).tar.gz \
		> jannovar-db-$(JANNOVAR_RELEASE).tar.gz.sha256

