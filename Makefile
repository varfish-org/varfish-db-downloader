DATA_RELEASE = 20200922


all: pack_server pack_annotator pack_jannovar


pack_server:
	tar \
		--owner=0 \
		--group=0 \
		-chzvf varfish-server-background-db-$(DATA_RELEASE).tar.gz \
		varfish-server-background-db-$(DATA_RELEASE)/
	sha256sum \
		varfish-server-background-db-$(DATA_RELEASE).tar.gz \
		> varfish-server-background-db-$(DATA_RELEASE).tar.gz.sha256


pack_annotator:
	tar \
		--owner=0 \
		--group=0 \
		-chzvf varfish-annotator-db-$(DATA_RELEASE).tar.gz \
		varfish-annotator-db-$(DATA_RELEASE)/
	sha256sum \
		varfish-annotator-db-$(DATA_RELEASE).tar.gz \
		> varfish-annotator-db-$(DATA_RELEASE).tar.gz.sha256


pack_jannovar:
	tar \
		--owner=0 \
		--group=0 \
		-czvf jannovar-db-$(DATA_RELEASE).tar.gz \
		jannovar-db-$(DATA_RELEASE)/*.ser
	sha256sum \
		jannovar-db-$(DATA_RELEASE).tar.gz \
		> jannovar-db-$(DATA_RELEASE).tar.gz.sha256

