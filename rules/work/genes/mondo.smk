## Rules related to MONDO data


rule genes_mondo_download:  # -- download orphadatas
    output:
        obo="work/genes/mondo/{version}/mondo.obo",
        omim_unmapped="work/genes/mondo/{version}/omim_unmapped_terms.tsv",
    shell:
        """
        wget -O {output.obo} \
            http://purl.obolibrary.org/obo/mondo.obo
        wget -O {output.omim_unmapped} \
            https://raw.githubusercontent.com/monarch-initiative/mondo-ingest/main/src/ontology/reports/omim_unmapped_terms.tsv
        """
