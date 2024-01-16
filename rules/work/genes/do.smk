## Rules related to DO download


rule work_disease_ontology_download:  # -- download DO files
    output:
        unmapped_csv="work/download/do/{date}/omim-unmapped.csv",
        omimindo_tsv="work/download/do/{date}/OMIMinDO.tsv",
        omim_import="work/download/do/{date}/omim_import.obo",
    shell:
        r"""
        wget -O {output.unmapped_csv} \
            https://github.com/DiseaseOntology/HumanDiseaseOntology/raw/main/src/deprecated/reports/omim-unmapped.csv
        wget -O {output.omimindo_tsv} \
            https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/DOreports/OMIMinDO.tsv
        wget -O {output.omim_import} \
            https://github.com/DiseaseOntology/HumanDiseaseOntology/raw/main/src/deprecated/DO_NON_Production_Files/omim_import.obo
        perl -p -i -e 's/(is_a:[^!]*?),.*/\1/g' {output.omim_import}
        perl -p -i -e 's/OMIM:PS/OMIM:/g' {output.omim_import}
        """
