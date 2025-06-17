## Rules to create build annonars regions annotation database..


rule work_annonars_regions_download:  # -- download clingen regions
    output:
        "work/download/clingen/{genome_release}/{today}/ClinGen_region_curation_list_{genome_release}.tsv",
    shell:
        r"""
        if [[ "{wildcards.genome_release}" == "grch38" ]]; then
            GENOME=GRCh37
        else
            GENOME=GRCh38
        fi

        wget -O {output} \
            https://ftp.clinicalgenome.org/ClinGen_region_curation_list_$GENOME.tsv
        """


rule output_annonars_regions:  # -- build annonars regions RocksDB file
    input:
        "work/download/clingen/{genome_release}/{date}/ClinGen_region_curation_list_{genome_release}.tsv",
    output:
        rocksdb_identity=(
            "output/full/annonars/regions-{genome_release}-{date}+{v_annonars}/" "rocksdb/IDENTITY"
        ),
        spec_yaml=("output/full/annonars/regions-{genome_release}-{date}+{v_annonars}/spec.yaml"),
        manifest=("output/full/annonars/regions-{genome_release}-{date}+{v_annonars}/MANIFEST.txt"),
    wildcard_constraints:
        v_refseq=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        annonars regions import -vvv \
            --genome-release {wildcards.genome_release} \
            --path-in-clingen {input} \
            --path-out-rocksdb $(dirname {output.rocksdb_identity})

        varfish-db-downloader tpl \
            --template rules/output/annonars/regions.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.date}+{wildcards.v_annonars} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """
