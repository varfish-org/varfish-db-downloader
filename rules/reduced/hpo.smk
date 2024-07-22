## Rules to create build viguno phenotypes database subsets (dev/exomes).
#
# We will copy the full HPO (text and binary) but reduce the simulation count.


rule subset_viguno_pheno_exomes:  # -- create exomes subset
    input:
        obo="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/full/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
    output:
        obo="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
        manifest="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/MANIFEST.txt",
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_viguno=RE_VERSION,
    shell:
        r"""
        cp -a {input.obo} {output.obo}
        cp -a {input.hpoa} {output.hpoa}
        cp -a {input.phenotype_to_genes} {output.phenotype_to_genes}
        cp -a {input.bin} {output.bin}
        cp -a {input.spec_yaml} {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """


rule subset_worker_pheno_dev:  # -- create development subset
    input:
        obo="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/full/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
    output:
        obo="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
        manifest="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/MANIFEST.txt",
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_viguno=RE_VERSION,
    threads: VIGUNO_SIMULATE_THREADS
    resources:
        mem_mb_per_cpu="2GB",
    shell:
        r"""
        cp -a {input.obo} {output.obo}
        cp -a {input.hpoa} {output.hpoa}
        cp -a {input.phenotype_to_genes} {output.phenotype_to_genes}
        cp -a {input.bin} {output.bin}
        cp -a {input.spec_yaml} {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """
