## Rules related to mehari frequencies.

import os


def input_gnomad_genomes_auto(wildcards):
    return expand(
        "work/download/annos/{g}/seqvars/gnomad_genomes/{v}/gnomad.genomes.{t}{v}.sites.chr{c}.vcf.bgz",
        g=wildcards.genome_release,
        v=gnomad_versions[wildcards.genome_release],
        t="r" if wildcards.genome_release == "grch37" else "v",
        c=CHROMS_AUTO,
    )


def input_gnomad_genomes_gono(wildcards):
    chroms = CHROMS_AUTO + ("X",)
    if wildcards.genome_release == "grch38":
        chroms += ("Y",)
    return expand(
        "work/download/annos/{g}/seqvars/gnomad_genomes/{v}/gnomad.genomes.{t}{v}.sites.chr{c}.vcf.bgz",
        g=wildcards.genome_release,
        v=gnomad_versions[wildcards.genome_release],
        t="r" if wildcards.genome_release == "grch37" else "v",
        c=chroms,
    )


def input_gnomad_genomes_tbi(wildcards):
    chroms = CHROMS_AUTO + ("X",)
    if wildcards.genome_release == "grch38":
        chroms += ("Y",)
    return expand(
        "work/download/annos/{g}/seqvars/gnomad_genomes/{v}/gnomad.genomes.{t}{v}.sites.chr{c}.vcf.bgz.tbi",
        g=wildcards.genome_release,
        v=gnomad_versions[wildcards.genome_release],
        t="r" if wildcards.genome_release == "grch37" else "v",
        c=chroms,
    )


def input_gnomad_exomes_auto(wildcards):
    return expand(
        "work/download/annos/{g}/seqvars/gnomad_exomes/{v}/gnomad.exomes.{t}{v}.sites.chr{c}.vcf.bgz",
        g=wildcards.genome_release,
        v=gnomad_versions[wildcards.genome_release],
        t="r" if wildcards.genome_release == "grch37" else "v",
        c=CHROMS_AUTO,
    )


def input_gnomad_exomes_gono(wildcards):
    return expand(
        "work/download/annos/{g}/seqvars/gnomad_exomes/{v}/gnomad.exomes.{t}{v}.sites.chr{c}.vcf.bgz",
        g=wildcards.genome_release,
        v=gnomad_versions[wildcards.genome_release],
        t="r" if wildcards.genome_release == "grch37" else "v",
        c=["X", "Y"],
    )


def input_gnomad_exomes_tbi(wildcards):
    return expand(
        "work/download/annos/{g}/seqvars/gnomad_exomes/{v}/gnomad.exomes.{t}{v}.sites.chr{c}.vcf.bgz.tbi",
        g=wildcards.genome_release,
        v=gnomad_versions[wildcards.genome_release],
        t="r" if wildcards.genome_release == "grch37" else "v",
        c=CHROMS,
    )


rule output_mehari_freqs_build:  # -- build frequency tables for mehari
    input:
        "work/download/annos/{genome_release}/seqvars/gnomad_genomes/{v_gnomad_genomes}/.done",
        "work/download/annos/{genome_release}/seqvars/gnomad_exomes/{v_gnomad_exomes}/.done",
        gnomad_genomes_auto=input_gnomad_genomes_auto,
        gnomad_genomes_gono=input_gnomad_genomes_gono,
        gnomad_genomes_tbi=input_gnomad_genomes_tbi,
        gnomad_exomes_auto=input_gnomad_exomes_auto,
        gnomad_exomes_gono=input_gnomad_exomes_gono,
        gnomad_exomes_tbi=input_gnomad_exomes_tbi,
        gnomad_mtdna="work/annos/{genome_release}/seqvars/gnomad_mtdna/{v_gnomad_mtdna}/gnomad_mtdna.vcf.gz",
        helixmtdb="work/annos/{genome_release}/seqvars/helixmtdb/{v_helixmtdb}/helixmtdb.vcf.gz",
        validate_script="scripts/validate_rocksdb.sh",
    output:
        rocksdb_dir=directory(
            "output/full/mehari/freqs-{genome_release}-{v_gnomad_genomes}+{v_gnomad_exomes}+"
            "{v_gnomad_mtdna}+{v_helixmtdb}+{v_annonars}/rocksdb"
        ),
        spec_yaml=(
            "output/full/mehari/freqs-{genome_release}-{v_gnomad_genomes}+{v_gnomad_exomes}+"
            "{v_gnomad_mtdna}+{v_helixmtdb}+{v_annonars}/spec.yaml"
        ),
        manifest=(
            "output/full/mehari/freqs-{genome_release}-{v_gnomad_genomes}+{v_gnomad_exomes}+"
            "{v_gnomad_mtdna}+{v_helixmtdb}+{v_annonars}/MANIFEST.txt"
        ),
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_gnomad_genomes=RE_VERSION,
        v_gnomad_exomes=RE_VERSION,
        v_gnomad_mtdna=RE_VERSION,
        v_helixmtdb=r"\d+",
        v_annonars=RE_VERSION,
    shell:
        r"""
        if [[ "${{CI:-false}}" == "true" ]]; then
            echo "Skipping rule output_mehari_freqs_build because CI=true"
            mkdir -p {output.rocksdb_dir}
            touch {output.spec_yaml} {output.manifest}
            exit 0
        fi
        
        build-args()
        {{
            arg=$1
            files=$2

            for file in $files; do
                echo $arg $file
            done
        }}

        annonars freqs import \
            --genome-release "{wildcards.genome_release}" \
            \
            --gnomad-genomes-version "{wildcards.v_gnomad_genomes}" \
            --gnomad-exomes-version "{wildcards.v_gnomad_exomes}" \
            --gnomad-mtdna-version "{wildcards.v_gnomad_mtdna}" \
            --helixmtdb-version "{wildcards.v_helixmtdb}" \
            \
            --path-out-rocksdb {output.rocksdb_dir} \
            \
            --path-gnomad-mtdna {input.gnomad_mtdna} \
            --path-helixmtdb {input.helixmtdb} \
            \
            build-args "--path-gnomad-genomes-auto" "{input.gnomad_genomes_auto}" \
            build-args "--path-gnomad-genomes-xy" "{input.gnomad_genomes_gono}" \
            \
            build-args "--path-gnomad-exomes-auto" "{input.gnomad_exomes_auto}" \
            build-args "--path-gnomad-exomes-xy" "{input.gnomad_exomes_gono}"

        bash {input.validate_script} "{output.rocksdb_dir}"

        varfish-db-downloader tpl \
            --template rules/output/mehari/freqs.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_gnomad_genomes}+{wildcards.v_gnomad_exomes}+{wildcards.v_gnomad_mtdna}+{wildcards.v_helixmtdb}+{wildcards.v_annonars} \
            --value v_gnomad_genomes={wildcards.v_gnomad_genomes} \
            --value v_gnomad_exomes={wildcards.v_gnomad_exomes} \
            --value v_gnomad_mtdna={wildcards.v_gnomad_mtdna} \
            --value v_helixmtdb={wildcards.v_helixmtdb} \
            \
            --value v_annovars={wildcards.v_annonars} \
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
