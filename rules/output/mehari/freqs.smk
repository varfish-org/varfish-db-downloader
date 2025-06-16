## Rules related to mehari frequencies.

import os


rule output_mehari_freqs_build:  # -- build frequency tables for mehari
    input:
        gnomad_genomes="work/download/annos/{genome_release}/seqvars/gnomad_genomes/{v_gnomad_genomes}/.done",
        gnomad_exomes="work/download/annos/{genome_release}/seqvars/gnomad_exomes/{v_gnomad_exomes}/.done",
        gnomad_mtdna="work/annos/{genome_release}/seqvars/gnomad_mtdna/{v_gnomad_mtdna}/gnomad_mtdna.vcf.gz",
        helixmtdb="work/annos/{genome_release}/seqvars/helixmtdb/{v_helixmtdb}/helixmtdb.vcf.gz",
    output:
        rocksdb_identity=(
            "output/full/mehari/freqs-{genome_release}-{v_gnomad_genomes}+{v_gnomad_exomes}+"
            "{v_gnomad_mtdna}+{v_helixmtdb}+{v_annonars}/rocksdb/IDENTITY"
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
        output_rocksdb=$(dirname {output.rocksdb_identity})

        build-args()
        {{
            path=$1
            arg=$2
            regex=$3

            for file in $(find $path -type f -and -name "*.vcf.bgz" \
                    | sort \
                    | egrep -i "$regex"); do
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
            --path-out-rocksdb "$output_rocksdb" \
            \
            --path-gnomad-mtdna {input.gnomad_mtdna} \
            --path-helixmtdb {input.helixmtdb} \
            \
            $(build-args $(dirname {input.gnomad_genomes}) --path-gnomad-genomes-auto "sites\.(chr)?[0-9]+\.") \
            $(build-args $(dirname {input.gnomad_genomes}) --path-gnomad-genomes-xy   "sites\.(chr)?[XY]\.") \
            \
            $(build-args $(dirname {input.gnomad_exomes})  --path-gnomad-exomes-auto  "sites\.(chr)?[0-9]+\.") \
            $(build-args $(dirname {input.gnomad_exomes})  --path-gnomad-exomes-xy    "sites\.(chr)?[XY]\.")

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
