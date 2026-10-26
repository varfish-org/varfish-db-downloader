rule generate_mehari_dbs_setup_script:
    output:
        script="output/full/mehari/setup_mehari_dbs.sh",
    params:
        template="rules/output/mehari/setup_mehari_dbs.sh.j2",
        primary_tx_source="refseq",  # Choice of: refseq, ensembl, ensembl-and-refseq
    run:
        from jinja2 import Template

        mehari_freq_versions = {
            genombuild: (
                f"freqs-{genombuild}-"
                f"{gnomad_versions[genombuild]}+"
                f"{gnomad_versions[genombuild]}+"
                f"{DV.gnomad_mtdna}+"
                f"{DV.helixmtdb}+"
                f"{PV.annonars}"
            )
            for genombuild in GENOMEBUILDS
        }

        mehari_tx_files = {
            genombuild: (
                f"mehari-data-txs-{genomebuild_cap[genombuild]}-"
                f"{mehari_tx_sources[genombuild][params.primary_tx_source]}-"
                f"{DV.mehari_tx}.bin.zst"
            )
            for genombuild in GENOMEBUILDS
        }

        with open(params.template) as tfh:
            template_content = tfh.read()

        with open(output.script, "w") as ofh:
            ofh.write(
                Template(template_content).render(
                    DV=DV,
                    PV=PV,
                    genomebuilds=GENOMEBUILDS,
                    mehari_tx_files=mehari_tx_files,
                    mehari_freq_versions=mehari_freq_versions,
                    downloader_base=workflow.basedir,
                )
            )
        shell("chmod +x {output.script}")
