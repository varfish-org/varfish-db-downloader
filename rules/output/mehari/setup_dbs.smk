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

        with open(params.template) as tfh:
            template_content = tfh.read()

        with open(output.script, "w") as ofh:
            ofh.write(
                Template(template_content).render(
                    DV=DV,
                    PV=PV,
                    genomebuilds=GENOMEBUILDS,
                    primary_tx_source=params.primary_tx_source,
                    mehari_freq_versions=mehari_freq_versions,
                    downloader_base=workflow.basedir,
                )
            )
        shell("chmod +x {output.script}")
