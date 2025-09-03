rule generate_mehari_dbs_setup_script:  # -- Generate setup script for Mehari CLI
    output:
        script="output/full/mehari/setup_mehari_dbs.sh",
    params:
        template="rules/output/mehari/setup_mehari_dbs.sh.j2",
        primary_tx_source='refseq', #Choice of: refseq, ensembl, ensembl-and-refseq
    run:
        from jinja2 import Template
        grch37 = "grch37" in GENOMEBUILDS
        grch38 = "grch38" in GENOMEBUILDS

        freq_version_37 = gnomad_versions['grch37']+gnomad_versions['grch37']+DV.gnomad_mtdna+DV.helixmtdb+PV.annonars
        freq_version_38 = gnomad_versions['grch38'] + gnomad_versions['grch38'] + DV.gnomad_mtdna + DV.helixmtdb + PV.annonars

        with open(params.template) as tfh:
            template_content = tfh.read()
        with open(output.script, "w") as ofh:
            ofh.write(Template(template_content).render(
                DV=DV, PV=PV, grch37=grch37, grch38=grch38, primary_tx_source=params.primary_tx_source,
                freq_version_37=freq_version_37, freq_version_38=freq_version_38,
                downloader_base=workflow.basedir)
            )
        shell("chmod +x {output.script}")
