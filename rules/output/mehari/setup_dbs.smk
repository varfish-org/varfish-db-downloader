rule generate_mehari_dbs_setup_script:  # -- Generate setup script for Mehari CLI
    output:
        script="output/full/mehari/setup_mehari_dbs.sh",
    params:
        template="rules/output/mehari/setup_mehari_dbs.sh.j2",
    run:
        from jinja2 import Template
        grch37 = "grch37" in GENOMEBUILDS
        grch38 = "grch38" in GENOMEBUILDS

        with open(params.template) as tfh:
            template_content = tfh.read()
        with open(output.script, "w") as ofh:
            ofh.write(Template(template_content).render(DV=DV, PV=PV, grch37=grch37, grch38=grch38, downloader_base=workflow.basedir))
        shell("chmod +x {output.script}")