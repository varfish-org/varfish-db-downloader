# Obtain current dump of HPO files


rule noref_hpo_download:
    output:
        txt="noref/hpo/latest/download/phenotype.hpoa",
    shell:
        r"""
        wget \
            -O {output.txt} \
            http://compbio.charite.de/jenkins/job/hpo.annotations.2018/lastSuccessfulBuild/artifact/misc_2018/phenotype.hpoa
        """


rule noref_hpo_name_mapping:
    output:
        tsv="noref/hpo/latest/HpoName.tsv",
        release_info="noref/hpo/latest/HpoName.release_info",
    run:
        import obonet
        import datetime

        graph = obonet.read_obo(
            "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo"
        )
        with open(output.tsv, "w") as fh:
            fh.write("hpo_id\tname\n")
            for _id, data in graph.nodes(data=True):
                fh.write("{}\t{}\n".format(_id, data["name"]))
        with open(output.release_info, "w") as fh:
            fh.write("table\tversion\tgenomebuild\tnull_value\n")
            fh.write(
                "HpoName\t{}\tGRCh38\t\n".format(datetime.datetime.today().strftime("%Y/%m/%d"))
            )


rule noref_hpo_to_tsv:
    input:
        header="header/hpo.txt",
        txt="noref/hpo/latest/download/phenotype.hpoa",
    output:
        tsv="noref/hpo/latest/Hpo.tsv",
        release_info="noref/hpo/latest/Hpo.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +6 {input.txt} \
            | awk -F $'\t' '
                BEGIN {{ OFS=FS }}
                {{
                    split($1,a,":")
                    print $0,a[1]=="OMIM"?a[2]:"",a[1]=="DECIPHER"?a[2]:"",a[1]=="ORPHA"?a[2]:""
                }}'
        ) \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nHpo\t$(date +%Y/%m/%d)\tGRCh38\t" > {output.release_info}
        """
