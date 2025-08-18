rule result_noref_hpo_name_mapping:
    input:
        obo="work/download/hpo/{v_hpo}/hp.obo",
    output:
        tsv="output/pre-mehari/noref/hpo/{v_hpo}/HpoName.tsv",
        release_info="output/pre-mehari/noref/hpo/{v_hpo}/HpoName.release_info",
    run:
        import obonet

        graph = obonet.read_obo(input.obo)
        with open(output.tsv, "w") as fh:
            fh.write("hpo_id\tname\n")
            for _id, data in graph.nodes(data=True):
                fh.write("{}\t{}\n".format(_id, data["name"]))
        with open(output.release_info, "w") as fh:
            fh.write("table\tversion\tgenomebuild\tnull_value\n")
            fh.write(f"HpoName\t{wildcards.v_hpo}\t\t\n")


rule result_noref_hpo_to_tsv:
    input:
        header="rules/pre-mehari/header/hpo.txt",
        txt="work/download/hpo/{v_hpo}/phenotype.hpoa",
    output:
        tsv="output/pre-mehari/noref/hpo/{v_hpo}/Hpo.tsv",
        release_info="output/pre-mehari/noref/hpo/{v_hpo}/Hpo.release_info",
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

        echo -e "table\tversion\tgenomebuild\tnull_value\nHpo\t{wildcards.v_hpo}\t\t" > {output.release_info}
        """
