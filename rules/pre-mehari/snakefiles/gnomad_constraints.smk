def input_result_grch3X_gnomad_constraints_tsv(wildcards):
    version = wildcards.version
    if wildcards.genome_build == "GRCh37":
        txt = f"work/download/genes/gnomad/{version}/gnomad.v{version}.lof_metrics.by_gene.txt.bgz",
    elif wildcards.genome_build == "GRCh38":
        txt = f"work/download/genes/gnomad/{version}/gnomad.v{version}.constraint_metrics.tsv",
    else:
        raise ValueError(f"Unsupported genome build: {wildcards.genome_build}")
    
    return {
        "header": "rules/pre-mehari/header/gnomadconstraints.txt",
        "txt": txt,
    }
    

rule result_grch3X_gnomad_constraints_tsv:
    input:
        unpack(input_result_grch3X_gnomad_constraints_tsv),
    output:
        tsv="output/pre-mehari/{genome_build}/gnomAD_constraints/v{version}/GnomadConstraints.tsv",
        release_info="output/pre-mehari/{genome_build}/gnomAD_constraints/v{version}/GnomadConstraints.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            zcat {input.txt} \
            | tail -n +2 \
            | sort -u \
            | awk -F $'\t' '
                BEGIN {{ OFS = FS }}
                {{
                    for (i=1; i<=NF; ++i) {{
                        if ($i == "NA") {{
                            $i = ""
                        }}
                    }}
                    print
                }}'
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nGnomadConstraints\tv{wildcards.version}\t{wildcards.genome_build}\t" > {output.release_info}
        """


# rule result_grch38_gnomad_constraints_tsv:
#     input:
#         txt="work/download/genes/gnomad/{version}/gnomad.v{version}.constraint_metrics.tsv",
#         header="rules/pre-mehari/header/gnomadconstraints.txt",
#     output:
#         tsv="output/pre-mehari/GRCh38/gnomAD_constraints/v{version}/GnomadConstraints.tsv",
#         release_info="output/pre-mehari/GRCh38/gnomAD_constraints/v{version}/GnomadConstraints.release_info",
#     shell:
#         r"""
#         (
#             cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
#             zcat {input.txt} \
#             | tail -n +2 \
#             | sort -u \
#             | awk -F $'\t' '
#                 BEGIN {{ OFS = FS }}
#                 {{
#                     for (i=1; i<=NF; ++i) {{
#                         if ($i == "NA") {{
#                             $i = ""
#                         }}
#                     }}
#                     print
#                 }}'
#         ) > {output.tsv}

#         echo -e "table\tversion\tgenomebuild\tnull_value\nGnomadConstraints\tv{wildcards.version}\tGRCh38\t" > {output.release_info}
#         """
