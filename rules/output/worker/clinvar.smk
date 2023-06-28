## Conversion of ClinVar files for worker.


def input_annos_strucvar_clinvar_convert(wildcards):
    """Return input files for ``rule annos_strucvar_clinvar_convert``."""
    clinvar_version = wildcards.clinvar_release.replace("-", "").split("+")[0]
    return {
        "tsv": (
            f"work/download/annos/{wildcards.genome_release}/strucvars/"
            f"clinvar/{clinvar_version}/clinvar_strucvar.tsv.gz"
        )
    }


rule annos_strucvar_clinvar_convert:
    input:
        unpack(input_annos_strucvar_clinvar_convert),
    output:
        bin=f"output/full/worker/clinvar-strucvars-{{genome_release}}-{{clinvar_release}}+{PV.worker}/clinvar-strucvars.bin",
    wildcard_constraints:
        genome_release=RE_GENOME,
        clinvar_release=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type clinvar-sv \
            --path-input {input.tsv} \
            --path-output-bin {output.bin}
        """
