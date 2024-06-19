## Conversion of ClinVar files for worker.


def input_annos_strucvar_clinvar_convert(wildcards):
    """Return input files for ``rule annos_strucvar_clinvar_convert``."""
    assert DV.clinvar_release.startswith(wildcards.clinvar_release), f"DV.clinvar_release={DV.clinvar_release}, wildcards.clinvar_release={wildcards.clinvar_release}"
    return {
        "jsonl": (
            f"work/download/annos/{wildcards.genome_release}/strucvars/"
            f"clinvar/{DV.clinvar_release}/clinvar-variants-{wildcards.genome_release}-strucvars.jsonl.gz"
        )
    }


rule annos_strucvar_clinvar_convert:
    input:
        unpack(input_annos_strucvar_clinvar_convert),
    output:
        bin=f"output/full/worker/clinvar-strucvars-{{genome_release}}-{{clinvar_release}}+{PV.worker}/clinvar-strucvars.bin",
        spec=f"output/full/worker/clinvar-strucvars-{{genome_release}}-{{clinvar_release}}+{PV.worker}/clinvar-strucvars.spec.yaml",
    wildcard_constraints:
        genome_release=RE_GENOME,
        clinvar_release=RE_VERSION,
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type clinvar-sv \
            --path-input {input.jsonl} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/clinvar_strucvars.spec.yaml \
            --value id_prefix=variant-server-worker/patho-mms \
            --value today={TODAY} \
            \
            --value genome_release={wildcards.genome_release} \
            --value clinvar_release={wildcards.clinvar_release} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """
