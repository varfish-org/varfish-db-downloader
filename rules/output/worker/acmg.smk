## Rules to prepare acmg listing for worker.


rule acmg_prepare_worker:
    input:
        tsv="data/acmg_sf/{v_acmg_sf}/acmg_sf.tsv",
        spec="data/acmg_sf/{v_acmg_sf}/acmg_sf.spec.yaml",
    output:
        tsv=f"output/full/worker/acmg-sf-{{v_acmg_sf}}+{PV.worker}/acmg_sf.tsv",
        spec=f"output/full/worker/acmg-sf-{{v_acmg_sf}}+{PV.worker}/acmg_sf.spec.yaml",
    shell:
        r"""
        cp {input.tsv} {output.tsv}
        cp {input.spec} {output.spec}
        """
