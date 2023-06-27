## Rules to prepare acmg listing for worker.


rule acmg_prepare_worker:
    input:
        acmg_sf="data/acmg/{v_acmg_sf}/acmg.tsv",
    output:
        tsv="output/full/worker/acmg-sf-{v_acmg_sf}/acmg_sf.tsv",
    shell:
        r"""
        cp {input.acmg_sf} {output.tsv}
        """
