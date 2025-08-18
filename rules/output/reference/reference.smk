
rule grch3x_reference:
    input:
        fa="work/reference/{genomebuild}/reference.fa",
        fai="work/reference/{genomebuild}/reference.fa.fai",
    output:
        fa="output/full/reference/{genomebuild}/reference.fa",
        fai="output/full/reference/{genomebuild}/reference.fa.fai",
    shell:
        r"""
        cp {input.fa} {output.fa}
        cp {input.fai} {output.fai}
        """
