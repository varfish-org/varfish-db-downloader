## Rules to create annonars RocksDB for (Wetzel & Darbro, 2022) data.


rule output_worker_patho_mms:
    output:
        bed="output/worker/annos/strucvars/patho-mms-{genome_release}-{v_patho_mms}/patho-mms.bed",
        bed_md5="output/worker/annos/strucvars/patho-mms-{genome_release}-{v_patho_mms}/patho-mms.bed.md5",
        spec="output/worker/annos/strucvars/patho-mms-{genome_release}-{v_patho_mms}/patho-smms.bed.spec.json",
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_patho_mms=RE_VERSION,
    shell:
        r"""
        base=data/patho-mms/{wildcards.v_patho_mms}/patho-mms-{wildcards.genome_release}

        cp $base.bed           {output.bed}
        cp $base.bed.md5       {output.bed_md5}
        cp $base.bed.spec.json {output.spec}
        """
