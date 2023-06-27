## Write mim2gene file for worker.


rule genes_ncbi_process_mim2gene:  # -- process NCBI MedGen mim2gene
    input:
        download="work/download/genes/ncbi/{date}/mim2gene_medgen",
    output:
        tsv="output/full/worker/mim2gene/{date}/mim2gene.tsv",
        tsv_md5="output/full/worker/mim2gene/{date}/mim2gene.tsv.md5",
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        awk -f scripts/genes-mim2gene.awk \
            -F $'\t' \
            {input.download} \
        > {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """
