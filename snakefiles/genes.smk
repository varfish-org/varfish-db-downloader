rule genes_xlink_ensembl:
    output:
        tsv="genes/xlink/ensembl.tsv",
        tsv_md5="genes/xlink/ensembl.tsv.md5",
    shell:
        r"""
        echo -e "ensembl_gene_id\tensembl_transcript_id\tentrez_id\tgene_symbol" >{output.tsv}

        wget --no-check-certificate \
            -O- \
            'https://ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "entrezgene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>' \
        | sort -u \
        >> {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """


rule genes_xlink_hgnc:
    output:
        download_json=temp("genes/xlink/download/hgnc/hgnc_complete_set.json"),
        tsv="genes/xlink/hgnc.tsv",
        tsv_md5="genes/xlink/hgnc.tsv.md5",
    shell:
        r"""
        set -x

        wget --no-check-certificate \
            -O {output.download_json} \
            https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc_complete_set.json

        jq \
            --raw-output \
            --from-file scripts/genes-xlink-hgnc.jq \
            {output.download_json} \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv_md5}
        """


rule genes_hgnc_info:
    input:
        download_json="genes/xlink/download/hgnc/hgnc_complete_set.json",
    output:
        tsv="genes/hgnc/hgnc_info.jsonl",
        tsv_md5="genes/hgnc/hgnc_info.jsonl.md5",
    shell:
        r"""
        set -x

        jq \
            --compact-output \
            --raw-output \
            --from-file scripts/genes-hgnc-info.jq \
            {input.download_json} \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv_md5}
        """


rule genes_gnomad_constraints_v2_1_1_download:
    output:
        bgz="genes/gnomad_constraints/download/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        bgz_md5="genes/gnomad_constraints/download/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz.md5",
    shell:
        r"""
        set -x

        wget --no-check-certificate \
            -O {output.bgz} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

        md5sum {output.bgz} >{output.bgz_md5}
        """


def run_genes_gnomad_constraints_v2_1_1_to_tsv(wildcards):
    """Extra function because of snakefmt issues."""
    columns_src = [
        "transcript",
        "exp_lof",
        "exp_mis",
        "exp_syn",
        "mis_z",
        "obs_lof",
        "obs_mis",
        "obs_syn",
        "oe_lof",
        "oe_lof_lower",
        "oe_lof_upper",
        "oe_mis",
        "oe_mis_lower",
        "oe_mis_upper",
        "oe_syn",
        "oe_syn_lower",
        "oe_syn_upper",
        "pLI",
        "syn_z",
        "exac_pLI",
        "exac_obs_lof",
        "exac_exp_lof",
        "exac_oe_lof",
    ]
    columns_src_str = ",".join(columns_src)
    columns_tmp = ["ensembl_transcript_id"] + columns_src[1:]
    columns_tmp_str = ",".join(columns_tmp)
    columns_dst = ["ensembl_gene_id", "entrez_id", "gene_symbol"] + columns_src[1:]
    columns_dst_str = ",".join(columns_dst)
    shell(
        r"""
                set -x

                zcat {input.bgz} \
                | tr '\t' ',' \
                > {output.txt_tmp}

                qsv select {columns_src_str} {output.txt_tmp} \
                | qsv rename {columns_tmp_str} \
                | qsv sort -u \
                | tr ',' '\t' \
                > {output.tsv_tmp}

                qsv join -d '\t' ensembl_transcript_id {output.tsv_tmp} ensembl_transcript_id {input.xlink_ensembl} \
                | qsv select {columns_dst_str} \
                | tr ',' '\t' \
                > {output.tsv}

                md5sum {output.tsv} >{output.tsv_md5}
            """
    )


rule genes_gnomad_constraints_v2_1_1_to_tsv:
    input:
        bgz="genes/gnomad_constraints/download/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        xlink_ensembl="genes/xlink/ensembl.tsv",
    output:
        txt_tmp=temp("genes/gnomad_constraints/download/gnomad.v2.1.1.lof_metrics.by_gene.txt"),
        tsv_tmp=temp("genes/gnomad_constraints/download/gnomad_constraints-subset.tsv"),
        tsv="genes/gnomad_constraints/gnomad_constraints.tsv",
        tsv_md5="genes/gnomad_constraints/gnomad_constraints.tsv.md5",
    run:
        run_genes_gnomad_constraints_v2_1_1_to_tsv(wildcards)


rule genes_mim2gene:
    output:
        download=temp("genes/mim2gene/download/mim2gene_medgen"),
        tsv="genes/mim2gene/mim2gene.tsv",
        tsv_md5="genes/mim2gene/mim2gene.tsv.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.download} \
            https://ftp.ncbi.nih.gov/gene/DATA/mim2gene_medgen

        awk -f scripts/genes-mim2gene.awk \
            -F $'\t' \
            {output.download} \
        > {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """


rule genes_gene_download:
    output:
        ags="genes/ncbi/download/Homo_sapiens.ags.gz",
        ags_md5="genes/ncbi/download/Homo_sapiens.ags.gz.md5",
        gene2xml="genes/ncbi/download/linux64.gene2xml",
        gene2xml_md5="genes/ncbi/download/linux64.gene2xml.md5",
    shell:
        r"""
        cd $(dirname {output.ags})

        if [[ ! -e Homo_sapiens.ags.gz.md5 ]]; then
            wget --no-check-certificate \
                -O Homo_sapiens.ags.gz \
                https://ftp.ncbi.nih.gov/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz
            md5sum Homo_sapiens.ags.gz >Homo_sapiens.ags.gz.md5
        fi

        if [[ ! -e linux64.gene2xml.md5 ]]; then
            wget --no-check-certificate \
                -O linux64.gene2xml.gz \
                https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/gene2xml/linux64.gene2xml.gz
            gzip -d -c linux64.gene2xml.gz > linux64.gene2xml
            chmod +x linux64.gene2xml
            md5sum linux64.gene2xml > linux64.gene2xml.md5
        fi
        """


rule result_noref_ncbi_gene_process:
    input:
        ags="genes/ncbi/download/Homo_sapiens.ags.gz",
        gene2xml="genes/ncbi/download/linux64.gene2xml",
    output:
        jsonl="genes/ncbi/gene_info.jsonl",
        jsonl_md5="genes/ncbi/gene_info.jsonl.md5",
    shell:
        r"""
        ./{input.gene2xml} -b T -c T -i {input.ags} \
        | python3 scripts/refseq_xml_to_json.py \
        > {output.jsonl}
        md5sum {output.jsonl} >{output.jsonl_md5}
        """


# For GRCh37, we use ucsc transcript ID instead of enst as the conservation
# file from UCSC uses these IDs.
rule genes_enst_ensg_grch37:
    output:
        download_txt=temp(
            "genes/enst_ensg/grch37/download/knowntoEnsembl.txt.gz",
        ),
        download_gtf=temp(
            "genes/enst_ensg/grch37/download/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
        ),
        tmp1="genes/enst_ensg/grch37/download/tmp1.txt",
        tmp2="genes/enst_ensg/grch37/download/tmp2.txt",
        tsv="genes/enst_ensg/grch37/enst_ensg.tsv",
        tsv_md5="genes/enst_ensg/grch37/enst_ensg.tsv.md5",
    shell:
        r"""
        set -x
        export LC_ALL=C

        wget --no-check-certificate \
            -O {output.download_gtf} \
            'https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'
        wget --no-check-certificate \
            -O {output.download_txt} \
            'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownToEnsembl.txt.gz'

        awk \
            -F $'\t' \
            -f scripts/genes-enst-ensg.awk \
            <(zcat {output.download_gtf}) \
        | sort \
        > {output.tmp1}
        zcat {output.download_txt} \
        | sed -e 's/\..//g' \
        | sort -k2,2 \
        >> {output.tmp2}

        echo -e "real_enst\tenst\tensg" > {output.tsv}
        join -t $'\t' -1 2 -2 1 {output.tmp2} {output.tmp1} \
        >> {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """


# We use the full dbNSFP genes information file.
rule genes_:
    input:
        tsv=f"annos/grch37/dbnsfp-{DBNSFP_VERSION}a/download/dbNSFP{DBNSFP_VERSION}_gene.complete.gz",
    output:
        tsv="genes/dbnsfp/genes.tsv.gz",
        tsv_md5="genes/dbnsfp/genes.tsv.gz.md5",
    shell:
        r"""
        set -x
        export LC_ALL=C

        zcat {input.tsv} \
        | pigz -c \
        > {output.tsv

        md5sum {output.tsv} >{output.tsv_md5}
        """
