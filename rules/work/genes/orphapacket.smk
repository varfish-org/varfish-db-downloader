## Rules related to annotating genes with ORDO terms


rule genes_orphapacket_download:  # -- download orphapacket file
    output:
        tar="work/download/genes/orphapacket/{version}/orphapacket.tar.gz",
    shell:
        r"""
        wget -O {output.tar} \
            https://github.com/Orphanet/orphapacket/archive/refs/tags/v10.1.tar.gz
        """


rule genes_orphapacket_diseases:  # -- postprocess file for HGNC gene IDs
    input:
        tar="work/download/genes/orphapacket/{version}/orphapacket.tar.gz",
        xlink="output/full/mehari/genes-xlink-{date}/genes-xlink.tsv",
    output:
        tsv="work/genes/orphapacket/{version}+{date}/orpha_diseases.tsv",
        tsv_md5="work/genes/orphapacket/{version}+{date}/orpha_diseases.tsv.md5",
    shell:
        """
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        tar -C $TMPDIR -xf $(readlink -f {input.tar})

        python ./scripts/genes-orpha-diseases.py {input.xlink} $TMPDIR/orphapacket-*/json \
        | qsv sort -d '\t' \
        | qsv fmt -t '\t' \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv}.md5
        """
