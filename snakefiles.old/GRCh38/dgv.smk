rule grch38_dgv_2020_download:
    output:
        gff="GRCh38/DGV/2020/download/DGV.GS.hg38.gff3",
        gff_md5="GRCh38/DGV/2020/download/DGV.GS.hg38.gff3.md5",
        txt="GRCh38/DGV/2020/download/GRCh38_hg38_variants_2020-02-25.txt",
        txt_md5="GRCh38/DGV/2020/download/GRCh38_hg38_variants_2020-02-25.txt.md5",
    log:
        "GRCh38/DGV/2020/download/wget.log",
    shell:
        r"""
        pushd $(dirname {output.gff})
        wget --no-check-certificate \
            http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3 \
            -o $(basename {log})
        md5sum $(basename {output.gff}) > $(basename {output.gff_md5})
        popd

        pushd $(dirname {output.txt})
        wget --no-check-certificate \
            http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt \
            -o $(basename {log})
        md5sum $(basename {output.txt}) > $(basename {output.txt_md5})
        popd
        """


rule result_grch38_dgv_2020_tsv:
    input:
        header="header/dgvsvs.txt",
        txt="GRCh38/DGV/2020/download/GRCh38_hg38_variants_2020-02-25.txt",
    output:
        tsv="GRCh38/DGV/2020/DgvSvs.tsv",
        release_info="GRCh38/DGV/2020/DgvSvs.release_info",
    run:
        to_tsv(input.txt, output.tsv, output.release_info, input.header)


rule result_grch38_dgv_goldstandard_2020_tsv:
    input:
        header="header/dgvgoldstandardsvs.txt",
        gff="GRCh38/DGV/2020/download/DGV.GS.hg38.gff3",
    output:
        tsv="GRCh38/DGV/2020/DgvGoldStandardSvs.tsv",
        release_info="GRCh38/DGV/2020/DgvGoldStandardSvs.release_info",
    run:
        to_tsv(input.gff, output.tsv, output.release_info, input.header)
