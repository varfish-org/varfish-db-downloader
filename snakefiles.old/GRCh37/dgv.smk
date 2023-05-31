rule grch37_dgv_2016_download:
    output:
        gff="GRCh37/DGV/2016/download/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
        gff_md5=(
            "GRCh37/DGV/2016/download/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3.md5"
        ),
        txt="GRCh37/DGV/2020/download/GRCh37_hg19_variants_2020-02-25.txt",
        txt_md5="GRCh37/DGV/2020/download/GRCh37_hg19_variants_2020-02-25.txt.md5",
    log:
        "GRCh37/DGV/2016/download/wget.log",
    shell:
        r"""
        pushd $(dirname {output.gff})
        wget --no-check-certificate \
            http://dgv.tcag.ca/dgv/docs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3 \
            -o $(basename {log})
        md5sum $(basename {output.gff}) > $(basename {output.gff_md5})
        popd

        pushd $(dirname {output.txt})
        wget --no-check-certificate \
            http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt \
            -o $(basename {log})
        md5sum $(basename {output.txt}) > $(basename {output.txt_md5})
        popd
        """


rule result_grch37_dgv_2020_tsv:
    input:
        header="header/dgvsvs.txt",
        txt="GRCh37/DGV/2020/download/GRCh37_hg19_variants_2020-02-25.txt",
    output:
        tsv="GRCh37/DGV/2020/DgvSvs.tsv",
        release_info="GRCh37/DGV/2020/DgvSvs.release_info",
    run:
        to_tsv(input.txt, output.tsv, output.release_info, input.header)


rule result_grch37_dgv_goldstandard_2016_tsv:
    input:
        header="header/dgvgoldstandardsvs.txt",
        gff="GRCh37/DGV/2016/download/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
    output:
        tsv="GRCh37/DGV/2016/DgvGoldStandardSvs.tsv",
        release_info="GRCh37/DGV/2016/DgvGoldStandardSvs.release_info",
    run:
        to_tsv(input.gff, output.tsv, output.release_info, input.header)
