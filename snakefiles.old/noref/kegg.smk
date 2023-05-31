# Obtain dump of KEGG information from UCSC.

# Note in retrieving data from hgTables with wget:
# Just passing the parameters from the last submit form results in an error.
# One has to pass one of the parameters from the first form as well, which is normally stored in the sessions.
# As we don't access the sessions, there is the undocumented parameter to be appended to the call to make it work:
# hgta_fieldSelectTable=hg19.<value from table selection>
# To find <value from table selection> go into to the source code and find the correct value for the dropdown list 'hgta_table'


rule grchXX_kegg_download_ensembltokegg:
    output:
        tsv="noref/kegg/april2011/download/EnsemblToKegg.tsv",
    shell:
        r"""
        wget --no-check-certificate -O - 'https://genome-euro.ucsc.edu/cgi-bin/hgTables' --post-data 'hgta_database=hg19&hgta_table=keggPathway&hgta_fieldSelectTable=hg19.keggPathway&boolshad.hgta_fs.check.hg19.keggPathway.kgID=0&boolshad.hgta_fs.check.hg19.keggPathway.locusID=0&hgta_fs.check.hg19.keggPathway.mapID=on&boolshad.hgta_fs.check.hg19.keggPathway.mapID=0&hgta_doPrintSelectedFields=get+output&jsh_pageVertPos=0&hgta_fs.check.hg19.ensGtp.gene=on&boolshad.hgta_fs.check.hg19.ensGtp.gene=0&boolshad.hgta_fs.check.hg19.ensGtp.transcript=0&boolshad.hgta_fs.check.hg19.ensGtp.protein=0&boolshad.hgta_fs.check.hg19.knownToEnsembl.name=0&boolshad.hgta_fs.check.hg19.knownToEnsembl.value=0&boolshad.hgta_fs.linked.hg19.bioCycPathway=0&boolshad.hgta_fs.linked.hg19.ccdsKgMap=0&boolshad.hgta_fs.linked.hg19.ceBlastTab=0&boolshad.hgta_fs.linked.hg19.dmBlastTab=0&boolshad.hgta_fs.linked.hg19.drBlastTab=0&boolshad.hgta_fs.linked.hg19.ensGene=0&hgta_fs.linked.hg19.ensGtp=on&boolshad.hgta_fs.linked.hg19.ensGtp=0&boolshad.hgta_fs.linked.hg19.ensPep=0&boolshad.hgta_fs.linked.hg19.ensemblSource=0&boolshad.hgta_fs.linked.hg19.ensemblToGeneName=0&boolshad.hgta_fs.linked.hg19.foldUtr3=0&boolshad.hgta_fs.linked.hg19.foldUtr5=0&boolshad.hgta_fs.linked.hg19.gnfAtlas2Distance=0&boolshad.hgta_fs.linked.hg19.gnfU95Distance=0&boolshad.hgta_fs.linked.hg19.humanHprdP2P=0&boolshad.hgta_fs.linked.hg19.humanVidalP2P=0&boolshad.hgta_fs.linked.hg19.humanWankerP2P=0&boolshad.hgta_fs.linked.hg19.keggMapDesc=0&boolshad.hgta_fs.linked.hg19.kg5ToKg6=0&boolshad.hgta_fs.linked.hg19.kg6ToKg7=0&boolshad.hgta_fs.linked.hg19.kgAlias=0&boolshad.hgta_fs.linked.hg19.kgColor=0&boolshad.hgta_fs.linked.hg19.kgProtAlias=0&boolshad.hgta_fs.linked.hg19.kgProtMap2=0&boolshad.hgta_fs.linked.hg19.kgSpAlias=0&boolshad.hgta_fs.linked.hg19.kgTargetAli=0&boolshad.hgta_fs.linked.hg19.kgTxInfo=0&boolshad.hgta_fs.linked.hg19.kgXref=0&boolshad.hgta_fs.linked.hg19.knownBlastTab=0&boolshad.hgta_fs.linked.hg19.knownCanonical=0&boolshad.hgta_fs.linked.hg19.knownGene=0&boolshad.hgta_fs.linked.hg19.knownGeneMrna=0&boolshad.hgta_fs.linked.hg19.knownGenePep=0&boolshad.hgta_fs.linked.hg19.knownGeneTxMrna=0&boolshad.hgta_fs.linked.hg19.knownGeneTxPep=0&boolshad.hgta_fs.linked.hg19.knownIsoforms=0&boolshad.hgta_fs.linked.hg19.knownToAllenBrain=0&hgta_fs.linked.hg19.knownToEnsembl=on&boolshad.hgta_fs.linked.hg19.knownToEnsembl=0&boolshad.hgta_fs.linked.hg19.knownToGnfAtlas2=0&boolshad.hgta_fs.linked.hg19.knownToHInv=0&boolshad.hgta_fs.linked.hg19.knownToHprd=0&boolshad.hgta_fs.linked.hg19.knownToKeggEntrez=0&boolshad.hgta_fs.linked.hg19.knownToLocusLink=0&boolshad.hgta_fs.linked.hg19.knownToLynx=0&boolshad.hgta_fs.linked.hg19.knownToNextProt=0&boolshad.hgta_fs.linked.hg19.knownToPfam=0&boolshad.hgta_fs.linked.hg19.knownToRefSeq=0&boolshad.hgta_fs.linked.hg19.knownToSuper=0&boolshad.hgta_fs.linked.hg19.knownToTreefam=0&boolshad.hgta_fs.linked.hg19.knownToU133=0&boolshad.hgta_fs.linked.hg19.knownToU133Plus2=0&boolshad.hgta_fs.linked.hg19.knownToU95=0&boolshad.hgta_fs.linked.hg19.knownToVisiGene=0&boolshad.hgta_fs.linked.hg19.knownToWikipedia=0&boolshad.hgta_fs.linked.hg19.mmBlastTab=0&boolshad.hgta_fs.linked.hg19.omim2gene=-2&boolshad.hgta_fs.linked.hg19.rnBlastTab=0&boolshad.hgta_fs.linked.hg19.scBlastTab=0&boolshad.hgta_fs.linked.hg19.ucscScop=0&boolshad.hgta_fs.linked.hgFixed.refLink=0' \
        > {output.tsv}
        """


rule grchXX_kegg_download_refseqtokegg:
    output:
        tsv="noref/kegg/april2011/download/RefseqToKegg.tsv",
    shell:
        r"""
        wget -O - 'https://genome-euro.ucsc.edu/cgi-bin/hgTables' --post-data 'hgta_database=hg19&hgta_table=keggPathway&boolshad.hgta_fs.check.hg19.keggPathway.kgID=0&hgta_fieldSelectTable=hg19.keggPathway&boolshad.hgta_fs.check.hg19.keggPathway.locusID=0&hgta_fs.check.hg19.keggPathway.mapID=on&boolshad.hgta_fs.check.hg19.keggPathway.mapID=0&hgta_doPrintSelectedFields=get+output&jsh_pageVertPos=0&boolshad.hgta_fs.check.hg19.kgXref.kgID=0&boolshad.hgta_fs.check.hg19.kgXref.mRNA=0&boolshad.hgta_fs.check.hg19.kgXref.spID=0&boolshad.hgta_fs.check.hg19.kgXref.spDisplayID=0&boolshad.hgta_fs.check.hg19.kgXref.geneSymbol=0&boolshad.hgta_fs.check.hg19.kgXref.refseq=0&boolshad.hgta_fs.check.hg19.kgXref.protAcc=0&boolshad.hgta_fs.check.hg19.kgXref.description=0&boolshad.hgta_fs.check.hg19.kgXref.rfamAcc=0&boolshad.hgta_fs.check.hg19.kgXref.tRnaName=0&boolshad.hgta_fs.check.hg19.knownToKeggEntrez.name=0&boolshad.hgta_fs.check.hg19.knownToKeggEntrez.value=0&boolshad.hgta_fs.check.hg19.knownToKeggEntrez.keggEntrez=0&boolshad.hgta_fs.check.hg19.knownToRefSeq.name=0&boolshad.hgta_fs.check.hg19.knownToRefSeq.value=0&boolshad.hgta_fs.check.hgFixed.refLink.name=0&boolshad.hgta_fs.check.hgFixed.refLink.product=0&boolshad.hgta_fs.check.hgFixed.refLink.mrnaAcc=0&boolshad.hgta_fs.check.hgFixed.refLink.protAcc=0&boolshad.hgta_fs.check.hgFixed.refLink.geneName=0&boolshad.hgta_fs.check.hgFixed.refLink.prodName=0&hgta_fs.check.hgFixed.refLink.locusLinkId=on&boolshad.hgta_fs.check.hgFixed.refLink.locusLinkId=0&boolshad.hgta_fs.check.hgFixed.refLink.omimId=0&boolshad.hgta_fs.linked.go.goaPart=0&boolshad.hgta_fs.linked.hg19.allenBrainGene=0&boolshad.hgta_fs.linked.hg19.bioCycPathway=0&boolshad.hgta_fs.linked.hg19.ccdsInfo=0&boolshad.hgta_fs.linked.hg19.ccdsKgMap=0&boolshad.hgta_fs.linked.hg19.ceBlastTab=0&boolshad.hgta_fs.linked.hg19.cgapAlias=0&boolshad.hgta_fs.linked.hg19.dmBlastTab=0&boolshad.hgta_fs.linked.hg19.drBlastTab=0&boolshad.hgta_fs.linked.hg19.foldUtr3=0&boolshad.hgta_fs.linked.hg19.foldUtr5=0&boolshad.hgta_fs.linked.hg19.gad=0&boolshad.hgta_fs.linked.hg19.gadAll=0&boolshad.hgta_fs.linked.hg19.geneNetworkId=0&boolshad.hgta_fs.linked.hg19.geneReviews=0&boolshad.hgta_fs.linked.hg19.geneReviewsDetail=0&boolshad.hgta_fs.linked.hg19.gnfAtlas2Distance=0&boolshad.hgta_fs.linked.hg19.gnfU95Distance=0&boolshad.hgta_fs.linked.hg19.humanHprdP2P=0&boolshad.hgta_fs.linked.hg19.humanVidalP2P=0&boolshad.hgta_fs.linked.hg19.humanWankerP2P=0&boolshad.hgta_fs.linked.hg19.keggMapDesc=0&boolshad.hgta_fs.linked.hg19.kg5ToKg6=0&boolshad.hgta_fs.linked.hg19.kg6ToKg7=0&boolshad.hgta_fs.linked.hg19.kgAlias=0&boolshad.hgta_fs.linked.hg19.kgColor=0&boolshad.hgta_fs.linked.hg19.kgProtAlias=0&boolshad.hgta_fs.linked.hg19.kgProtMap2=0&boolshad.hgta_fs.linked.hg19.kgSpAlias=0&boolshad.hgta_fs.linked.hg19.kgTargetAli=0&boolshad.hgta_fs.linked.hg19.kgTxInfo=0&hgta_fs.linked.hg19.kgXref=on&boolshad.hgta_fs.linked.hg19.kgXref=0&boolshad.hgta_fs.linked.hg19.knownBlastTab=0&boolshad.hgta_fs.linked.hg19.knownCanonical=0&boolshad.hgta_fs.linked.hg19.knownGene=0&boolshad.hgta_fs.linked.hg19.knownGeneMrna=0&boolshad.hgta_fs.linked.hg19.knownGenePep=0&boolshad.hgta_fs.linked.hg19.knownGeneTxMrna=0&boolshad.hgta_fs.linked.hg19.knownGeneTxPep=0&boolshad.hgta_fs.linked.hg19.knownIsoforms=0&boolshad.hgta_fs.linked.hg19.knownToAllenBrain=0&boolshad.hgta_fs.linked.hg19.knownToEnsembl=0&boolshad.hgta_fs.linked.hg19.knownToGnfAtlas2=0&boolshad.hgta_fs.linked.hg19.knownToHInv=0&boolshad.hgta_fs.linked.hg19.knownToHprd=0&hgta_fs.linked.hg19.knownToKeggEntrez=on&boolshad.hgta_fs.linked.hg19.knownToKeggEntrez=0&boolshad.hgta_fs.linked.hg19.knownToLocusLink=0&boolshad.hgta_fs.linked.hg19.knownToLynx=0&boolshad.hgta_fs.linked.hg19.knownToNextProt=0&boolshad.hgta_fs.linked.hg19.knownToPfam=0&hgta_fs.linked.hg19.knownToRefSeq=on&boolshad.hgta_fs.linked.hg19.knownToRefSeq=0&boolshad.hgta_fs.linked.hg19.knownToSuper=0&boolshad.hgta_fs.linked.hg19.knownToTreefam=0&boolshad.hgta_fs.linked.hg19.knownToU133=0&boolshad.hgta_fs.linked.hg19.knownToU133Plus2=0&boolshad.hgta_fs.linked.hg19.knownToU95=0&boolshad.hgta_fs.linked.hg19.knownToVisiGene=0&boolshad.hgta_fs.linked.hg19.knownToWikipedia=0&boolshad.hgta_fs.linked.hg19.lsSnpPdb=0&boolshad.hgta_fs.linked.hg19.mmBlastTab=0&boolshad.hgta_fs.linked.hg19.mrnaOrientInfo=0&boolshad.hgta_fs.linked.hg19.omim2gene=-2&boolshad.hgta_fs.linked.hg19.omimAv=-2&boolshad.hgta_fs.linked.hg19.omimGene2=-2&boolshad.hgta_fs.linked.hg19.omimGeneMap=-2&boolshad.hgta_fs.linked.hg19.omimGeneSymbol=-2&boolshad.hgta_fs.linked.hg19.omimPhenotype=-2&boolshad.hgta_fs.linked.hg19.refFlat=0&boolshad.hgta_fs.linked.hg19.refGene=0&boolshad.hgta_fs.linked.hg19.refSeqAli=0&boolshad.hgta_fs.linked.hg19.rnBlastTab=0&boolshad.hgta_fs.linked.hg19.scBlastTab=0&boolshad.hgta_fs.linked.hg19.spMrna=0&boolshad.hgta_fs.linked.hg19.ucscScop=0&boolshad.hgta_fs.linked.hgFixed.ctdSorted=0&boolshad.hgta_fs.linked.hgFixed.gbCdnaInfo=0&boolshad.hgta_fs.linked.hgFixed.gbSeq=0&boolshad.hgta_fs.linked.hgFixed.geneName=0&boolshad.hgta_fs.linked.hgFixed.ggGeneClass=0&boolshad.hgta_fs.linked.hgFixed.ggGeneName=0&boolshad.hgta_fs.linked.hgFixed.ggLink=0&boolshad.hgta_fs.linked.hgFixed.imageClone=0&boolshad.hgta_fs.linked.hgFixed.productName=0&hgta_fs.linked.hgFixed.refLink=on&boolshad.hgta_fs.linked.hgFixed.refLink=0&boolshad.hgta_fs.linked.hgFixed.refSeqStatus=0&boolshad.hgta_fs.linked.hgFixed.refSeqSummary=0&boolshad.hgta_fs.linked.proteome.hgncXref=0&boolshad.hgta_fs.linked.proteome.spOldNew=0&boolshad.hgta_fs.linked.proteome.spReactomeEvent=0&boolshad.hgta_fs.linked.proteome.spReactomeId=0&boolshad.hgta_fs.linked.uniProt.accToKeyword=0&boolshad.hgta_fs.linked.uniProt.accToTaxon=0&boolshad.hgta_fs.linked.uniProt.citation=0&boolshad.hgta_fs.linked.uniProt.comment=0&boolshad.hgta_fs.linked.uniProt.description=0&boolshad.hgta_fs.linked.uniProt.displayId=0&boolshad.hgta_fs.linked.uniProt.extDbRef=0&boolshad.hgta_fs.linked.uniProt.feature=0&boolshad.hgta_fs.linked.uniProt.gene=0&boolshad.hgta_fs.linked.uniProt.geneLogic=0&boolshad.hgta_fs.linked.uniProt.info=0&boolshad.hgta_fs.linked.uniProt.otherAcc=0&boolshad.hgta_fs.linked.uniProt.protein=0&boolshad.hgta_fs.linked.visiGene.gene=0' \
        > {output.tsv}
        """


rule grchXX_kegg_download_kegginfo:
    output:
        tsv="noref/kegg/april2011/download/KeggInfo.tsv",
    shell:
        r"""
        wget -O - 'https://genome-euro.ucsc.edu/cgi-bin/hgTables' --post-data 'hgta_database=hg19&hgta_table=keggMapDesc&hgta_fieldSelectTable=hg19.keggMapDesc&hgta_fs.check.hg19.keggMapDesc.mapID=on&boolshad.hgta_fs.check.hg19.keggMapDesc.mapID=0&hgta_fs.check.hg19.keggMapDesc.description=on&boolshad.hgta_fs.check.hg19.keggMapDesc.description=0&hgta_doPrintSelectedFields=get+output&jsh_pageVertPos=0&boolshad.hgta_fs.linked.hg19.keggPathway=0&boolshad.hgta_fs.linked.hg19.knownToKeggEntrez=0' \
        > {output.tsv}
        """


rule result_grchXX_kegg_to_tsv:
    input:
        header_ensembl_to_kegg="header/ensembltokegg.txt",
        header_refseq_to_kegg="header/refseqtokegg.txt",
        header_kegg_info="header/kegginfo.txt",
        ensembl_to_kegg="noref/kegg/april2011/download/EnsemblToKegg.tsv",
        refseq_to_kegg="noref/kegg/april2011/download/RefseqToKegg.tsv",
        kegg_info="noref/kegg/april2011/download/KeggInfo.tsv",
    output:
        ensembl_to_kegg="noref/kegg/april2011/EnsemblToKegg.tsv",
        refseq_to_kegg="noref/kegg/april2011/RefseqToKegg.tsv",
        kegg_info="noref/kegg/april2011/KeggInfo.tsv",
        release_info_kegginfo="noref/kegg/april2011/KeggInfo.release_info",
        release_info_ensembltokegg="noref/kegg/april2011/EnsemblToKegg.release_info",
        release_info_refseqtokegg="noref/kegg/april2011/RefseqToKegg.release_info",
    shell:
        r"""
        (
            cat {input.header_ensembl_to_kegg} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.ensembl_to_kegg} \
            | sort -u \
            | awk -F $'\t' '
                BEGIN {{ OFS=FS; }}
                {{
                    if ($2 == "n/a") {{
                        next;
                    }}
                    print $2, $1;
                }}'
        ) \
        > {output.ensembl_to_kegg}

        (
            cat {input.header_refseq_to_kegg} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.refseq_to_kegg} \
            | sort -u \
            | sed 's/,$//' \
            | awk -F $'\t' '
                BEGIN {{ OFS=FS; }}
                {{
                    if ($2 == "n/a") {{
                        next;
                    }}
                    print $2, $1;
                }}'
        ) \
        > {output.refseq_to_kegg}

        (
            cat {input.header_kegg_info} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.kegg_info}
        ) \
        > {output.kegg_info}

        echo -e "table\tversion\tgenomebuild\tnull_value\nKeggInfo\tApril 2011 (via UCSC Genome Browser)\t\t" > {output.release_info_kegginfo}
        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblToKegg\tApril 2011 (via UCSC Genome Browser)\t\t" > {output.release_info_ensembltokegg}
        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToKegg\tApril 2011 (via UCSC Genome Browser)\t\t" > {output.release_info_refseqtokegg}
        """
