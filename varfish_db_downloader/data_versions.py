"""Declaration of data versions."""

import attrs


@attrs.frozen()
class DataVersions:
    #: String to use for GRCh37 ENSEMBL version.
    ensembl_37: str
    #: String to use for ENSEMBL version.
    ensembl: str
    #: String to use for current date.
    today: str
    #: Version of dbNSFP.
    dbnsfp: str
    #: Version of dbscSNV.
    dbscsnv: str
    #: Version of CADD.
    cadd: str
    #: Version of gnomAD for constraints.
    gnomad_constraints: str
    #: Version of gnomAD mtDNA.
    gnomad_mtdna: str
    #: Version of gnomAD v2.
    gnomad_v2: str
    #: Version of gnomAD v3.
    gnomad_v3: str
    #: Version of gnomAD SV.
    gnomad_sv: str
    #: Version of dbVar.
    dbvar: str
    #: Version of DGV.
    dgv: str
    #: Version of DGV Gold Standard.
    dgv_gs: str
    #: ExAC CNVs.
    exac_cnv: str
    #: Thousand Genomes SVs.
    g1k_svs: str
    #: UCSC conservation.
    ucsc_cons: str
    #: UCSC repeat masker.
    ucsc_rmsk: str
    #: UCSC genomicSuperDups
    ucsc_genomic_super_dups: str
    #: UCSC genome browser altSeqLiftOverPsl.
    ucsc_alt_seq_liftover: str
    #: UCSC genome browser fixSeqLiftOverPsl.
    ucsc_fix_seq_liftover: str
    #: RefSeq version (GRCh37)
    refseq_37: str
    #: ENSEMBL version (GRCh37)
    ensembl_37: str


#: The data versions to use.
DATA_VERSIONS = DataVersions(
    ensembl_37="87",
    ensembl="109",
    today="2023-06-01",
    dbnsfp="4.4",
    dbscsnv="1.1",
    cadd="1.6",
    gnomad_constraints="2.1.1",
    gnomad_mtdna="3.1",
    gnomad_v2="2.1.1",
    gnomad_v3="3.1.2",
    gnomad_sv="2.1.1",
    dbvar="2023-05-16",
    dgv="2020-02-25",
    dgv_gs="2016-05-15",
    exac_cnv="0.3.1",
    g1k_svs="phase3-v2",
    ucsc_cons="2016-10-07",
    ucsc_rmsk="2020-03-22",
    ucsc_genomic_super_dups="2011-10-25",
    ucsc_alt_seq_liftover="2020-03-22",
    ucsc_fix_seq_liftover="2020-05-24",
    refseq_37="105",
    ensembl_37="87",
)
