#: Maximal distance (in bp) from exon to be considered "near coding".
NEAR_CODING_DIST = 1000


# Create BED file with "near coding regions", based on RefSeq.
rule annos_near_coding_regions:
    input:
        bed="features/{genome_build}/gene_regions/refseq.bed.gz",
    output:
        bed="annos/{genome_build}/near_coding/near_coding.bed",
        bed_gz="annos/{genome_build}/near_coding/near_coding.bed.bgz",
        bed_gz_tbi="annos/{genome_build}/near_coding/near_coding.bed.bgz.tbi",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        zcat {input.bed} \
        | tail -n +2 \
        | bedops --range {NEAR_CODING_DIST} --everything - \
        > /tmp/regions.bed

        bedops --merge /tmp/regions.bed \
        | sort-bed - \
        > {output.bed}

        bgzip -c {output.bed} >{output.bed_gz}
        tabix -p bed -f {output.bed_gz}
        """


rule annos_helixmtdb_download:
    output:
        tsv="annos/{genome_build}/helixmtdb/download/helixmtdb.tsv",
    shell:
        r"""
        wget \
            --no-check-certificate \
            -O {output} \
            https://helix-research-public.s3.amazonaws.com/mito/HelixMTdb_20200327.tsv
        """


rule annos_helixmtdb_convert:
    input:
        tsv="annos/{genome_build}/helixmtdb/download/helixmtdb.tsv",
    output:
        vcf="annos/{genome_build}/helixmtdb/helixmtdb.vcf.gz",
        vcf_tbi="annos/{genome_build}/helixmtdb/helixmtdb.vcf.gz.tbi",
    shell:
        r"""
        cat {input.tsv} \
        | python3  scripts/helix-to-vcf.py \
        > {output.vcf}.tmp

        if [[ {wildcards.genome_build} == GRCh37 ]]; then
            sed -e 's/chrM/MT/g' {output.vcf}.tmp \
            | bgzip -c >
            > {output.vcf}
        else
            bgzip -c {output.vcf}.tmp >{output.vcf}
        fi

        tabix -f {output.vcf}

        rm -f {output.vcf}.tmp
        """


rule annos_gnomad_mtdna:
    output:
        dl="annos/{genome_build}/gnomad_mtdna/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
        vcf="annos/{genome_build}/gnomad_mtdna/gnomad_mtdna.vcf.gz",
        vcf_tbi="annos/{genome_build}/gnomad_mtdna/gnomad_mtdna.vcf.gz.tbi",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.dl} \
            https://datasetgnomad.blob.core.windows.net/dataset/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz

        if [[ {wildcards.genome_build} == grch37 ]]; then
            zcat {output.dl} \
            | sed \
                -e 's/chrM/MT/g' \
                -e 's/GRCh38_MT/GRCh37/g'  \
                -e 's/,GERP_DIST/\&GERP_DIST/g' \
                -e 's/,BP_DIST/\&BP_DIST/g' \
                -e 's/,DIST_FROM_LAST_EXON/\&DIST_FROM_LAST_EXON/g' \
                -e 's/,50_BP_RULE/\&50_BP_RULE/g' \
                -e 's/,PHYLOCSF_TOO_SHORT/\&PHYLOCSF_TOO_SHORT/g' \
            | bgzip -c \
            > {output.vcf}
        else
            cp {output.dl} {output.vcf}
        fi

        tabix -f {output.vcf}
        """


# GNOMAD_PREFIX = "https://datasetgnomad.blob.core.windows.net/dataset/release"
GNOMAD_PREFIX = "https://gnomad-public-us-east-1.s3.amazonaws.com/release"
GNOMAD_V3 = "3.1.2"
GNOMAD_V2 = "2.1.1"


rule annos_gnomad_nuclear_download_2:
    output:
        vcf="annos/grch37/gnomad_{kind}/download/gnomad.{kind}.r{version}.sites.{chrom}.vcf.bgz",
        vcf_tbi="annos/grch37/gnomad_{kind}/download/gnomad.{kind}.r{version}.sites.{chrom}.vcf.bgz.tbi",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.vcf.bgz
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf_tbi} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.vcf.bgz.tbi
        """


rule annos_gnomad_nuclear_download_liftover_2:
    output:
        vcf="annos/grch38/gnomad_{kind}/download/gnomad.{kind}.r{version}.sites.{chrom}.liftover_grch38.vcf.bgz",
        vcf_tbi="annos/grch38/gnomad_{kind}/download/gnomad.{kind}.r{version}.sites.{chrom}.liftover_grch38.vcf.bgz.tbi",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/liftover_grch38/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.liftover_grch38.vcf.bgz
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf_tbi} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/liftover_grch38/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.liftover_grch38.vcf.bgz.tbi
        """


rule annos_gnomad_nuclear_download_3:
    output:
        vcf="annos/grch38/gnomad_{kind}/download/gnomad.{kind}.v{version}.sites.chr{chrom}.vcf.bgz",
        vcf_tbi="annos/grch38/gnomad_{kind}/download/gnomad.{kind}.v{version}.sites.chr{chrom}.vcf.bgz.tbi",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.v{wildcards.version}.sites.chr{wildcards.chrom}.vcf.bgz
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf_tbi} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.v{wildcards.version}.sites.chr{wildcards.chrom}.vcf.bgz.tbi
        """


def input_annos_gnomad_grch37(wildcards):
    chroms = list(range(1, 23)) + ["X"]
    # chrY is only available for GRCh37 genomes
    if wildcards.kind == "exomes":
        chroms.append("Y")
    tpl = "annos/grch37/gnomad_{kind}/download/gnomad.{kind}.r{version}.sites.{chrom}.vcf.bgz"
    return [tpl.format(kind=wildcards.kind, version=GNOMAD_V2, chrom=chrom) for chrom in chroms]


rule annos_gnomad_grch37:
    input:
        input_annos_gnomad_grch37,
    output:
        touch("annos/grch37/gnomad_{kind}/.done"),


def input_annos_gnomad_grch38(wildcards):
    chroms = list(range(1, 23)) + ["X", "Y"]
    if wildcards.kind == "exomes":
        tpl = "annos/grch38/gnomad_{kind}/download/gnomad.{kind}.r{version}.sites.{chrom}.liftover_grch38.vcf.bgz"
        return [tpl.format(kind=wildcards.kind, version=GNOMAD_V2, chrom=chrom) for chrom in chroms]
    else:
        tpl = (
            "annos/grch38/gnomad_{kind}/download/gnomad.{kind}.v{version}.sites.chr{chrom}.vcf.bgz"
        )
        return [tpl.format(kind=wildcards.kind, version=GNOMAD_V3, chrom=chrom) for chrom in chroms]


rule annos_gnomad_grch38:
    input:
        input_annos_gnomad_grch38,
    output:
        touch("annos/grch38/gnomad_{kind}/.done"),


rule annos_ucsc_conservation_download:
    output:
        fa="annos/{genome_build}/ucsc_conservation/download/knownGene.exonAA.fa.gz",
    shell:
        r"""
        if [[ {wildcards.genome_build} == grch37 ]]; then
            ucsc_name=hg19
        else
            ucsc_name=hg38
        fi

        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.fa} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            "https://hgdownload.cse.ucsc.edu/goldenpath/${{ucsc_name}}/multiz100way/alignments/knownGene.exonAA.fa.gz"
        """


rule annos_ucsc_conservation_to_vcf:
    input:
        hgnc="genes/hgnc/hgnc_info.jsonl",
        enst_ensg="genes/enst_ensg/{genome_build}/enst_ensg.tsv",
        reference="reference/{genome_build}/reference/reference.fa",
        fa="annos/{genome_build}/ucsc_conservation/download/knownGene.exonAA.fa.gz",
    output:
        vcf="annos/{genome_build}/ucsc_conservation/ucsc_conservation.vcf.gz",
        tbi="annos/{genome_build}/ucsc_conservation/ucsc_conservation.vcf.gz.tbi",
    shell:
        r"""
        python scripts/knowngeneaa.py \
            {input.hgnc} \
            {input.enst_ensg} \
            {input.reference} \
            {input.fa} \
            --output /dev/stdout \
        | bcftools sort \
            -O z \
            -o {output.vcf}
        tabix -f {output.vcf}
        """


rule annos_ucsc_conservation_to_tsv:
    input:
        header="header/knowngeneaa.txt",
        vcf="annos/{genome_build}/ucsc_conservation/ucsc_conservation.vcf.gz",
    output:
        tsv="annos/{genome_build}/ucsc_conservation/ucsc_conservation.tsv",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools query \
                -f "%CHROM\t%POS\t%END\t%HGNC_ID\t%ENST_ID\t%EXON\t%EXON_COUNT\t%ALIGNMENT\n" \
                {input.vcf} \
            | uniq
        ) \
        > {output.tsv}
        """


rule annos_dbsnp_download:
    output:
        vcf="annos/{genome_build}/dbsnp/dbsnp.vcf.gz",
        vcf_tbi="annos/{genome_build}/dbsnp/dbsnp.vcf.gz.tbi",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
        tabix -f {output.vcf}
        """


CADD_VERSION = "1.6"
CADD_PREFIX = f"https://kircherlab.bihealth.org/download/CADD/v{CADD_VERSION}"


rule annos_cadd_download:
    output:
        tsv="annos/{genome_release}/cadd/download/{filename}.tsv.gz",
        tsv_tbi="annos/{genome_release}/cadd/download/{filename}.tsv.gz.tbi",
    shell:
        r"""
        for path in {output};
        do
            aria2c \
                --check-certificate=false \
                --file-allocation=trunc \
                --out=$path \
                --split=16 \
                --max-concurrent-downloads=16 \
                --max-connection-per-server=16 \
                {CADD_PREFIX}/$(echo {wildcards.genome_release} | sed -e 's/grch/GRCh/')/$(basename $path)
        done
        """


rule annos_cadd_process_37:
    input:
        "annos/grch37/cadd/download/whole_genome_SNVs_inclAnno.tsv.gz",
        "annos/grch37/cadd/download/whole_genome_SNVs_inclAnno.tsv.gz.tbi",
        "annos/grch37/cadd/download/InDels_inclAnno.tsv.gz",
        "annos/grch37/cadd/download/InDels_inclAnno.tsv.gz.tbi",
    output:
        touch("annos/grch37/cadd/.done"),


rule annos_cadd_process_38:
    input:
        "annos/grch38/cadd/download/whole_genome_SNVs_inclAnno.tsv.gz",
        "annos/grch38/cadd/download/whole_genome_SNVs_inclAnno.tsv.gz.tbi",
        "annos/grch38/cadd/download/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz",
        "annos/grch38/cadd/download/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz.tbi",
    output:
        touch("annos/grch38/cadd/.done"),


DBNSFP_VERSION = "4.4"
DBNSFP_ACADEMIC_URL = "https://usf.box.com/shared/static/bvfzmkpgtphvbmmrvb2iyl2jl21o49kc"
DBNSFP_COMMMERCIAL_URL = "https://usf.box.com/shared/static/a84zcdlkx2asq2nxh6xr2gdb4csmyvhk"


def files_dbnsfp():
    lst = [
        "dbNSFP{version}{variant}.readme.txt",
        "dbNSFP{version}{variant}_variant.chr1.gz",
        "dbNSFP{version}{variant}_variant.chr10.gz",
        "dbNSFP{version}{variant}_variant.chr11.gz",
        "dbNSFP{version}{variant}_variant.chr12.gz",
        "dbNSFP{version}{variant}_variant.chr13.gz",
        "dbNSFP{version}{variant}_variant.chr14.gz",
        "dbNSFP{version}{variant}_variant.chr15.gz",
        "dbNSFP{version}{variant}_variant.chr16.gz",
        "dbNSFP{version}{variant}_variant.chr17.gz",
        "dbNSFP{version}{variant}_variant.chr18.gz",
        "dbNSFP{version}{variant}_variant.chr19.gz",
        "dbNSFP{version}{variant}_variant.chr2.gz",
        "dbNSFP{version}{variant}_variant.chr20.gz",
        "dbNSFP{version}{variant}_variant.chr21.gz",
        "dbNSFP{version}{variant}_variant.chr22.gz",
        "dbNSFP{version}{variant}_variant.chr3.gz",
        "dbNSFP{version}{variant}_variant.chr4.gz",
        "dbNSFP{version}{variant}_variant.chr5.gz",
        "dbNSFP{version}{variant}_variant.chr6.gz",
        "dbNSFP{version}{variant}_variant.chr7.gz",
        "dbNSFP{version}{variant}_variant.chr8.gz",
        "dbNSFP{version}{variant}_variant.chr9.gz",
        "dbNSFP{version}{variant}_variant.chrM.gz",
        "dbNSFP{version}{variant}_variant.chrX.gz",
        "dbNSFP{version}{variant}_variant.chrY.gz",
        "dbNSFP{version}_gene.complete.gz",
        "dbNSFP{version}_gene.gz",
        "LICENSE.txt",
        "try.vcf",
        "tryhg18.in",
        "tryhg19.in",
        "tryhg38.in",
    ]
    return ["annos/grch37/dbnsfp-{version}{variant}/download/%s" % e for e in lst]


rule annos_dbnsfp_download:
    output:
        files_dbnsfp(),
        zip="annos/grch37/dbnsfp-{version}{variant}/download/dbNSFP{version}{variant}.zip",
    wildcard_constraints:
        version=r"\d\.\d",
    shell:
        r"""
        if [[ "{wildcards.variant}" == a ]]; then
            url={DBNSFP_ACADEMIC_URL}
        else
            url={DBNSFP_COMMMERCIAL_URL}
        fi

        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.zip} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            $url
        unzip -d $(dirname {output.zip}) {output.zip}
        """


rule annos_dbnsfp_process:
    input:
        zip="annos/grch37/dbnsfp-{version}{variant}/download/dbNSFP{version}{variant}.zip",
    wildcard_constraints:
        version=r"\d\.\d",
    output:
        touch("annos/{genome_release}/dbnsfp-{version}{variant}/.done"),


def files_dbscsnv(version: str = "1.1"):
    """Files contained in the dbscSNV ZIP file."""
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    return [f"annos/grch37/dbscsnv/download/dbscSNV{version}.chr{chrom}" for chrom in chroms]


rule annos_dbscsnv_download:
    output:
        files_dbscsnv(),
        zip="annos/grch37/dbscsnv/download/dbscSNV1.1.zip",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.zip} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip
        unzip -d $(dirname {output.zip}) {output.zip}
        """


rule annos_dbscsnv_process:
    input:
        zip="annos/grch37/dbscsnv/download/dbscSNV1.1.zip",
    output:
        touch("annos/{genome_release}/dbscsnv/.done"),
