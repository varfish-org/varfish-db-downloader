## Rules related to UCSC conservation track.


rule annos_features_cons_download:
    output:
        fa="work/download/annos/{genome_release}/features/cons/{version}/knownGene.exonAA.fa.gz",
    shell:
        r"""
        if [[ {wildcards.genome_release} == grch37 ]]; then
            ucsc_name=hg19
        else
            ucsc_name=hg38
        fi

        # Check version.
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        wget -O $TMPDIR/alignment.html \
            https://hgdownload.cse.ucsc.edu/goldenpath/${{ucsc_name}}/multiz100way/alignments
        version=$(grep knownGene.exonAA.fa.gz $TMPDIR/listing.html | awk '{{ print $$2 }}')
        if [[ "$version" != "{wildcards.version}" ]]; then
            >&2 echo "Version mismatch for knownGene.exonAA.fa.gz: expected {version}, got $version"
            exit 1
        fi


        # Actually download the file.
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.fa} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            "https://hgdownload.cse.ucsc.edu/goldenpath/${{ucsc_name}}/multiz100way/alignments/knownGene.exonAA.fa.gz"
        """


rule annos_features_cons_to_vcf:
    input:
        hgnc="work/genes/hgnc/hgnc_info.jsonl",
        enst_ensg="work/genes/enst_ensg/{genome_release}/enst_ensg.tsv",
        reference="work/reference/{genome_release}/reference.fa",
        fa="work/download/annos/{genome_release}/features/cons/knownGene.exonAA.fa.gz",
    output:
        vcf="work/download/annos/{genome_release}/features/cons/ucsc_conservation.vcf.gz",
        tbi="work/download/annos/{genome_release}/features/cons/ucsc_conservation.vcf.gz.tbi",
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


rule annos_features_cons_to_tsv:
    input:
        vcf="work/download/annos/{genome_release}/features/cons/{version}/ucsc_conservation.vcf.gz",
    output:
        tsv="work/annos/{genome_release}/features/cons/{version}/ucsc_conservation.tsv",
        tsv_md5="work/annos/{genome_release}/features/cons/{version}/ucsc_conservation.tsv.md5",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        echo -e "chromosome\tstart\tstop\thgnc_id\tenst_id\texon_num\texon_count\talignment" >$TMPDIR/header.txt

        cat $TMPDIR/header.txt \
        | grep -v '^$' \
        | tr '\n' '\t' \
        | sed -e 's/\t*$/\n/g' \
        > {output.tsv}

        bcftools query \
            -f "%CHROM\t%POS\t%END\t%HGNC_ID\t%ENST_ID\t%EXON\t%EXON_COUNT\t%ALIGNMENT\n" \
            {input.vcf} \
        | uniq
        >> {output.tsv}

        md5sum {output.tsv} > {output.tsv_md5}
        """
