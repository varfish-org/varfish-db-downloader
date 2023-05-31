rule features_grch37_tad_domains:
    output:
        download_imr90="features/grch37/tads/download/IMR90_domains_hg19.bed",
        download_hesc="features/grch37/tads/download/hESC_domains_hg19.bed",
        bed_imr90="features/grch37/tads/imr90.bed",
        bed_imr90_md5="features/grch37/tads/imr90.bed.md5",
        bed_hesc="features/grch37/tads/hesc.bed",
        bed_hesc_md5="features/grch37/tads/hesc.bed.md5",
    shell:
        r"""
        set -x

        wget --no-check-certificate \
            -O {output.download_imr90} \
            http://compbio.med.harvard.edu/modencode/webpage/hic/IMR90_domains_hg19.bed
        wget --no-check-certificate \
            -O {output.download_hesc} \
            http://compbio.med.harvard.edu/modencode/webpage/hic/hESC_domains_hg19.bed

        echo -e "#chrom\tbegin\tend" >{output.bed_imr90}
        sed -e 's/^chr//' {output.download_imr90} >>{output.bed_imr90}

        echo -e "#chrom\tbegin\tend" >{output.bed_hesc}
        sed -e 's/^chr//' {output.download_hesc} >>{output.bed_hesc}

        md5sum {output.bed_imr90} >{output.bed_imr90_md5}
        md5sum {output.bed_hesc} >{output.bed_hesc_md5}
        """


rule features_grch37_refseq_gene_regions:
    output:
        download_acc="features/grch37/gene_regions/download/chr_accessions_GRCh37.p13",
        download_gtf="features/grch37/gene_regions/download/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
        tsv="features/grch37/gene_regions/refseq.bed.gz",
        tsv_md5="features/grch37/gene_regions/refseq.bed.gz.md5",
        tsv_tbi="features/grch37/gene_regions/refseq.bed.gz.tbi",
        tsv_tbi_md5="features/grch37/gene_regions/refseq.bed.gz.tbi.md5",
    shell:
        r"""
        set -x

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        wget --no-check-certificate \
            -O {output.download_acc} \
            'http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/Assembled_chromosomes/chr_accessions_GRCh37.p13'

        wget --no-check-certificate \
            -O {output.download_gtf} \
            'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz'

        awk \
            -F $'\t' \
            -f scripts/features-refseq-gene-regions.awk \
            {output.download_acc} \
            <(zcat {output.download_gtf}) \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """


rule features_grch37_ensembl_gene_regions:
    output:
        download_gtf="features/grch37/gene_regions/download/Homo_sapiens.GRCh37.87.gtf.gz",
        tsv="features/grch37/gene_regions/ensembl.bed.gz",
        tsv_md5="features/grch37/gene_regions/ensembl.bed.gz.md5",
        tsv_tbi="features/grch37/gene_regions/ensembl.bed.gz.tbi",
        tsv_tbi_md5="features/grch37/gene_regions/ensembl.bed.gz.tbi.md5",
    shell:
        r"""
        set -x

        wget --no-check-certificate \
            -O {output.download_gtf} \
            'http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'

        awk \
            -F $'\t' \
            -f scripts/features-ensembl-gene-regions.awk \
            <(zcat {output.download_gtf}) \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | awk -F $'\t' \
            -f scripts/sort-bed.awk \
        | bgzip -c \
        > {output.tsv}

        tabix -f {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        md5sum {output.tsv_tbi} >{output.tsv_tbi_md5}
        """


rule features_grch37_masked_repeat:
    input:
        bed="tracks/grch37/ucsc_rmsk.bed.gz",
    output:
        bed="features/grch37/masked/repeat.bed.gz",
        bed_md5="features/grch37/masked/repeat.bed.gz.md5",
        bed_tbi="features/grch37/masked/repeat.bed.gz.tbi",
        bed_tbi_md5="features/grch37/masked/repeat.bed.gz.tbi.md5",
    shell:
        r"""
        set -x

        zcat {input.bed} \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | egrep -v '^Un|_random|_fix|_alt|_hap' \
        | bgzip -c \
        > {output.bed}
        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule features_grch37_masked_segdup:
    input:
        bed="tracks/grch37/ucsc_genomicSuperDups.bed.gz",
    output:
        bed="features/grch37/masked/segdup.bed.gz",
        bed_md5="features/grch37/masked/segdup.bed.gz.md5",
        bed_tbi="features/grch37/masked/segdup.bed.gz.tbi",
        bed_tbi_md5="features/grch37/masked/segdup.bed.gz.tbi.md5",
    shell:
        r"""
        set -x

        zcat {input.bed} \
        | egrep '^#|^X|^Y|^M|^[1-9]' \
        | egrep -v '^Un|_random|_fix|_alt|_hap' \
        | bgzip -c \
        > {output.bed}
        tabix -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
