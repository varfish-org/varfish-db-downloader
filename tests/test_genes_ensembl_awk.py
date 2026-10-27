import pathlib
import subprocess
import tempfile
import unittest

#: ENSEMBL entrez.tsv excerpt; ENST00000000002 has two Entrez IDs, the trans_name row is ignored.
ENTREZ_TSV = """\
gene_stable_id	transcript_stable_id	protein_stable_id	xref	db_name	info_type	source_identity	xref_identity	linkage_type
ENSG00000000001	ENST00000000001	ENSP00000000001	1001	EntrezGene	DEPENDENT	-	-	-
ENSG00000000001	ENST00000000002	-	1001	EntrezGene	DEPENDENT	-	-	-
ENSG00000000001	ENST00000000002	-	1002	EntrezGene	DEPENDENT	-	-	-
ENSG00000000001	ENST00000000001	-	GENE1-201	EntrezGene_trans_name	DIRECT	-	-	-
"""

#: ENSEMBL GTF excerpt; ENSG00000000002 has no gene_name.
GTF = "\n".join(
    "\t".join(fields)
    for fields in [
        ["1", "ensembl", "gene", "100", "900", ".", "+", ".", 'gene_id "ENSG00000000001"; gene_version "1"; gene_name "GENE1"; gene_biotype "protein_coding";'],
        ["1", "ensembl", "transcript", "100", "900", ".", "+", ".", 'gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; gene_name "GENE1";'],
        ["1", "ensembl", "exon", "100", "200", ".", "+", ".", 'gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; gene_name "GENE1";'],
        ["1", "ensembl", "transcript", "150", "900", ".", "+", ".", 'gene_id "ENSG00000000001"; transcript_id "ENST00000000002"; gene_name "GENE1";'],
        ["1", "ensembl", "transcript", "300", "800", ".", "+", ".", 'gene_id "ENSG00000000001"; transcript_id "ENST00000000003"; gene_name "GENE1";'],
        ["1", "ensembl", "gene", "1000", "1900", ".", "-", ".", 'gene_id "ENSG00000000002"; gene_version "1"; gene_biotype "lncRNA";'],
        ["1", "ensembl", "transcript", "1000", "1900", ".", "-", ".", 'gene_id "ENSG00000000002"; transcript_id "ENST00000000004";'],
    ]
) + "\n"  # fmt: skip


class GenesEnsemblAwkTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        tmp = pathlib.Path(self.tmp_dir.name)
        self.entrez_tsv = tmp / "entrez.tsv"
        self.entrez_tsv.write_text(ENTREZ_TSV)
        self.gtf = tmp / "genes.gtf"
        self.gtf.write_text(GTF)

    def tearDown(self):
        self.tmp_dir.cleanup()

    def run_awk(self, script: str, *inputs: pathlib.Path) -> str:
        return subprocess.run(
            ["awk", "-F", "\t", "-f", script, *map(str, inputs)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout

    def test_xlink_has_one_line_per_transcript_and_entrez_id(self):
        output = self.run_awk("scripts/genes-ensembl-xlink.awk", self.entrez_tsv, self.gtf)

        self.assertEqual(
            output,
            "ENSG00000000001\tENST00000000001\t1001\tGENE1\n"
            "ENSG00000000001\tENST00000000002\t1001\tGENE1\n"
            "ENSG00000000001\tENST00000000002\t1002\tGENE1\n"
            "ENSG00000000001\tENST00000000003\t\tGENE1\n"
            "ENSG00000000002\tENST00000000004\t\t\n",
        )

    def test_genesymbol_skips_genes_without_symbol(self):
        output = self.run_awk("scripts/genes-ensembl-genesymbol.awk", self.gtf)

        self.assertEqual(output, "ENSG00000000001\tGENE1\n")


if __name__ == "__main__":
    unittest.main()
