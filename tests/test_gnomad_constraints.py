import ast
import gzip
import pathlib
import re
import tempfile
import unittest

import polars as pl

from varfish_db_downloader import wget

RULES = pathlib.Path("rules/work/genes/gnomad.smk")
VERSIONS = pathlib.Path("varfish_db_downloader/versions.py")
CONVERTER = pathlib.Path("rules/work/genes/scripts/gnomad_constraints_v4_to_tsv.py")


def pinned_version(name: str) -> str:
    """Read a keyword argument of ``DataVersions(...)`` without importing ``versions.py``."""
    tree = ast.parse(VERSIONS.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "DataVersions":
            for keyword in node.keywords:
                if keyword.arg == name:
                    return ast.literal_eval(keyword.value)
    raise KeyError(name)


def constraints_download_url() -> str:
    """Return the URL that ``genes_gnomad_download_v4`` fetches for the pinned version."""
    rule = RULES.read_text().split("rule genes_gnomad_download_v4:")[1].split("\nrule ")[0]
    template = re.search(r"https://\S+", rule).group(0)
    return template.replace("{wildcards.version}", pinned_version("gnomad_constraints"))


class GnomadConstraintsDownloadTest(unittest.TestCase):
    def test_download_url_has_test_mode_excerpt(self):
        url = constraints_download_url()
        entries = {entry.url: entry for entry in wget.load_urls_yaml("download_urls.yml")}

        self.assertIn(url, entries)
        excerpt = pathlib.Path("excerpt-data") / entries[url].hash / url.split("/")[-1]
        self.assertTrue(excerpt.exists(), excerpt)


class GnomadConstraintsConverterTest(unittest.TestCase):
    def test_converts_the_test_mode_excerpt(self):
        url = constraints_download_url()
        entry = {e.url: e for e in wget.load_urls_yaml("download_urls.yml")}[url]
        excerpt = pathlib.Path("excerpt-data") / entry.hash / url.split("/")[-1]

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = pathlib.Path(tmp_dir)
            constraints_tsv = tmp / "constraints.tsv"
            constraints_tsv.write_bytes(gzip.decompress(excerpt.read_bytes()))
            # A1BG: MANE/canonical ENST00000263100 (see the excerpt).
            xlink = tmp / "xlink.tsv"
            xlink.write_text(
                "ensembl_gene_id\tensembl_transcript_id\tentrez_id\tgene_symbol\n"
                "ENSG00000121410\tENST00000263100\t1\tA1BG\n"
            )
            output = tmp / "out.tsv"

            namespace = {"__name__": "converter"}
            exec(compile(CONVERTER.read_text(), str(CONVERTER), "exec"), namespace)
            namespace["main"](str(constraints_tsv), str(xlink), str(output), str(tmp / "out.md5"))

            df = pl.read_csv(output, separator="\t", infer_schema=False)
            self.assertEqual(df.get_column("ensembl_gene_id").to_list(), ["ENSG00000121410"])
            row = df.row(0, named=True)
            self.assertEqual(row["gene_symbol"], "A1BG")
            self.assertNotEqual(row["oe_lof_upper"], "NA")
            self.assertNotEqual(row["pLI"], "NA")


if __name__ == "__main__":
    unittest.main()
