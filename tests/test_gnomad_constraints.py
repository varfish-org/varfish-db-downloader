import ast
import pathlib
import re
import tempfile
import unittest

import polars as pl

from varfish_db_downloader import wget

RULES = pathlib.Path("rules/work/genes/gnomad.smk")
VERSIONS = pathlib.Path("varfish_db_downloader/versions.py")
CONVERTER = pathlib.Path("rules/work/genes/scripts/gnomad_constraints_v4_to_tsv.py")

#: The A1BG rows of gnomAD v4.1.2 ``constraint_metrics.tsv.bgz``, reduced to the columns the
#: converter reads (inline, because CI's lint job checks out without Git LFS).
CONSTRAINTS_TSV = "\n".join(
    "\t".join(fields)
    for fields in [
        ["gene", "gene_id", "transcript", "canonical", "mane_select", "cds_length", "lof_hc_lc.exp", "mis.exp", "syn.exp", "mis.z_score", "lof_hc_lc.obs", "mis.obs", "syn.obs", "lof.oe", "lof.oe_ci.lower", "lof.oe_ci.upper", "mis.oe", "mis.oe_ci.lower", "mis.oe_ci.upper", "syn.oe", "syn.oe_ci.lower", "syn.oe_ci.upper", "lof.pLI", "syn.z_score"],
        ["A1BG", "1", "NM_130786.4", "true", "true", "NA", "4.6940e+01", "6.8884e+02", "3.1470e+02", "-2.9516e-01", "45", "707", "316", "9.5603e-01", "7.5098e-01", "1.2292e+00", "1.0264e+00", "9.6512e-01", "1.0922e+00", "1.0041e+00", "9.1610e-01", "1.1021e+00", "6.7352e-15", "-4.3795e-02"],
        ["A1BG", "ENSG00000121410", "ENST00000263100", "true", "true", "1485", "4.6940e+01", "6.8884e+02", "3.1470e+02", "-2.9516e-01", "45", "707", "316", "9.5603e-01", "7.5098e-01", "1.2292e+00", "1.0264e+00", "9.6512e-01", "1.0922e+00", "1.0041e+00", "9.1610e-01", "1.1021e+00", "6.7352e-15", "-4.3795e-02"],
        ["A1BG", "ENSG00000121410", "ENST00000600966", "false", "false", "917", "2.8743e+01", "4.1368e+02", "1.8031e+02", "3.0779e-01", "24", "399", "166", "7.7585e-01", "5.5025e-01", "1.1172e+00", "9.6452e-01", "8.8881e-01", "1.0478e+00", "9.2063e-01", "8.1153e-01", "1.0471e+00", "2.1620e-06", "6.3588e-01"],
    ]
) + "\n"  # fmt: skip


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
    def test_keeps_the_mane_transcript_of_v4_1_2(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = pathlib.Path(tmp_dir)
            constraints_tsv = tmp / "constraints.tsv"
            constraints_tsv.write_text(CONSTRAINTS_TSV)
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
            self.assertEqual(row["oe_lof_upper"], "1.2292")
            self.assertEqual(row["pLI"], "6.7352e-15")


if __name__ == "__main__":
    unittest.main()
