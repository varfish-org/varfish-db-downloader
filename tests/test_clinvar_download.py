import ast
import pathlib
import re
import unittest

from varfish_db_downloader import wget

RULES = pathlib.Path("rules/work/misc/clinvar_download.smk")
VERSIONS = pathlib.Path("varfish_db_downloader/versions.py")


def module_constant(name: str) -> str:
    """Read a module-level string constant of ``versions.py`` without importing it."""
    for node in ast.parse(VERSIONS.read_text()).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == name:
            return ast.literal_eval(node.value)
    raise KeyError(name)


def clinvar_download_url() -> str:
    """Return the URL that ``clinvar_download.smk`` fetches for the pinned ClinVar release."""
    template = re.search(r"https://\S+", RULES.read_text()).group(0)
    return template.replace(
        "{wildcards.clinvar_version}", module_constant("CLINVAR_VERSION")
    ).replace("{wildcards.annonars_version}", module_constant("CLINVAR_THIS"))


class ClinvarDownloadTest(unittest.TestCase):
    def test_download_url_has_test_mode_excerpt(self):
        url = clinvar_download_url()
        entries = {entry.url: entry for entry in wget.load_urls_yaml("download_urls.yml")}

        self.assertIn(url, entries)
        excerpt = pathlib.Path("excerpt-data") / entries[url].hash / url.split("/")[-1]
        self.assertTrue(excerpt.exists(), excerpt)


if __name__ == "__main__":
    unittest.main()
