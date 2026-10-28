import pathlib
import re
import unittest

#: Markers of downloads from NCBI, which answers concurrent requests with HTTP 503.
NCBI_PATTERN = re.compile(r"ftp(-trace)?\.ncbi\.(nlm\.)?nih\.gov|refseq_base_url")
#: Top-level blocks of a Snakemake file (rules, functions, assignments).
BLOCK_PATTERN = re.compile(r"^\S.*?(?=^\S|\Z)", re.S | re.M)


def ncbi_download_rules() -> dict[str, str]:
    """Return the rules that download from NCBI, mapped to their source."""
    result = {}
    for path in sorted(pathlib.Path("rules").rglob("*.smk")):
        blocks = BLOCK_PATTERN.findall(path.read_text())
        ncbi_constants = [
            m.group(1)
            for block in blocks
            if (m := re.match(r"([A-Z_]+)\s*=", block)) and NCBI_PATTERN.search(block)
        ]
        for block in blocks:
            m = re.match(r"rule (\w+):", block)
            uses_constant = any(re.search(rf"\b{name}\b", block) for name in ncbi_constants)
            if m and (NCBI_PATTERN.search(block) or uses_constant):
                result[m.group(1)] = block
    return result


class NcbiDownloadsTest(unittest.TestCase):
    def test_ncbi_download_rules_use_the_ncbi_downloads_resource(self):
        rules = ncbi_download_rules()

        self.assertIn("reference_download", rules)
        self.assertIn("annos_dbsnp_download", rules)
        missing = sorted(name for name, block in rules.items() if "ncbi_downloads=1" not in block)
        self.assertEqual(missing, [])

    def test_slurm_run_allows_one_ncbi_download_at_a_time(self):
        self.assertIn("--resources ncbi_downloads=1", pathlib.Path("run-slurm.sh").read_text())


if __name__ == "__main__":
    unittest.main()
