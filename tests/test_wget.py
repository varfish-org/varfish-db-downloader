import pathlib
import tempfile
import unittest
from unittest import mock

from varfish_db_downloader.wget import UrlEntry, copy_excerpt

#: A BioMart query URL, whose excerpt is stored under the name from the URL path.
BIOMART_URL = (
    "https://grch37.ensembl.org/biomart/martservice?query=<Query><Dataset "
    'name = "hsapiens_gene_ensembl" /></Query>'
)


class CopyExcerptTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.data_dir = pathlib.Path(self.tmp_dir.name) / "excerpt-data"
        self.output = pathlib.Path(self.tmp_dir.name) / "output.tsv"

    def tearDown(self):
        self.tmp_dir.cleanup()

    def write_excerpt(self, url: UrlEntry, name: str, content: str):
        excerpt_dir = self.data_dir / url.hash
        excerpt_dir.mkdir(parents=True)
        (excerpt_dir / name).write_text(content)

    def test_copies_excerpt_named_after_url_path_when_server_unreachable(self):
        url = UrlEntry(url=BIOMART_URL)
        self.write_excerpt(url, "martservice", "ENSG00000261657\tSLC25A26\n")

        with mock.patch("varfish_db_downloader.wget.get_filename_from_url", return_value=None):
            copy_excerpt(url, str(self.data_dir), str(self.output))

        self.assertEqual(self.output.read_text(), "ENSG00000261657\tSLC25A26\n")

    def test_copies_excerpt_named_by_server_when_url_path_has_no_excerpt(self):
        url = UrlEntry(url="https://usf.box.com/shared/static/2hzcx5s6p1xui7oen16xqzndfrkt8l9l")
        self.write_excerpt(url, "dbNSFP4.5a.zip", "zip content")

        with mock.patch(
            "varfish_db_downloader.wget.get_filename_from_url", return_value="dbNSFP4.5a.zip"
        ):
            copy_excerpt(url, str(self.data_dir), str(self.output))

        self.assertEqual(self.output.read_text(), "zip content")


if __name__ == "__main__":
    unittest.main()
