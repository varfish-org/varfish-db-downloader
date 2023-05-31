import io
import unittest
from unittest import mock as mock

from ..ucsc_binning import main

expected = """\
release	chromosome	start	end	bin	reference	alternative	value
GRCh37	1	10	1000000	73	A	C	something
GRCh37	1	100	100	585	A	C	something
GRCh37	1	200	200	585	A	C	something
GRCh37	1	1000000	1000000	592	A	C	something
"""


class TestUcscBinning(unittest.TestCase):
    @mock.patch("sys.stdout", new_callable=io.StringIO)
    def test_main_valid_input(self, mock_stdout):
        main(["--input", "tools/tests/input_valid.tsv"])
        self.assertEqual(mock_stdout.getvalue(), expected)

    @mock.patch("sys.stdout", new_callable=io.StringIO)
    def test_main_invalid_input_bin_column_missing(self, mock_stdout):
        with self.assertRaises(Exception):
            main(["--input", "tools/tests/input_invalid.tsv"])
