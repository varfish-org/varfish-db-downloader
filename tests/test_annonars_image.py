import os
import pathlib
import subprocess
import tempfile
import unittest

from click.testing import CliRunner

from varfish_db_downloader.cli import main
from varfish_db_downloader.rocksdb import max_table_format_version

WRAPPER = pathlib.Path("annonars-image-bin/annonars")


class AnnonarsImageWrapperTest(unittest.TestCase):
    def test_runs_annonars_from_pinned_docker_image(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake_apptainer = pathlib.Path(tmp) / "apptainer"
            fake_apptainer.write_text('#!/usr/bin/bash\nprintf "%s\\n" "$@"\n')
            fake_apptainer.chmod(0o755)
            env = {**os.environ, "PATH": f"{tmp}:{os.environ['PATH']}"}

            result = subprocess.run(
                [str(WRAPPER), "genes", "import", "--path-out-rocksdb", "out dir"],
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )

        self.assertEqual(
            result.stdout.splitlines(),
            [
                "exec",
                "docker://ghcr.io/varfish-org/annonars:0.45.2",
                "annonars",
                "genes",
                "import",
                "--path-out-rocksdb",
                "out dir",
            ],
        )


def write_options(rocksdb_dir: pathlib.Path, name: str, format_version: int):
    """Write a RocksDB OPTIONS file with one column family in ``rocksdb_dir``."""
    rocksdb_dir.mkdir(parents=True, exist_ok=True)
    (rocksdb_dir / name).write_text(
        "[Version]\n"
        "  rocksdb_version=11.8.1\n"
        '[TableOptions/BlockBasedTable "default"]\n'
        f"  format_version={format_version}\n"
    )


class TableFormatVersionTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.rocksdb_dir = pathlib.Path(self.tmp_dir.name) / "genes" / "rocksdb"

    def tearDown(self):
        self.tmp_dir.cleanup()

    def test_max_table_format_version_over_all_options_files(self):
        write_options(self.rocksdb_dir, "OPTIONS-000005", 5)
        write_options(self.rocksdb_dir, "OPTIONS-000007", 7)

        self.assertEqual(max_table_format_version(self.rocksdb_dir), 7)

    def test_check_format_fails_above_max_format_version(self):
        write_options(self.rocksdb_dir, "OPTIONS-000005", 7)

        result = CliRunner().invoke(main, ["rocksdb", "check-format", str(self.rocksdb_dir)])

        self.assertEqual(result.exit_code, 1)
        self.assertIn(f"{self.rocksdb_dir}: format_version 7", result.output)

    def test_check_format_passes_readable_format_version(self):
        write_options(self.rocksdb_dir, "OPTIONS-000005", 5)

        result = CliRunner().invoke(main, ["rocksdb", "check-format", str(self.rocksdb_dir)])

        self.assertEqual(result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()
