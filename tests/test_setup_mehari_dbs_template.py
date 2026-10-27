import pathlib
import subprocess
import types
import unittest

from jinja2 import Template

TEMPLATE = pathlib.Path("rules/output/mehari/setup_mehari_dbs.sh.j2")


def render(dv: types.SimpleNamespace, pv: types.SimpleNamespace) -> str:
    return Template(TEMPLATE.read_text()).render(
        DV=dv,
        PV=pv,
        genomebuilds=["grch37"],
        mehari_tx_files={},
        mehari_freq_versions={},
        downloader_base="/data",
    )


def clinvar_file_name(script: str, assembly: str) -> str:
    """Evaluate the ClinVar file name as the rendered script does."""
    lines = script.splitlines()
    assignments = [line for line in lines if line[:1].isupper() and "=" in line.split()[0]]
    file_line = next(
        line.strip() for line in lines if line.strip().startswith('file="annonars-clinvar-minimal')
    )
    return subprocess.run(
        [
            "bash",
            "-c",
            "\n".join([*assignments, f"assembly={assembly}", file_line, 'echo "$file"']),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


class SetupMehariDbsTemplateTest(unittest.TestCase):
    def test_clinvar_file_name_uses_annonars_version_of_clinvar_data(self):
        dv = types.SimpleNamespace(
            clinvar_version="20260226",
            annonars_clinvar="0.44.0",
            hgnc_quarterly="2026-10-06",
            mehari_tx="0.14.1",
        )
        pv = types.SimpleNamespace(annonars="0.45.1")

        script = render(dv, pv)

        self.assertEqual(
            clinvar_file_name(script, "grch37"),
            "annonars-clinvar-minimal-grch37-20260226+0.44.0",
        )


if __name__ == "__main__":
    unittest.main()
