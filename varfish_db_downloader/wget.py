"""Implementation of the stub ``wget`` and supporting tool."""

import hashlib
import os
import pathlib
import shutil
import subprocess
import tempfile
import typing
import urllib.parse

import attrs
import cattrs
import requests
import yaml
from loguru import logger


def excerpt_manual(url: str, path_out: str, count: int):
    """Do not download, assume output is already there."""
    _ = url
    _ = count
    if not os.path.exists(path_out):
        raise RuntimeError(f"File {path_out} does not exist")


def no_excerpt(url: str, path_out: str, count: int):
    """Do not excerpt, use all."""
    _ = count
    with open(path_out, 'wb') as outputf:
        r = requests.get(url, allow_redirects=True)
        outputf.write(r.content)


def excerpt_head(url: str, path_out: str, count: int):
    """Excerpt a plaint-text file by copying lines."""
    with open(path_out, "wb") as f_out:
        r = requests.get(url, stream=True)
        for i, line in enumerate(r.iter_lines()):
            if i >= count:
                break
            f_out.write(line)


def excerpt_vcf_head(url: str, path_out: str, count: int):
    """Excerpt a VCF file by applying the following steps:

    First, the chromosomes are listed.

    Then, the following steps will be performed:

    1. extract the VCF header with ``tabix`` to a temporary file
    2. append the first ``count`` non-header lines to the temporary file
    3. compress the temporary file with ``bgzip`` to the output file
    4. index the output file with ``tabix``
    """
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        # Obtain the list of chromosomes.
        lst_cmd = f"tabix --list-chroms {url}"
        logger.debug("    + {}", lst_cmd)
        chroms_out = subprocess.check_output(lst_cmd, shell=True, cwd=tmpdir)
        chroms_lst = [l.strip() for l in chroms_out.decode("utf-8").split("\n") if l.strip()]
        chroms_arg = " ".join(chroms_lst)
        # Perform the extraction steps.
        tmpfile = tmpdir + "/tmp.vcf"
        cmds = [
            f"tabix --only-header {url} > {tmpfile}",
            f"tabix {url} {chroms_arg} | head -n {count} >> {tmpfile}",
            f"bgzip -c {tmpfile} > {cwd}/{path_out}",
            f"tabix --preset vcf {cwd}/{path_out}",
        ]
        for cmd in cmds:
            logger.debug("    + {}", cmd)
            subprocess.check_call(cmd, shell=True, cwd=tmpdir)


#: Length of the hash to use.  The hash is created from the underlying URL.
HASH_LENGTH = 16

#: Mapping from strategy name to implementation.
STRATEGY_MAP = {
    "head": excerpt_head,
    "gz-head": excerpt_head,
    "vcf-head": excerpt_vcf_head,
    "manual": excerpt_manual,
    "no-excerpt": no_excerpt,
}


@attrs.frozen()
class ExcerptStrategy:
    """Defines a strategy to use for excerpting a file."""

    #: The strategy to use.
    strategy: str
    #: The number of records to excerpt.
    count: typing.Optional[int]

    @staticmethod
    def default_for(url: str) -> "ExcerptStrategy":
        """Return the default strategy."""
        if url.endswith(".vcf.gz") or url.endswith(".vcf.bgz"):
            return ExcerptStrategy(strategy="vcf-head", count=100)
        elif url.endswith(".gz") or url.endswith(".bgz"):
            return ExcerptStrategy(strategy="gz-head", count=100)
        else:
            return ExcerptStrategy(strategy="head", count=100)


@attrs.frozen()
class UrlEntry:
    """Entry in the ``downloads_urls.yml`` file."""

    #: The URL of the original file download.
    url: str
    #: The hash to use (computed from ``url``).
    hash: typing.Optional[str] = None
    #: The excerpt strategy to use.
    excerpt_strategy: typing.Optional[ExcerptStrategy] = None

    def __attrs_post_init__(self):
        object.__setattr__(
            self, "hash", hashlib.sha256(self.url.encode("utf-8")).hexdigest()[:HASH_LENGTH]
        )


def load_urls_yaml(yaml_path: str) -> typing.List[UrlEntry]:
    """Load the URLs YAML file."""
    with open(yaml_path, "rt") as f:
        result = []
        raw = yaml.safe_load(f)
        for url in cattrs.structure(raw, typing.List[UrlEntry]):
            if url.excerpt_strategy is None:
                url = attrs.evolve(url, excerpt_strategy=ExcerptStrategy.default_for(url.url))
            result.append(url)
        return result


def download_excerpt(url: UrlEntry, data_dir: str, force: bool):
    """Download excerpt for the given ``url``."""
    logger.info("  Downloading excerpt for {}...", url.url)
    logger.info("    hash is {}", url.hash)
    logger.info("    strategy is {}", url.excerpt_strategy.strategy)

    out_path = pathlib.Path(data_dir) / url.hash
    if out_path.exists() and not force and url.excerpt_strategy.strategy != "manual":
        logger.warning("  ...skipping, already exists and not --force")
        return

    out_path.mkdir(parents=True, exist_ok=True)
    out_path_url = out_path / "url.txt"
    logger.info("    writing URL to {}", out_path_url)
    with out_path_url.open("wt") as f:
        print(url.url, file=f)

    excerpt_fun = STRATEGY_MAP[url.excerpt_strategy.strategy]
    parsed = urllib.parse.urlparse(url.url)
    basename = parsed.path.split("/")[-1]
    out_path_data = str(out_path / basename)
    logger.info("    getting excerpt to {}", out_path_data)
    excerpt_fun(url.url, str(out_path_data), url.excerpt_strategy.count)

    logger.info(" ... done downloading excerpt")


def copy_excerpt(url: UrlEntry, data_dir: str, output_document: str):
    """Copy downloaded excerpt to output file."""
    logger.info("    copying excerpt for {} to {}", url.url, output_document)
    in_path = pathlib.Path(data_dir) / url.hash
    parsed = urllib.parse.urlparse(url.url)
    basename = parsed.path.split("/")[-1]
    excerpt_path = in_path / basename
    shutil.copy(excerpt_path, output_document)
