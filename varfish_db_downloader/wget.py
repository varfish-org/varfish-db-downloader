"""Implementation of the stub ``wget`` and supporting tool."""

import gzip
import hashlib
import itertools
import os
import pathlib
import shutil
import subprocess
import tempfile
import typing
import urllib.parse
import zlib

import attrs
import cattrs
import click
import requests
import requests_ftp
import yaml
from loguru import logger


def excerpt_manual(url: str, path_out: str, count: int):
    """Do not download, assume output is already there."""
    _ = url
    _ = count
    logger.info("    (strategy MANUAL from {} to {})", url, path_out)
    if not os.path.exists(path_out):
        raise RuntimeError(f"File {path_out} does not exist (strategy MANUAL)")


def no_excerpt(url: str, path_out: str, count: int):
    """Do not excerpt, use all."""
    logger.info("    (strategy no-excerpt from {} to {})", url, path_out)
    _ = count
    with open(path_out, "wb") as outputf:
        requests_ftp.monkeypatch_session()
        s = requests.Session()
        r = s.get(url, allow_redirects=True)
        outputf.write(r.content)


def decompress_stream(stream):
    """Helper for decompressing data from a stream."""
    o = zlib.decompressobj(16 + zlib.MAX_WBITS)
    for chunk in stream:
        yield o.decompress(chunk)
    yield o.flush()


def excerpt_head(url: str, path_out: str, count: int):
    """Excerpt a plaint-text file by copying lines."""
    logger.info("    (strategy head from {} to {})", url, path_out)
    try_gzip = url.endswith(".gz") or url.endswith(".bgz")
    if try_gzip:
        opener = gzip.open
    else:
        opener = open
    with opener(path_out, "wb") as f_out:
        requests_ftp.monkeypatch_session()
        s = requests.Session()
        r = s.get(url, stream=True)
        # First, attempt to read line by line which should work if requests is
        # correctly identifying gzip compression.
        is_raw_gzip = False
        for i, line in enumerate(r.iter_lines()):
            if i == 0 and line.startswith(b"\x1f\x8b"):
                is_raw_gzip = True
                break
            if i >= count:
                break
            f_out.write(line)
            f_out.write(b"\n")
        if is_raw_gzip:
            logger.info("    falling back to raw gzip decompression")
            # If requests failed to detect gzip compression, we need to decompress
            # manually.
            r = s.get(url, stream=True)
            parseable_data = decompress_stream(r.iter_content(1024))
            hacky_n_chunks = 100  # we will read only 10 chunks for now
            chunks_byte = list(itertools.islice(parseable_data, hacky_n_chunks))
            collapsed_byte = b"".join(chunks_byte)
            collapsed = collapsed_byte.decode("utf-8")
            for i, line in enumerate(collapsed.split("\n")):
                if i >= count:
                    break
                f_out.write(line.encode("utf-8"))
                f_out.write(b"\n")
    if try_gzip:
        # Fixup resulting compressed files so they are proper bgzip.
        subprocess.check_call(["gzip", "-d", path_out])
        subprocess.check_call(["bgzip", path_out.replace(".gz", "")])


def excerpt_copy_tbi(url: str, path_out: str, count: int):
    """Copy ``.tbi`` file from the (previously downloaded) ``.vcf.gz`` file.)"""
    _ = count
    vcf_url = url[:-4]
    vcf_file = vcf_url.split("/")[-1]
    vcf_url_hash = hash_url(vcf_url)
    base_path = os.path.dirname(os.path.dirname(path_out))
    vcf_tbi_path = f"{base_path}/{vcf_url_hash}/{vcf_file}.tbi"
    logger.info("    (strategy COPY_TBI from {} to {})", vcf_tbi_path, path_out)
    if not os.path.exists(vcf_tbi_path):
        raise RuntimeError(f"File {vcf_tbi_path} does not exist (strategy COPY-TBI)")
    shutil.copy(vcf_tbi_path, path_out)


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
        logger.info("    + {}", lst_cmd)
        chroms_out = subprocess.check_output(lst_cmd, shell=True, cwd=tmpdir)
        chroms_lst = [
            line.strip() for line in chroms_out.decode("utf-8").split("\n") if line.strip()
        ]
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
            logger.info("    + {}", cmd)
            subprocess.check_call(cmd, shell=True, cwd=tmpdir)


#: Length of the hash to use.  The hash is created from the underlying URL.
HASH_LENGTH = 16

#: Mapping from strategy name to implementation.
STRATEGY_MAP = {
    "head": excerpt_head,
    "gz-head": excerpt_head,
    "vcf-head": excerpt_vcf_head,
    "copy-tbi": excerpt_copy_tbi,
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
        elif url.endswith(".tbi"):
            return ExcerptStrategy(strategy="copy-tbi", count=None)
        elif "." not in url.split("/")[-1]:
            return ExcerptStrategy(strategy="no-excerpt", count=None)
        else:
            return ExcerptStrategy(strategy="head", count=100)


def hash_url(url: str) -> str:
    """Return hashed URL string."""
    return hashlib.sha256(url.encode("utf-8")).hexdigest()[:HASH_LENGTH]


@attrs.frozen()
class UrlEntry:
    """Entry in the ``downloads_urls.yml`` file."""

    #: The URL of the original file download.
    url: str
    #: The hash to use (computed from ``url``).
    hash: typing.Optional[str] = None
    #: The excerpt strategy to use.
    excerpt_strategy: typing.Optional[ExcerptStrategy] = None
    #: Whether or not to skip the upstream check.
    skip_upstream_check: bool = False

    def __attrs_post_init__(self):
        object.__setattr__(self, "hash", hash_url(self.url))


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
    logger.info("    + mkdir -p {}", out_path)
    out_path_url = out_path / "url.txt"
    logger.info("    writing URL to {}", out_path_url)
    with out_path_url.open("wt") as f:
        print(url.url, file=f)

    excerpt_fun = STRATEGY_MAP[url.excerpt_strategy.strategy]
    parsed = urllib.parse.urlparse(url.url)
    basename = parsed.path.split("/")[-1] or "__index__"
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
    click.echo(err=True, message="copying {} => {}".format(excerpt_path, output_document))
    if os.path.isdir(excerpt_path):
        shutil.copy(f"{excerpt_path}/__index__", output_document)
    else:
        shutil.copy(excerpt_path, output_document)
