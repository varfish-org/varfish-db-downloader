#!/usr/bin/env python
"""Helper script to retrieve ORDO gene-disease associations from OrphaData.

When run in CI mode environment variable "CI" is "true" then we will copy
the manually created file from excerpts/__orphadata__.  Note that this is
a deviation from the usual approach of putting the URL into `download_urls`.
"""

import json
import os
import sys

import httpx
import trio
from loguru import logger

#: URL for listing ORPHAcode.
URL_ORPHACODE_LIST = "https://api.orphadata.com/rd-cross-referencing/orphacodes"
#: URL template for getting information on one ORPHAcode.
URL_ORPHACODE_GET = "https://api.orphadata.com/rd-cross-referencing/orphacodes/{}?lang=en"
#: URL template for getting enes on one ORPHAcode.
URL_ORPHACODE_GET_GENE = "https://api.orphadata.com/rd-associated-genes/orphacodes/{}"

#: Whether we run in test mode.
RUNS_IN_CI = os.environ.get("CI", "false") == "true"


async def main():
    async with httpx.AsyncClient() as client:
        logger.info("Fetching ORPHAcode list...")
        lst = await client.get(URL_ORPHACODE_LIST)
        logger.info("...done")
        disease_ids = {disease["ORPHAcode"] for disease in lst.json()["data"]["results"]}

    async def work(no: int, orpha_id: int, limiter: trio.CapacityLimiter):
        async with limiter:
            async with httpx.AsyncClient() as client:
                try:
                    cross_references = (await client.get(URL_ORPHACODE_GET.format(orpha_id))).json()
                    disease_genes = (
                        await client.get(URL_ORPHACODE_GET_GENE.format(orpha_id), timeout=60)
                    ).json()
                except Exception as e:
                    logger.error(f"Error fetching {orpha_id}: {e}")
                    raise
                finally:
                    if no % 100 == 0:
                        logger.info(f"done fetching ORPHAcode details {no}/{len(disease_ids)}")
        json.dump(
            {
                "orpha_id": f"ORPHA:{orpha_id}",
                "cross_references": cross_references["data"]["results"],
                "disease_genes": disease_genes.get("data", {}).get("results"),
            },
            sys.stdout,
        )
        sys.stdout.write("\n")

    logger.info("Fetching ORPHAcode details...")
    limiter = trio.CapacityLimiter(10)
    async with trio.open_nursery() as nursery:
        for no, disease_id in enumerate(disease_ids):
            nursery.start_soon(work, no, disease_id, limiter)
    logger.info("...done")


if __name__ == "__main__":
    if RUNS_IN_CI:
        with open("excerpt-data/__orphadata__/orphadata.jsonl", "rt") as inputf:
            sys.stdout.write(inputf.read())
    else:
        trio.run(main)
