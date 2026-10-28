import importlib.util
import io
import json
import unittest

import httpx
import trio

spec = importlib.util.spec_from_file_location(
    "genes_orpha_diseases", "scripts/genes-orpha-diseases.py"
)
genes_orpha_diseases = importlib.util.module_from_spec(spec)
spec.loader.exec_module(genes_orpha_diseases)


def orphadata_api(failing_urls: set[str]) -> httpx.MockTransport:
    """Mock api.orphadata.com; the first request to each URL in ``failing_urls`` is dropped."""

    def handler(request: httpx.Request) -> httpx.Response:
        url = str(request.url)
        if url in failing_urls:
            failing_urls.remove(url)
            raise httpx.ConnectError("EOF occurred in violation of protocol", request=request)
        if url == genes_orpha_diseases.URL_ORPHACODE_LIST:
            return httpx.Response(200, json={"data": {"results": [{"ORPHAcode": 5}]}})
        if url == genes_orpha_diseases.URL_ORPHACODE_GET.format(5):
            return httpx.Response(200, json={"data": {"results": {"Preferred term": "Disease 5"}}})
        if url == genes_orpha_diseases.URL_ORPHACODE_GET_GENE.format(5):
            return httpx.Response(200, json={"data": {"results": {"Gene": "GENE5"}}})
        return httpx.Response(404)

    return httpx.MockTransport(handler)


class DownloadTest(unittest.TestCase):
    def setUp(self):
        genes_orpha_diseases.RETRY_DELAY = 0.0

    def test_dropped_connections_are_retried(self):
        failing_urls = {
            genes_orpha_diseases.URL_ORPHACODE_GET.format(5),
            genes_orpha_diseases.URL_ORPHACODE_GET_GENE.format(5),
        }
        out = io.StringIO()

        async def run():
            async with httpx.AsyncClient(transport=orphadata_api(failing_urls)) as client:
                await genes_orpha_diseases.download(client, out)

        trio.run(run)

        self.assertEqual(
            [json.loads(line) for line in out.getvalue().splitlines()],
            [
                {
                    "orpha_id": "ORPHA:5",
                    "cross_references": {"Preferred term": "Disease 5"},
                    "disease_genes": {"Gene": "GENE5"},
                }
            ],
        )


if __name__ == "__main__":
    unittest.main()
