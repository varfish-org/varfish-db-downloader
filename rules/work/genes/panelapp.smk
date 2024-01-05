## Rules related to AlphaMissense per-gene scores


import os
import subprocess
import sys
import tempfile


rule genes_panelapp_download:  # -- download AlphaMissense per-gene scores
    output:
        jsonl="work/download/genes/panelapp/{date}/panelapp.jsonl",
    run:
        base_url = "https://panelapp.genomicsengland.co.uk/api/v1"

        pages = []
        page_no = 1
        page_count = None
        url = f"{base_url}/entities/"
        with tempfile.TemporaryDirectory() as tmpdir:
            while url:
                print(
                    f"downloading page {page_no}/{page_count if page_count else '?'}...",
                    file=sys.stderr,
                )
                subprocess.check_call(["wget", "-O", f"{tmpdir}/page.json", url])
                with open(f"{tmpdir}/page.json", "rt") as f:
                    page = json.load(f)
                pages.append(page)
                url = page.get("next")
                page_no += 1
                if not page_count:
                    per_page = len(page.get("results", [None]))
                    page_count = (page.get("count") + per_page - 1) // per_page
                if os.environ.get("CI", None) == "true" and page_no > 2:
                    print("CI mode: only downloading first 2 pages", file=sys.stderr)
                    break

        os.makedirs(f"work/download/genes/panelapp/{wildcards.date}", exist_ok=True)
        with open(output.jsonl, "wt") as f:
            for page in pages:
                for result in page.get("results", []):
                    print(json.dumps(result), file=f)
