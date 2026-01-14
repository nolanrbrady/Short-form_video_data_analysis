"""
Check that links in CITATIONS.md resolve to a real resource.

This script is intentionally conservative for scientific integrity:
- It *fails* on hard errors (404/410/5xx, DNS failures, timeouts).
- It *warns* (but does not fail) on 403, since many publisher pages block
  automated requests while still being valid for human readers.
- For DOI links (https://doi.org/<doi>), it additionally validates the DOI exists
  via Crossref, which is less susceptible to publisher bot-blocking.

Usage:
  python check_citations_links.py
  python check_citations_links.py --file CITATIONS.md
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.parse
import urllib.error
import urllib.request
from dataclasses import dataclass
from typing import Iterable, Sequence


URL_RE = re.compile(r"https?://[^\s)]+")


@dataclass(frozen=True)
class LinkResult:
    url: str
    ok: bool
    status: int | None
    final_url: str | None
    note: str | None = None


def extract_urls(text: str) -> list[str]:
    return URL_RE.findall(text)


def _request(url: str, timeout_s: float) -> tuple[int, str]:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (citation-link-checker)"})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        return int(resp.status), str(resp.geturl())


def _is_doi_url(url: str) -> bool:
    parsed = urllib.parse.urlparse(url)
    return parsed.scheme in {"http", "https"} and parsed.netloc.lower() in {"doi.org", "dx.doi.org"}


def _extract_doi(url: str) -> str | None:
    if not _is_doi_url(url):
        return None
    parsed = urllib.parse.urlparse(url)
    doi = parsed.path.lstrip("/")
    if not doi:
        return None
    return urllib.parse.unquote(doi)


def _crossref_check(doi: str, timeout_s: float) -> tuple[bool, str | None]:
    api = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(api, headers={"User-Agent": "Mozilla/5.0 (citation-link-checker)"})
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            if int(resp.status) != 200:
                return False, f"Crossref returned HTTP {resp.status}"
            payload = json.loads(resp.read().decode("utf-8", errors="replace"))
            title = (payload.get("message", {}) or {}).get("title", [])
            if title:
                return True, str(title[0])
            return True, None
    except Exception as e:  # noqa: BLE001 - surfaced to caller for hard-fail.
        return False, f"Crossref request failed: {type(e).__name__}: {e}"


def check_urls(urls: Sequence[str], timeout_s: float) -> list[LinkResult]:
    results: list[LinkResult] = []

    for url in urls:
        note_parts: list[str] = []

        doi = _extract_doi(url)
        if doi is not None:
            cr_ok, cr_note = _crossref_check(doi, timeout_s=timeout_s)
            if not cr_ok:
                results.append(LinkResult(url=url, ok=False, status=None, final_url=None, note=cr_note))
                continue
            if cr_note:
                note_parts.append(f"Crossref: {cr_note}")

        try:
            status, final_url = _request(url, timeout_s=timeout_s)
        except urllib.error.HTTPError as e:
            status = int(e.code)
            final_url = getattr(e, "url", None) or url
        except Exception as e:  # noqa: BLE001 - surfaced to caller for hard-fail.
            results.append(
                LinkResult(url=url, ok=False, status=None, final_url=None, note=f"{type(e).__name__}: {e}")
            )
            continue

        if status == 403:
            results.append(
                LinkResult(
                    url=url,
                    ok=True,
                    status=status,
                    final_url=final_url,
                    note="HTTP 403 (likely publisher bot-blocking); treat as warning",
                )
            )
            continue

        if 200 <= status < 400:
            note = "; ".join(note_parts) if note_parts else None
            results.append(LinkResult(url=url, ok=True, status=status, final_url=final_url, note=note))
            continue

        results.append(LinkResult(url=url, ok=False, status=status, final_url=final_url, note="HTTP error"))

    return results


def _print_results(results: Iterable[LinkResult]) -> tuple[int, int, int]:
    ok = 0
    warn = 0
    fail = 0

    for r in results:
        if r.ok:
            ok += 1
            if r.status == 403:
                warn += 1
        else:
            fail += 1

    for r in results:
        if r.ok and r.status == 403:
            print(f"WARN {r.status} {r.url} -> {r.final_url}")
            continue
        if r.ok:
            if r.note:
                print(f"OK   {r.status} {r.url} ({r.note})")
            else:
                print(f"OK   {r.status} {r.url}")
        else:
            detail = f"{r.status}" if r.status is not None else "ERR"
            print(f"FAIL {detail} {r.url} ({r.note})")

    return ok, warn, fail


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default="CITATIONS.md", help="Path to CITATIONS markdown file.")
    parser.add_argument("--timeout", type=float, default=30.0, help="Per-link timeout in seconds.")
    args = parser.parse_args(argv)

    try:
        text = open(args.file, "r", encoding="utf-8").read()
    except FileNotFoundError:
        print(f"FAIL ERR {args.file} (file not found)")
        return 2

    urls = extract_urls(text)
    if not urls:
        print("FAIL ERR (no URLs found)")
        return 2

    results = check_urls(urls, timeout_s=float(args.timeout))
    ok, warn, fail = _print_results(results)
    print(f"\nSummary: {ok} ok, {warn} warn (403), {fail} fail")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
