"""
Validate that each CITATIONS.md entry points to the intended paper.

Each citation entry in CITATIONS.md must include these required fields:
  - **Title:** <expected paper title>
  - **First Author:** <expected first-author family name>
  - **Year:** <expected publication year>
  - **Link:** <URL to the paper/landing page>

Checks performed:
1) Link resolution
   - FAIL on hard errors (404/410/5xx, DNS failures, timeouts).
   - WARN (but do not fail) on HTTP 403 since many publisher pages block automated requests.

2) Intended-paper validation (metadata match)
   - Title: fuzzy match (default threshold 0.90) against authoritative metadata.
   - Year: exact match.
   - First author family name: must match an author in authoritative metadata (case-insensitive).

Metadata sources:
  - DOI links (doi.org): Crossref is treated as authoritative for title/year/authors.
  - HTML pages (PMC/PubMed/Frontiers/Springer/etc): citation_* / dc.* / og:* meta tags.
  - PDF links: title/author/year heuristics from first-page text (pypdf).

Usage:
  python check_citations_links.py
  python check_citations_links.py --file CITATIONS.md --threshold 0.9
"""

from __future__ import annotations

import argparse
import difflib
import io
import json
import re
import sys
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from html.parser import HTMLParser
from typing import Sequence

from pypdf import PdfReader


REQUIRED_FIELDS = ("Title", "First Author", "Year", "Link")


@dataclass(frozen=True)
class CitationEntry:
    header: str
    fields: dict[str, str]


@dataclass(frozen=True)
class PaperMetadata:
    title: str | None
    year: int | None
    authors: list[str]
    doi: str | None = None
    source: str | None = None


class _MetaTagParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.meta: dict[str, list[str]] = {}

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "meta":
            return
        attr = {k.lower(): (v or "") for k, v in attrs}
        key = (attr.get("name") or attr.get("property") or "").strip().lower()
        val = (attr.get("content") or "").strip()
        if not key or not val:
            return
        self.meta.setdefault(key, []).append(val)


def _request(url: str, timeout_s: float) -> tuple[int, str]:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (citation-checker)"})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        return int(resp.status), str(resp.geturl())


def _fetch(url: str, timeout_s: float, max_bytes: int = 2_000_000) -> tuple[int, str, str, bytes]:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (citation-checker)"})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        status = int(resp.status)
        final_url = str(resp.geturl())
        ctype = (resp.headers.get("content-type") or "").lower()
        body = resp.read(max_bytes)
    return status, final_url, ctype, body


def _is_doi_url(url: str) -> bool:
    parsed = urllib.parse.urlparse(url)
    return parsed.scheme in {"http", "https"} and parsed.netloc.lower() in {"doi.org", "dx.doi.org"}


def _extract_doi(url: str) -> str | None:
    if not _is_doi_url(url):
        return None
    parsed = urllib.parse.urlparse(url)
    doi = parsed.path.lstrip("/")
    return urllib.parse.unquote(doi) if doi else None


def _extract_doi_from_url(url: str) -> str | None:
    """
    Best-effort DOI extraction from non-doi.org URLs.

    Examples:
      - https://www.frontiersin.org/articles/10.3389/fnhum.2018.00505/full -> 10.3389/fnhum.2018.00505
      - https://link.springer.com/book/10.1007/b98882 -> 10.1007/b98882
    """
    cleaned = urllib.parse.unquote(url)
    cleaned = cleaned.split("?", 1)[0].split("#", 1)[0]
    m = re.search(r"(10\.\d{4,9}/[-._;()/:A-Z0-9]+)", cleaned, flags=re.I)
    if not m:
        return None
    doi = m.group(1)
    # Heuristic trimming for publisher URL wrappers (e.g., Frontiers appends /full).
    for suffix in ("/full", "/abstract", "/pdf", "/epdf", "/html"):
        if suffix in doi.lower():
            doi = doi[: doi.lower().index(suffix)]
    return doi.rstrip("/")


def _normalize(s: str) -> str:
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    s = s.lower()
    s = re.sub(r"[^a-z0-9\s]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _similarity(a: str, b: str) -> float:
    return difflib.SequenceMatcher(None, _normalize(a), _normalize(b)).ratio()


def parse_citations_md(text: str) -> list[CitationEntry]:
    entries: list[CitationEntry] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r"^- \*\*(.+?)\*\*\s*$", lines[i])
        if not m:
            i += 1
            continue
        header = m.group(1)
        fields: dict[str, str] = {}
        i += 1
        while i < len(lines):
            if re.match(r"^- \*\*(.+?)\*\*\s*$", lines[i]) or re.match(r"^#\s", lines[i]):
                break
            fm = re.match(r"^\s{2}- \*\*(.+?):\*\*\s*(.*)\s*$", lines[i])
            if fm:
                key = fm.group(1).strip()
                val = fm.group(2).strip()
                if key and val:
                    fields[key] = val
            i += 1
        entries.append(CitationEntry(header=header, fields=fields))
    return entries


def _crossref_fetch(doi: str, timeout_s: float) -> PaperMetadata | None:
    api = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(api, headers={"User-Agent": "Mozilla/5.0 (citation-checker)"})
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            if int(resp.status) != 200:
                return None
            payload = json.loads(resp.read().decode("utf-8", errors="replace"))
    except Exception:
        return None

    msg = payload.get("message", {}) or {}
    title = (msg.get("title") or [None])[0]
    issued = (msg.get("issued") or {}).get("date-parts") or []
    year = int(issued[0][0]) if issued and issued[0] and issued[0][0] else None
    authors = msg.get("author") or []
    families = [a.get("family") for a in authors if a.get("family")]
    return PaperMetadata(
        title=str(title) if title else None,
        year=year,
        authors=[str(a) for a in families],
        doi=doi,
        source="crossref",
    )


def _html_meta_extract(body: bytes) -> dict[str, list[str]]:
    parser = _MetaTagParser()
    parser.feed(body.decode("utf-8", errors="ignore"))
    return parser.meta


def _paper_from_html_meta(url: str, timeout_s: float) -> PaperMetadata | None:
    try:
        _, _, _, body = _fetch(url, timeout_s=timeout_s, max_bytes=2_000_000)
    except Exception:
        return None

    meta = _html_meta_extract(body)

    def first(*keys: str) -> str | None:
        for k in keys:
            vals = meta.get(k.lower())
            if vals:
                return vals[0]
        return None

    title = first("citation_title", "dc.title", "og:title")
    doi = first("citation_doi", "dc.identifier")

    date = first("citation_date", "citation_publication_date", "citation_online_date", "dc.date")
    year = None
    if date:
        m = re.search(r"(19|20)\d{2}", date)
        if m:
            year = int(m.group(0))

    authors: list[str] = []
    if (vals := meta.get("citation_author")):
        authors = vals
    elif (vals := meta.get("citation_authors")):
        if len(vals) == 1 and ";" in vals[0]:
            authors = [a.strip() for a in vals[0].split(";") if a.strip()]
        else:
            authors = vals
    elif (vals := meta.get("dc.creator")):
        authors = vals

    return PaperMetadata(title=title, year=year, authors=authors, doi=doi, source="html-meta")


def _paper_from_jsonld(url: str, timeout_s: float) -> PaperMetadata | None:
    try:
        _, _, _, body = _fetch(url, timeout_s=timeout_s, max_bytes=5_000_000)
    except Exception:
        return None

    html = body.decode("utf-8", errors="ignore")
    scripts = re.findall(r'<script[^>]+type=["\\\']application/ld\+json["\\\'][^>]*>(.*?)</script>', html, re.I | re.S)
    for raw in scripts:
        raw = raw.strip()
        if not raw:
            continue
        try:
            data = json.loads(raw)
        except Exception:
            continue

        # JSON-LD can be a list or a dict.
        candidates = data if isinstance(data, list) else [data]
        for obj in candidates:
            if not isinstance(obj, dict):
                continue
            title = obj.get("name") or obj.get("headline")
            date = obj.get("datePublished") or obj.get("dateCreated") or obj.get("dateModified")
            year = None
            if isinstance(date, str):
                m = re.search(r"(19|20)\d{2}", date)
                if m:
                    year = int(m.group(0))

            authors: list[str] = []
            author_obj = obj.get("author")
            if isinstance(author_obj, list):
                for a in author_obj:
                    if isinstance(a, dict) and isinstance(a.get("name"), str):
                        authors.append(a["name"])
                    elif isinstance(a, str):
                        authors.append(a)
            elif isinstance(author_obj, dict) and isinstance(author_obj.get("name"), str):
                authors.append(author_obj["name"])
            elif isinstance(author_obj, str):
                authors.append(author_obj)

            if title or year or authors:
                return PaperMetadata(
                    title=str(title) if title else None,
                    year=year,
                    authors=authors,
                    doi=None,
                    source="json-ld",
                )
    return None


def _paper_from_pdf(url: str, timeout_s: float) -> PaperMetadata | None:
    try:
        _, _, _, body = _fetch(url, timeout_s=timeout_s, max_bytes=50_000_000)
    except Exception:
        return None
    try:
        reader = PdfReader(io.BytesIO(body))
        text = reader.pages[0].extract_text() or ""
    except Exception:
        return None

    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    title = lines[0] if lines else None
    authors_line = lines[1] if len(lines) > 1 else ""
    authors = [a.strip() for a in re.split(r",|\\sand\\s", authors_line) if a.strip()]
    year = None
    m = re.search(r"(19|20)\d{2}", text)
    if m:
        year = int(m.group(0))
    return PaperMetadata(title=title, year=year, authors=authors, doi=None, source="pdf-first-page")


def fetch_paper_metadata(url: str, timeout_s: float) -> PaperMetadata | None:
    doi = _extract_doi(url) or _extract_doi_from_url(url)
    crossref_meta: PaperMetadata | None = None
    if doi is not None:
        crossref_meta = _crossref_fetch(doi, timeout_s=timeout_s)
        # If Crossref has enough metadata, return early.
        if crossref_meta and crossref_meta.title and crossref_meta.year and crossref_meta.authors:
            return crossref_meta

    # Determine if PDF via content-type, url shape, or known sources.
    try:
        _, final_url, ctype, _ = _fetch(url, timeout_s=timeout_s, max_bytes=1024 * 64)
    except Exception:
        return None

    if (
        "application/pdf" in ctype
        or final_url.lower().endswith(".pdf")
        or url.startswith("https://par.nsf.gov/servlets/purl/")
    ):
        return _paper_from_pdf(url, timeout_s=timeout_s)

    # Try HTML meta and JSON-LD, and merge with Crossref if present.
    html_meta = _paper_from_html_meta(url, timeout_s=timeout_s)
    jsonld_meta = _paper_from_jsonld(url, timeout_s=timeout_s)

    # Merge priority: title/year prefer Crossref when available; authors prefer Crossref if present, else JSON-LD, else HTML meta.
    title = (crossref_meta.title if crossref_meta and crossref_meta.title else None) or (html_meta.title if html_meta else None) or (
        jsonld_meta.title if jsonld_meta else None
    )
    year = (crossref_meta.year if crossref_meta and crossref_meta.year else None) or (html_meta.year if html_meta else None) or (
        jsonld_meta.year if jsonld_meta else None
    )
    authors = []
    if crossref_meta and crossref_meta.authors:
        authors = crossref_meta.authors
    elif jsonld_meta and jsonld_meta.authors:
        authors = jsonld_meta.authors
    elif html_meta and html_meta.authors:
        authors = html_meta.authors

    doi_final = crossref_meta.doi if crossref_meta and crossref_meta.doi else (html_meta.doi if html_meta else None)
    source = "merged"
    return PaperMetadata(title=title, year=year, authors=authors, doi=doi_final, source=source)


def validate_citations(entries: Sequence[CitationEntry], timeout_s: float, threshold: float) -> tuple[int, int, int]:
    ok = 0
    warn = 0
    fail = 0

    for entry in entries:
        missing = [f for f in REQUIRED_FIELDS if f not in entry.fields]
        if missing:
            fail += 1
            print(f"FAIL (missing fields) {entry.header}: {', '.join(missing)}")
            continue

        expected_title = entry.fields["Title"]
        expected_author = entry.fields["First Author"]
        try:
            expected_year = int(entry.fields["Year"])
        except ValueError:
            fail += 1
            print(f"FAIL (bad year) {entry.header}: Year='{entry.fields['Year']}'")
            continue

        url = entry.fields["Link"]

        # Link resolution (warn on 403)
        try:
            status, final_url = _request(url, timeout_s=timeout_s)
        except urllib.error.HTTPError as e:
            status = int(e.code)
            final_url = getattr(e, "url", None) or url
        except Exception as e:
            fail += 1
            print(f"FAIL (link) {entry.header}: {type(e).__name__}: {e}")
            continue

        if status == 403:
            warn += 1
            print(f"WARN (403) {entry.header}: {url} -> {final_url}")
        elif not (200 <= status < 400):
            fail += 1
            print(f"FAIL (link) {entry.header}: HTTP {status} {url} -> {final_url}")
            continue

        meta = fetch_paper_metadata(url, timeout_s=timeout_s)
        if meta is None or not meta.title or meta.year is None:
            fail += 1
            print(f"FAIL (metadata) {entry.header}: unable to extract title/year from link")
            continue

        title_score = _similarity(expected_title, meta.title)
        if title_score < threshold:
            fail += 1
            print(
                f"FAIL (title) {entry.header}: score={title_score:.3f} < {threshold:.2f}\n"
                f"  expected: {expected_title}\n"
                f"  actual:   {meta.title} ({meta.source})"
            )
            continue

        if meta.year != expected_year:
            fail += 1
            print(f"FAIL (year) {entry.header}: expected {expected_year}, got {meta.year} ({meta.source})")
            continue

        exp = _normalize(expected_author)
        author_norms = [_normalize(a) for a in meta.authors]
        author_ok = any(exp == a for a in author_norms) or any(exp in a for a in author_norms)
        if not author_ok:
            fail += 1
            preview = ", ".join(meta.authors[:8]) if meta.authors else "(no authors parsed)"
            print(f"FAIL (author) {entry.header}: expected '{expected_author}', got [{preview}] ({meta.source})")
            continue

        ok += 1
        print(f"OK   {entry.header}: title={title_score:.3f}, year={meta.year}, author={expected_author}")

    return ok, warn, fail


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default="CITATIONS.md", help="Path to CITATIONS markdown file.")
    parser.add_argument("--timeout", type=float, default=30.0, help="Per-link timeout in seconds.")
    parser.add_argument("--threshold", type=float, default=0.90, help="Title fuzzy-match threshold (0–1).")
    args = parser.parse_args(argv)

    try:
        text = open(args.file, "r", encoding="utf-8").read()
    except FileNotFoundError:
        print(f"FAIL ERR {args.file} (file not found)")
        return 2

    entries = parse_citations_md(text)
    if not entries:
        print("FAIL ERR (no citation entries found)")
        return 2

    ok, warn, fail = validate_citations(entries, timeout_s=float(args.timeout), threshold=float(args.threshold))
    print(f"\nSummary: {ok} ok, {warn} warn (403), {fail} fail")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
