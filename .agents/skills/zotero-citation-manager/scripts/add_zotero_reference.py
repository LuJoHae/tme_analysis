#!/usr/bin/env python3
"""
Pure functional utility to convert Zotero item metadata into clean BibTeX
and append it to the article bibliography without duplicates.

Adheres strictly to the project's functional programming principles:
- Immutable data structures (frozen models)
- Railway-oriented error handling via Result monad (bind, map)
- Pure transformations without hidden side effects
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Generic, Optional, TypeVar, Union

T = TypeVar("T")
E = TypeVar("E")
U = TypeVar("U")


# 1. Functional Result Monad
@dataclass(frozen=True)
class Success(Generic[T]):
    value: T

    def is_success(self) -> bool:
        return True

    def map(self, fn: Callable[[T], U]) -> Success[U]:
        return Success(fn(self.value))

    def bind(self, fn: Callable[[T], Result[U, E]]) -> Result[U, E]:
        return fn(self.value)


@dataclass(frozen=True)
class Failure(Generic[E]):
    error: E

    def is_success(self) -> bool:
        return False

    def map(self, fn: Callable[..., Any]) -> Failure[E]:
        return self

    def bind(self, fn: Callable[..., Any]) -> Failure[E]:
        return self


Result = Union[Success[T], Failure[E]]


# 2. Immutable Data Models
@dataclass(frozen=True)
class Creator:
    first_name: str
    last_name: str
    creator_type: str


@dataclass(frozen=True)
class ZoteroItem:
    key: str
    item_type: str
    title: str
    creators: tuple[Creator, ...]
    publication_title: str
    publisher: str
    date: str
    volume: str
    issue: str
    pages: str
    doi: str
    url: str
    citation_key: str


@dataclass(frozen=True)
class BibTeXEntry:
    citekey: str
    entry_text: str
    doi: str


# 3. Pure Parsing & Transformation Functions
def parse_creator(raw: dict[str, Any]) -> Creator:
    creator_type = str(raw.get("creatorType", "author"))
    if "lastName" in raw and "firstName" in raw:
        return Creator(
            first_name=str(raw["firstName"]),
            last_name=str(raw["lastName"]),
            creator_type=creator_type
        )
    if "name" in raw:
        parts = str(raw["name"]).split(" ", 1)
        return Creator(
            first_name=parts[0] if len(parts) > 1 else "",
            last_name=parts[1] if len(parts) > 1 else parts[0],
            creator_type=creator_type
        )
    return Creator(first_name="", last_name="", creator_type=creator_type)


def parse_zotero_json(raw: dict[str, Any]) -> Result[ZoteroItem, str]:
    if not isinstance(raw, dict):
        return Failure("Provided payload is not a valid dictionary.")

    key = str(raw.get("key", ""))
    if not key:
        return Failure("Missing 'key' in Zotero payload.")

    creators_list = raw.get("creators", [])
    parsed_creators: tuple[Creator, ...] = tuple(
        parse_creator(c) for c in creators_list if isinstance(c, dict)
    )

    return Success(
        ZoteroItem(
            key=key,
            item_type=str(raw.get("itemType", "journalArticle")),
            title=str(raw.get("title", "")),
            creators=parsed_creators,
            publication_title=str(raw.get("publicationTitle", "")),
            publisher=str(raw.get("publisher", "")),
            date=str(raw.get("date", "")),
            volume=str(raw.get("volume", "")),
            issue=str(raw.get("issue", raw.get("number", ""))),
            pages=str(raw.get("pages", "")),
            doi=str(raw.get("DOI", raw.get("doi", ""))),
            url=str(raw.get("url", "")),
            citation_key=str(raw.get("citationKey", ""))
        )
    )


def extract_year(date_str: str) -> str:
    m = re.search(r"\b(19\d\d|20\d\d)\b", date_str)
    return m.group(1) if m else ""


def derive_citekey(item: ZoteroItem) -> str:
    if item.citation_key and re.match(r"^[a-zA-Z0-9_\-]+$", item.citation_key):
        return item.citation_key

    first_author_last = item.creators[0].last_name if item.creators else "Unknown"
    clean_author = re.sub(r"[^a-zA-Z]", "", first_author_last).lower()
    year = extract_year(item.date) or "nd"

    # Extract first significant word from title
    title_words = [w for w in re.sub(r"[^a-zA-Z0-9\s]", "", item.title).split() if len(w) > 3]
    keyword = title_words[0].lower() if title_words else "ref"

    return f"{clean_author}{year}{keyword}"


def format_authors_bibtex(creators: tuple[Creator, ...]) -> str:
    formatted: list[str] = []
    for c in creators:
        if c.last_name and c.first_name:
            formatted.append(f"{c.last_name}, {c.first_name}")
        elif c.last_name:
            formatted.append(c.last_name)
    return " and ".join(formatted) if formatted else "Unknown"


def item_type_to_bibtex_type(item_type: str) -> str:
    type_mapping = {
        "journalArticle": "article",
        "book": "book",
        "bookSection": "incollection",
        "conferencePaper": "inproceedings",
        "thesis": "phdthesis",
    }
    return type_mapping.get(item_type, "misc")


def build_bibtex_entry(item: ZoteroItem) -> Result[BibTeXEntry, str]:
    citekey = derive_citekey(item)
    bib_type = item_type_to_bibtex_type(item.item_type)
    authors_str = format_authors_bibtex(item.creators)
    year = extract_year(item.date)

    fields: list[str] = [
        f"  title = {{{item.title}}}",
        f"  author = {{{authors_str}}}"
    ]

    if item.publication_title:
        fields.append(
            f"  journal = {{{item.publication_title}}}"
            if bib_type == "article"
            else f"  booktitle = {{{item.publication_title}}}"
        )
    if item.publisher:
        fields.append(f"  publisher = {{{item.publisher}}}")
    if year:
        fields.append(f"  year = {{{year}}}")
    if item.volume:
        fields.append(f"  volume = {{{item.volume}}}")
    if item.issue:
        fields.append(f"  number = {{{item.issue}}}")
    if item.pages:
        fields.append(f"  pages = {{{item.pages}}}")
    if item.doi:
        fields.append(f"  doi = {{{item.doi}}}")
    if item.url:
        fields.append(f"  url = {{{item.url}}}")

    entry_lines = [f"@{bib_type}{{{citekey},"] + [",\n".join(fields)] + ["}\n"]
    entry_text = "\n".join(entry_lines)

    return Success(BibTeXEntry(citekey=citekey, entry_text=entry_text, doi=item.doi))


def check_existing_citekey(bib_content: str, citekey: str, doi: str) -> Optional[str]:
    # Check citekey exact match
    key_pattern = rf"@\w+\s*\{{\s*{re.escape(citekey)}\s*,"
    if re.search(key_pattern, bib_content, re.IGNORECASE):
        return f"Citation key '{citekey}' already exists in bibliography."

    # Check DOI match if present
    if doi:
        clean_doi = re.escape(doi.strip())
        doi_pattern = rf"doi\s*=\s*\{{\s*{clean_doi}\s*\}}"
        if re.search(doi_pattern, bib_content, re.IGNORECASE):
            return f"Entry with DOI '{doi}' already exists in bibliography."

    return None


def append_entry_to_bib_file(bib_path: Path, entry: BibTeXEntry) -> Result[str, str]:
    if not bib_path.exists():
        return Failure(f"Bibliography file does not exist: {bib_path}")

    try:
        content = bib_path.read_text(encoding="utf-8")
    except Exception as exc:
        return Failure(f"Failed to read {bib_path}: {exc}")

    existing_warning = check_existing_citekey(content, entry.citekey, entry.doi)
    if existing_warning:
        return Success(f"EXISTS: {existing_warning} Use citekey: @{entry.citekey}")

    # Pure append
    new_content = content.rstrip() + "\n\n" + entry.entry_text.strip() + "\n"
    try:
        bib_path.write_text(new_content, encoding="utf-8")
    except Exception as exc:
        return Failure(f"Failed to write to {bib_path}: {exc}")

    return Success(f"ADDED: Successfully added @{entry.citekey} to {bib_path}")


# 4. Pipeline Execution via Monadic Bind
def process_zotero_payload(raw_json: dict[str, Any], bib_path: Path) -> Result[dict[str, str], str]:
    return (
        parse_zotero_json(raw_json)
        .bind(lambda item: build_bibtex_entry(item))
        .bind(
            lambda entry: append_entry_to_bib_file(bib_path, entry).map(
                lambda msg: {
                    "status": msg,
                    "citekey": entry.citekey,
                    "typst_citation": f"@{entry.citekey}",
                    "bibtex": entry.entry_text
                }
            )
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Add Zotero item to article references.bib")
    parser.add_argument("--json", type=str, help="Raw JSON string from Zotero get_item_details")
    parser.add_argument("--file", type=Path, help="JSON file containing Zotero item details")
    parser.add_argument(
        "--bib-path",
        type=Path,
        default=Path("article/references.bib"),
        help="Path to target references.bib"
    )

    args = parser.parse_args()

    # Load input JSON
    raw_str = ""
    if args.json:
        raw_str = args.json
    elif args.file:
        raw_str = args.file.read_text(encoding="utf-8")
    elif not sys.stdin.isatty():
        raw_str = sys.stdin.read()
    else:
        print("Error: No input JSON provided via --json, --file, or stdin.", file=sys.stderr)
        sys.exit(1)

    try:
        raw_dict = json.loads(raw_str)
    except Exception as exc:
        print(f"Error: Invalid JSON payload: {exc}", file=sys.stderr)
        sys.exit(1)

    result = process_zotero_payload(raw_dict, args.bib_path)
    if isinstance(result, Failure):
        print(f"FAILED: {result.error}", file=sys.stderr)
        sys.exit(1)
    else:
        print(json.dumps(result.value, indent=2))
        sys.exit(0)


if __name__ == "__main__":
    main()
