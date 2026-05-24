#!/usr/bin/env python3
"""Export numbered references from manuscript.md to RIS format.

The parser is intentionally conservative and is tailored to the manuscript's
Vancouver-style references:

    1. Authors. Title. Journal. Year;volume(issue):pages. doi:...

Conference references with "In:" are exported as CPAPER entries.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
INPUT = BASE / "manuscript.md"
OUTPUT = BASE / "references.ris"

content = INPUT.read_text(encoding="utf-8")
if "## References" not in content:
    raise SystemExit("No References section found")

ref_section = content.split("## References", 1)[1]
ref_lines = [
    line.strip()
    for line in ref_section.splitlines()
    if re.match(r"^\d+\.\s+", line.strip())
]


def split_authors(authors_str: str) -> list[str]:
    authors_str = authors_str.replace(" Jr,", " Jr,")
    authors_str = re.sub(r",?\s+et al\.?", "", authors_str).strip()
    raw_authors = [a.strip() for a in authors_str.split(",") if a.strip()]
    authors: list[str] = []
    i = 0
    while i < len(raw_authors):
        raw = raw_authors[i]
        # Already appears as "Family Initials" in Vancouver style.
        parts = raw.rsplit(" ", 1)
        if len(parts) == 2:
            family, initials = parts
            authors.append(f"{family}, {initials}")
        else:
            authors.append(raw)
        i += 1
    return authors


def parse_ref(line: str) -> dict[str, object]:
    line = re.sub(r"^\d+\.\s*", "", line).strip()
    doi = ""
    doi_match = re.search(r"\sdoi:([^\s.]+(?:\.[^\s.]+)*)\.?$", line)
    if doi_match:
        doi = doi_match.group(1).rstrip(".")
        line_wo_doi = line[: doi_match.start()].rstrip(". ")
    else:
        line_wo_doi = line.rstrip(".")

    # Conference/proceedings reference.
    if " In: " in line_wo_doi:
        authors_part, rest = line_wo_doi.split(". ", 1)
        title, conf_rest = rest.split(" In: ", 1)
        year_match = re.search(r";\s*(\d{4})\.\s*p\.\s*([0-9A-Za-z\-–]+)", conf_rest)
        year = year_match.group(1) if year_match else ""
        pages = year_match.group(2).replace("–", "-") if year_match else ""
        booktitle = conf_rest[: year_match.start()].rstrip(". ") if year_match else conf_rest
        return {
            "type": "CPAPER",
            "authors": split_authors(authors_part),
            "title": title,
            "journal": "",
            "booktitle": booktitle,
            "year": year,
            "volume": "",
            "issue": "",
            "pages": pages,
            "doi": doi,
        }

    # Split at final bibliographic year pattern.
    year_match = re.search(r"(\d{4});([^.]*)$", line_wo_doi)
    if not year_match:
        # Handle preprint style: bioRxiv. 2025:identifier
        preprint_match = re.search(r"(\d{4}):(.+)$", line_wo_doi)
        if preprint_match:
            before = line_wo_doi[: preprint_match.start()].rstrip(". ")
            year = preprint_match.group(1)
            pages = preprint_match.group(2)
            parts = before.split(". ")
            authors_part = parts[0]
            journal = parts[-1]
            title = ". ".join(parts[1:-1])
            return {
                "type": "JOUR",
                "authors": split_authors(authors_part),
                "title": title,
                "journal": journal,
                "booktitle": "",
                "year": year,
                "volume": "",
                "issue": "",
                "pages": pages,
                "doi": doi,
            }
        raise ValueError(f"Could not parse reference: {line}")

    before_year = line_wo_doi[: year_match.start()].rstrip(". ")
    year = year_match.group(1)
    after_year = year_match.group(2).strip()
    parts = before_year.split(". ")
    authors_part = parts[0]
    journal = parts[-1]
    title = ". ".join(parts[1:-1])

    volume = issue = pages = ""
    vol_match = re.match(r"([^:(]+)(?:\(([^)]*)\))?(?::(.+))?", after_year)
    if vol_match:
        volume = vol_match.group(1).strip()
        issue = (vol_match.group(2) or "").strip()
        pages = (vol_match.group(3) or "").strip()

    return {
        "type": "JOUR",
        "authors": split_authors(authors_part),
        "title": title,
        "journal": journal,
        "booktitle": "",
        "year": year,
        "volume": volume,
        "issue": issue,
        "pages": pages,
        "doi": doi,
    }


def ris_entry(ref: dict[str, object]) -> str:
    lines = [f"TY  - {ref['type']}"]
    if ref.get("title"):
        lines.append(f"TI  - {ref['title']}")
    for author in ref.get("authors", []):
        lines.append(f"AU  - {author}")
    if ref.get("journal"):
        lines.append(f"JO  - {ref['journal']}")
        lines.append(f"JF  - {ref['journal']}")
    if ref.get("booktitle"):
        lines.append(f"T2  - {ref['booktitle']}")
    if ref.get("year"):
        lines.append(f"PY  - {ref['year']}")
        lines.append(f"Y1  - {ref['year']}")
    if ref.get("volume"):
        lines.append(f"VL  - {ref['volume']}")
    if ref.get("issue"):
        lines.append(f"IS  - {ref['issue']}")
    pages = str(ref.get("pages") or "")
    if pages:
        if "-" in pages:
            sp, ep = pages.split("-", 1)
            lines.append(f"SP  - {sp}")
            lines.append(f"EP  - {ep}")
        else:
            lines.append(f"SP  - {pages}")
    if ref.get("doi"):
        lines.append(f"DO  - {ref['doi']}")
    lines.append("ER  -")
    return "\n".join(lines)

entries = []
for ref_line in ref_lines:
    entries.append(ris_entry(parse_ref(ref_line)))

OUTPUT.write_text("\n\n".join(entries) + "\n", encoding="utf-8")
print(f"Exported {len(entries)} references to {OUTPUT}")
print(f"File size: {os.path.getsize(OUTPUT)} bytes")
