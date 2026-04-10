#!/usr/bin/env python3
"""PubMed daily paper fetcher.
Searches for recent papers by keyword via E-utilities,
outputs YAML-frontmatter Markdown files into papers/.
"""

import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError
from urllib.parse import quote

REPO_ROOT = Path(__file__).resolve().parent.parent
PAPERS_DIR = REPO_ROOT / "papers"
STATE_FILE = REPO_ROOT / ".fetch-state.json"

KEYWORDS = [
    ("cancer immunotherapy", "cancer"),
    ("gene therapy", "gene-therapy"),
    ("clinical trial artificial intelligence", "clinical-ai"),
]
DATE_RANGE = "2025:2026[dp]"
PER_KEYWORD = 10
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PMC_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Optional: set NCBI_API_KEY env var for higher rate limits
API_KEY = os.environ.get("NCBI_API_KEY", "")


def api_get(url, retries=3, as_xml=False):
    """GET with retry."""
    for attempt in range(retries):
        try:
            req = Request(url)
            with urlopen(req, timeout=30) as resp:
                data = resp.read()
                return data if as_xml else json.loads(data)
        except HTTPError as e:
            if e.code == 429:
                wait = 3 * (attempt + 1)
                print(f"  Rate limited, waiting {wait}s...")
                time.sleep(wait)
            else:
                raise
    return None


def slug(title):
    """Title -> filename-safe slug."""
    s = re.sub(r'[^\w\s-]', '', title.lower().strip())
    s = re.sub(r'[\s_]+', '-', s)
    return s[:80]


def load_state():
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    return {"fetched_pmids": [], "last_run": None}


def save_state(state):
    state["last_run"] = datetime.now().isoformat()
    STATE_FILE.write_text(json.dumps(state, indent=2, ensure_ascii=False))


def get_text(el, path, default=""):
    """Safe XML text extraction."""
    node = el.find(path)
    if node is not None and node.text:
        return node.text.strip()
    return default


def get_abstract(article):
    """Extract full abstract text from PubMed XML."""
    abstract_el = article.find(".//Abstract")
    if abstract_el is None:
        return ""
    parts = []
    for text_el in abstract_el.findall("AbstractText"):
        label = text_el.get("Label", "")
        # AbstractText may contain mixed content
        text = ET.tostring(text_el, encoding="unicode", method="text").strip()
        if label:
            parts.append(f"**{label}**: {text}")
        else:
            parts.append(text)
    return "\n\n".join(parts)


def get_authors(article):
    """Extract author names."""
    authors = []
    for author in article.findall(".//Author"):
        last = get_text(author, "LastName")
        first = get_text(author, "ForeName")
        if last:
            authors.append(f"{first} {last}".strip())
    return authors


def get_pmc_id(article_wrap) -> str:
    """Extract PMC ID if available."""
    for aid in article_wrap.findall(".//ArticleId"):
        if aid.get("IdType") == "pmc" and aid.text:
            return aid.text.strip()
    return ""


def fetch_pmc_fulltext(pmc_id: str) -> str:
    """Fetch full text from PMC as structured Markdown."""
    api_key_param = f"&api_key={API_KEY}" if API_KEY else ""
    url = f"{PMC_EFETCH}?db=pmc&id={pmc_id}&retmode=xml{api_key_param}"
    xml_data = api_get(url, as_xml=True)
    if not xml_data:
        return ""

    try:
        root = ET.fromstring(xml_data)
    except ET.ParseError:
        return ""

    sections = []
    body = root.find(".//body")
    if body is None:
        return ""

    for sec in body.findall(".//sec"):
        title_el = sec.find("title")
        sec_title = title_el.text.strip() if title_el is not None and title_el.text else ""

        paragraphs = []
        for p in sec.findall("p"):
            text = ET.tostring(p, encoding="unicode", method="text").strip()
            if text:
                paragraphs.append(text)

        if sec_title and paragraphs:
            sections.append(f"### {sec_title}\n\n" + "\n\n".join(paragraphs))
        elif paragraphs:
            sections.append("\n\n".join(paragraphs))

    return "\n\n".join(sections)


def get_doi(article):
    """Extract DOI from ArticleIdList."""
    for aid in article.findall(".//ArticleId"):
        if aid.get("IdType") == "doi" and aid.text:
            return aid.text.strip()
    # fallback: ELocationID
    for eloc in article.findall(".//ELocationID"):
        if eloc.get("EIdType") == "doi" and eloc.text:
            return eloc.text.strip()
    return ""


def fetch_keyword(keyword, tag, state):
    """Fetch papers for one keyword."""
    print(f"\n🔍 Searching: {keyword}")
    api_key_param = f"&api_key={API_KEY}" if API_KEY else ""

    # Step 1: esearch
    url = f"{ESEARCH}?db=pubmed&term={quote(keyword)}+AND+{DATE_RANGE}&retmax={PER_KEYWORD}&retmode=json&sort=relevance{api_key_param}"
    data = api_get(url)
    if not data:
        print("  ❌ esearch failed")
        return []

    result = data.get("esearchresult", {})
    total = result.get("count", "0")
    pmids = result.get("idlist", [])
    print(f"  Found {total} total, got {len(pmids)} PMIDs")

    if not pmids:
        return []

    # Filter already fetched
    new_pmids = [p for p in pmids if p not in state["fetched_pmids"]]
    if not new_pmids:
        print("  All already fetched")
        return []

    time.sleep(1)

    # Step 2: efetch XML for full details
    ids_str = ",".join(new_pmids)
    url = f"{EFETCH}?db=pubmed&id={ids_str}&retmode=xml{api_key_param}"
    xml_data = api_get(url, as_xml=True)
    if not xml_data:
        print("  ❌ efetch failed")
        return []

    root = ET.fromstring(xml_data)
    new_papers = []

    for article_wrap in root.findall(".//PubmedArticle"):
        article = article_wrap.find("MedlineCitation/Article")
        if article is None:
            continue

        pmid = get_text(article_wrap, "MedlineCitation/PMID")
        title = get_text(article, "ArticleTitle", "Untitled")
        # Clean HTML tags from title
        title = re.sub(r'<[^>]+>', '', title)
        abstract = get_abstract(article)
        authors = get_authors(article)
        doi = get_doi(article_wrap)
        journal = get_text(article, "Journal/Title", "")
        journal_short = get_text(article, "Journal/ISOAbbreviation", "")

        # Year
        pub_date = article.find("Journal/JournalIssue/PubDate")
        year = get_text(pub_date, "Year", "") if pub_date is not None else ""
        if not year:
            medline_date = get_text(pub_date, "MedlineDate", "") if pub_date is not None else ""
            year = medline_date[:4] if medline_date else "2025"

        # MeSH terms
        mesh_terms = []
        for mesh in article_wrap.findall(".//MeshHeading/DescriptorName"):
            if mesh.text:
                mesh_terms.append(mesh.text)

        # PMC full text
        pmc_id = get_pmc_id(article_wrap)
        fulltext = ""
        has_fulltext = False
        if pmc_id:
            print(f"    📄 Fetching PMC full text ({pmc_id})...")
            time.sleep(1)
            fulltext = fetch_pmc_fulltext(pmc_id)
            has_fulltext = bool(fulltext)
            if has_fulltext:
                print(f"    ✅ Full text: {len(fulltext)} chars")
            else:
                print(f"    ⚠️ No full text body in PMC")

        # Write markdown
        filename = f"{year}-pmid{pmid}-{slug(title)}.md"
        filepath = PAPERS_DIR / filename

        fulltext_section = ""
        if has_fulltext:
            fulltext_section = f"\n## Full Text\n\n{fulltext}\n"

        content = f"""---
title: "{title.replace('"', "'")}"
source: pubmed
keyword: "{keyword}"
pmid: "{pmid}"
pmc_id: "{pmc_id}"
year: {year}
authors: [{', '.join(f'"{a}"' for a in authors[:5])}]
doi: "{doi}"
journal: "{journal}"
journal_short: "{journal_short}"
tags: [{tag}]
mesh_terms: [{', '.join(f'"{m}"' for m in mesh_terms[:10])}]
has_fulltext: {str(has_fulltext).lower()}
content_layer: L1
fetched: "{datetime.now().strftime('%Y-%m-%d')}"
---

## Abstract

{abstract if abstract else '_No abstract available._'}
{fulltext_section}
## Metadata

- **PMID**: [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)
- **PMC**: {f'[{pmc_id}](https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/)' if pmc_id else 'N/A (not in PMC)'}
- **DOI**: {f'https://doi.org/{doi}' if doi else 'N/A'}
- **Journal**: {journal} ({journal_short})
- **Authors**: {', '.join(authors)}
- **MeSH Terms**: {', '.join(mesh_terms) if mesh_terms else 'N/A'}
- **Full Text**: {'✅ Included below' if has_fulltext else '❌ Abstract only'}
"""
        filepath.write_text(content, encoding="utf-8")
        state["fetched_pmids"].append(pmid)
        new_papers.append({"title": title, "keyword": keyword, "pmid": pmid, "file": filename})
        print(f"  ✅ {title[:60]}...")

    return new_papers


def main():
    PAPERS_DIR.mkdir(exist_ok=True)
    state = load_state()
    all_new = []

    for keyword, tag in KEYWORDS:
        new = fetch_keyword(keyword, tag, state)
        all_new.extend(new)
        time.sleep(2)  # respect rate limit

    save_state(state)
    print(f"\n📊 Summary: {len(all_new)} new papers fetched")
    for p in all_new:
        print(f"  [PMID {p['pmid']}] {p['title'][:60]}")

    return len(all_new)


if __name__ == "__main__":
    n = main()
    sys.exit(0)
