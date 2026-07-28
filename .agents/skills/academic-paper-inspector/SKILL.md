---
name: academic-paper-inspector
description: Skill for extracting key scientific findings and metadata from academic papers, particularly those behind publisher paywalls (e.g. Cell, Nature, Science).
---

# Academic Paper Inspector Skill

When the user asks you to inspect, read, or summarize an academic paper (especially from a major publisher URL), follow these instructions strictly to avoid paywalls and anti-bot measures.

## 1. Avoid Raw Scraping
**DO NOT** use standard HTTP request tools like `read_url_content`, `curl`, or Python `requests` directly against publisher URLs (e.g., `cell.com`, `nature.com`). You will almost always receive a `403 Forbidden` error because these sites use Cloudflare and other anti-bot protections to enforce their paywalls.

## 2. Use Semantic Web Search
Instead, rely exclusively on your native **`search_web`** tool. The semantic search index contains comprehensive metadata and abstracts for almost all published academic literature.

**How to query:**
- Pass the full URL as the query (e.g., `query="https://www.cell.com/cell/fulltext/S0092-8674(18)31394-1"`).
- If the URL search returns irrelevant results (due to identical publisher PII formats), extract the title of the paper or its DOI from the URL, and search for that instead (e.g. `query="Defining T cell states associated with response to checkpoint immunotherapy in melanoma"`).

## 3. Structured Extraction
When reporting back to the user, always structure your summary to include:
- **Paper Details**: Journal, Publication Date, DOI, and Authors.
- **Methodology**: Experimental techniques (e.g., scRNA-seq, CRISPR screens, bulk RNA-seq).
- **Key Findings**: The primary biological discoveries, mechanisms, or predictive biomarkers identified in the study.
- **Data Availability**: Look closely for accession numbers for public data repositories (e.g., GEO accession `GSE...`, SRA, dbGaP) so the user can download the raw data if needed.
