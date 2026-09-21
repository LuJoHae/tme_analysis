---
name: zotero-citation-manager
description: Skill for searching the local Zotero library, retrieving bibliographic metadata, automatically formatting BibTeX entries into article/references.bib, and citing them in Typst.
---

# Zotero Citation Manager Skill

Use this skill when you need to find literature in the user's local Zotero database, add bibliographic citations to the article's bibliography (`article/references.bib`), or cite papers in Typst (`@citationKey`).

## 1. Finding Literature in Zotero

The Zotero MCP server (`zotero`) provides tools to query the local library:

### A. General Keyword / Title / Author Search
Use the `search_library` MCP tool via `call_mcp_tool`:
```json
{
  "ServerName": "zotero",
  "ToolName": "search_library",
  "Arguments": {
    "q": "immunotherapy CD8 T cells",
    "limit": 5
  }
}
```
This returns a list of candidate items with their `key`, `title`, `creators`, and `date`.

### B. Full-Text Search
To search inside PDF attachments or notes, use `search_fulltext`:
```json
{
  "ServerName": "zotero",
  "ToolName": "search_fulltext",
  "Arguments": {
    "q": "TGM6 melanoma",
    "limit": 5
  }
}
```

---

## 2. Retrieving Full Item Metadata

Once an `itemKey` is identified, retrieve its detailed bibliographic record:
```json
{
  "ServerName": "zotero",
  "ToolName": "get_item_details",
  "Arguments": {
    "itemKey": "CKV6ST8A"
  }
}
```
The returned JSON structure contains:
- `key` (e.g. `"CKV6ST8A"`)
- `citationKey` (if Better BibTeX is active, e.g. `"vanderleunCD8CellStates2020"`)
- `title`, `creators` (authors list)
- `publicationTitle`, `date`, `volume`, `issue`, `pages`
- `DOI`, `url`

---

## 3. Adding the Reference to `article/references.bib`

Run the functional helper script `.agents/skills/zotero-citation-manager/scripts/add_zotero_reference.py` using `run_command`.

### Method: Passing JSON via CLI or stdin
```bash
python3 .agents/skills/zotero-citation-manager/scripts/add_zotero_reference.py \
  --bib-path article/references.bib \
  --json '<zotero_item_details_json>'
```

Or pass via pipe:
```bash
echo '<zotero_item_details_json>' | python3 .agents/skills/zotero-citation-manager/scripts/add_zotero_reference.py --bib-path article/references.bib
```

### Output Format
The script prints a structured JSON response:
```json
{
  "status": "ADDED: Successfully added @vanderleunCD8CellStates2020 to article/references.bib",
  "citekey": "vanderleunCD8CellStates2020",
  "typst_citation": "@vanderleunCD8CellStates2020",
  "bibtex": "@article{vanderleunCD8CellStates2020,\n  title = {CD8+ T cell states in human cancer: insights from single-cell analysis},\n  author = {van der Leun, Anne M. and Thommen, Daniela S. and Schumacher, Ton N.},\n  journal = {Nature Reviews Cancer},\n  year = {2020},\n  volume = {20},\n  number = {4},\n  pages = {218-232},\n  doi = {10.1038/s41568-019-0235-4}\n}\n"
}
```

### Built-in Idempotency & Duplicate Prevention
- If the citation key or the DOI already exists in `references.bib`, the script will **not** duplicate the entry.
- It returns status `EXISTS: ...` along with the existing `citekey` so you can immediately use it in the manuscript text.

---

## 4. Citing References in Typst

In Typst, citations use the `@` symbol followed by the citation key:
- Simple citation: `@vanderleunCD8CellStates2020`
- Multiple citations: `[@vanderleunCD8CellStates2020; @hugo2016genomic]`
- Integrated in text: `As demonstrated by @vanderleunCD8CellStates2020, ...`

Typst will automatically compile citations according to the bibliography style configured in `article/main.typ` (`#bibliography("references.bib", style: "nature")`).

---

## 5. End-to-End Workflow Example for Agents

1. **User asks**: "Can you cite the 2020 Nature Reviews Cancer paper on CD8 T cell states from my Zotero library?"
2. **Step 1**: Search library:
   `call_mcp_tool(ServerName="zotero", ToolName="search_library", Arguments={"q": "CD8 T cell states human cancer", "limit": 3})`
3. **Step 2**: Inspect details of matching result:
   `call_mcp_tool(ServerName="zotero", ToolName="get_item_details", Arguments={"itemKey": "CKV6ST8A"})`
4. **Step 3**: Append to bibliography:
   Run `add_zotero_reference.py` with the returned payload.
5. **Step 4**: Update Typst source:
   Insert `@vanderleunCD8CellStates2020` into the appropriate chapter file (e.g., `article/chapters/07_results_tcell_myeloid.typ`).
6. **Step 5**: Recompile article:
   Run `make -C article compile-article`.
