import json
import html
import re
import requests
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional
from urllib.parse import quote

from mcp.server.fastmcp import FastMCP

PROJECT_ROOT = Path(__file__).resolve().parents[2]

mcp = FastMCP("ReportGenerationTools")

RAW_REPORTS_PATH = str(PROJECT_ROOT / "generated_artifacts" / "reports" / "raw_reports")
HTML_REPORTS_PATH = str(PROJECT_ROOT / "generated_artifacts" / "reports" / "html_reports")


# ---------------------------
# Helpers (internal only)
# ---------------------------

def _h(x: Any) -> str:
    """HTML-escape for safe insertion."""
    if x is None:
        return ""
    return html.escape(str(x), quote=True)


def _ref_key(ref: Dict[str, Any]) -> tuple:
    """Stable key for de-duplication and numbering.

    Prefer ref_id (because VERA often reuses IDs across sections), then URL, then title.
    """
    if not isinstance(ref, dict):
        return ("__invalid__", "")
    rid = str(ref.get("ref_id") or "").strip()
    if rid:
        return ("id", rid.lower())
    url = str(ref.get("url") or "").strip()
    if url:
        return ("url", url)
    title = str(ref.get("title") or "").strip()
    if title:
        return ("title", title.lower())
    return ("__empty__", "")


def _build_citation_index(
    references: List[Dict[str, Any]],
    ref_map: Dict[tuple, int],
) -> tuple[List[Dict[str, Any]], Dict[int, str]]:
    """Build patterns to turn inline reference IDs into numbered citations."""
    patterns: List[Dict[str, Any]] = []
    num_to_title: Dict[int, str] = {}
    seen_pat = set()

    def add(pattern: str, ref_num: int) -> None:
        if not pattern or not ref_num:
            return
        key = (pattern, int(ref_num))
        if key in seen_pat:
            return
        seen_pat.add(key)
        patterns.append({"pattern": pattern, "ref_num": int(ref_num)})

    for ref in references:
        if not isinstance(ref, dict):
            continue
        ref_num = _get_reference_number(ref, ref_map)
        if not ref_num:
            continue

        title = str(ref.get("title") or "").strip()
        if title and ref_num not in num_to_title:
            num_to_title[ref_num] = title

        rid = str(ref.get("ref_id") or "").strip()
        url = str(ref.get("url") or "").strip()

        # Exact ref_id (case-insensitive)
        if rid and len(rid) >= 4:
            add(r"(?i)\b" + re.escape(rid) + r"\b", ref_num)

        # PMID variants
        if rid:
            m = re.search(r"PMID[:\s]*([0-9]{4,10})", rid, flags=re.I)
            if m:
                pmid = m.group(1)
                add(r"(?i)PMID[:\s]*" + re.escape(pmid) + r"\b", ref_num)
        
        # PMC variants
        if rid:
            m = re.search(r"PMC([0-9]{4,10})", rid, flags=re.I)
            if m:
                pmc = m.group(1)
                add(r"(?i)PMC[:\s]*" + re.escape(pmc) + r"\b", ref_num)
                add(r"(?i)PMC" + re.escape(pmc) + r"\b", ref_num)

        # DOI variants
        doi = str(ref.get("doi") or "").strip()
        if not doi:
            if rid:
                m2 = re.search(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", rid, flags=re.I)
                if m2:
                    doi = m2.group(0)
            if (not doi) and url:
                m3 = re.search(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", url, flags=re.I)
                if m3:
                    doi = m3.group(0)
        if doi and len(doi) >= 8:
            add(r"(?i)\b" + re.escape(doi) + r"\b", ref_num)
            add(r"(?i)doi[:\s]*" + re.escape(doi) + r"\b", ref_num)

        # Exact URL match
        if url and len(url) >= 12:
            add(re.escape(url), ref_num)

    patterns.sort(key=lambda d: len(d.get("pattern", "")), reverse=True)
    return patterns, num_to_title


def _strip_citation_patterns(
    text: Any,
    cite_patterns: List[Dict[str, Any]],
) -> str:
    """Remove recognized reference IDs from text (for headings where citations shouldn't appear)."""
    if text is None:
        return ""
    s = str(text).strip()
    if not s:
        return ""

    out = s
    for item in cite_patterns:
        pat = item.get("pattern")
        if not pat:
            continue
        try:
            out = re.sub(pat, "", out)
        except re.error:
            continue

    out = re.sub(r'\[\s*(?:,\s*)*\]', '', out)
    out = re.sub(r'\[\s*[^\]]*\]', '', out)
    out = re.sub(r'\s*,\s*,', ',', out)
    out = re.sub(r'^\s*,\s*', '', out)
    out = re.sub(r'\s*,\s*$', '', out)
    out = re.sub(r'\s{2,}', ' ', out)
    return out.strip()


def _render_text_with_cites(
    text: Any,
    cite_patterns: List[Dict[str, Any]],
    num_to_title: Optional[Dict[int, str]] = None,
) -> str:
    """Escape text and replace any recognized reference IDs with numbered [n] links."""
    if text is None:
        return ""
    s = str(text)
    if not s.strip():
        return ""

    out = s
    for item in cite_patterns:
        pat = item.get("pattern")
        n = item.get("ref_num")
        if not pat or not n:
            continue
        placeholder = f"__CITE_{int(n)}__"
        try:
            out = re.sub(pat, placeholder, out)
        except re.error:
            continue

    esc = html.escape(out, quote=True)

    def replace_citations(m: re.Match) -> str:
        full_match = m.group(0)
        has_brackets = full_match.startswith('[') and full_match.endswith(']')
        inner = full_match[1:-1] if has_brackets else full_match
        nums = re.findall(r'__CITE_(\d+)__', inner)
        if not nums:
            return full_match
        ref_nums = sorted(set(int(n) for n in nums))
        links = []
        for n in ref_nums:
            title = (num_to_title or {}).get(n, "")
            title_attr = f' title="{_h(title)}"' if title else ""
            links.append(f'<a href="#ref-{n}" class="ref-link"{title_attr}>{n}</a>')
        result = ', '.join(links)
        return '[' + result + ']' if not has_brackets else '[' + result + ']'

    result = re.sub(r'\[?__CITE_\d+__(?:[\s,;]*__CITE_\d+__)*\]?', replace_citations, esc)
    result = re.sub(r'\[\s*(?:,\s*)*\]', '', result)
    return result


def _as_list(x: Any) -> List[Any]:
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


def _extract_int(s: Any) -> Optional[int]:
    if s is None:
        return None
    m = re.search(r"(\d+)", str(s))
    return int(m.group(1)) if m else None


def _pct(x: Any, digits: int = 2) -> str:
    """0–1 -> 0–100 percent string. If x is not numeric, return empty."""
    try:
        v = float(x)
        if v > 1:
            return f"{v:.{digits}f}%"
        v *= 100.0
        return f"{v:.{digits}f}%"
    except Exception:
        return ""


def _signed_pct(x: Any, digits: int = 2) -> str:
    """Like _pct but with explicit +/− sign prefix."""
    try:
        v = float(x)
        raw = v if abs(v) > 1 else v * 100.0
        sign = "+" if raw > 0 else "\u2212" if raw < 0 else ""
        return f"{sign}{abs(raw):.{digits}f}%"
    except Exception:
        return ""


def _format_scientific(v: float) -> str:
    """Format a float as HTML scientific notation (e.g. 3.34 × 10⁻³)."""
    import math
    if v == 0:
        return "0"
    exp = int(math.floor(math.log10(abs(v))))
    mantissa = v / (10 ** exp)
    if exp < 0:
        return f"{abs(mantissa):.2f} &times; 10<sup>&minus;{abs(exp)}</sup>"
    elif exp == 0:
        return f"{v:.3f}"
    else:
        return f"{abs(mantissa):.2f} &times; 10<sup>{exp}</sup>"


def _format_ge_cell(col: str, value: Any) -> str:
    """Format a gene-expression result cell for display."""
    if value is None or str(value).strip() == "":
        return ""
    try:
        v = float(value)
    except (ValueError, TypeError):
        return _h(str(value))

    col_norm = col.lower().replace("-", "").replace("_", "").replace(" ", "")

    if any(tok in col_norm for tok in ["pvalue", "qvalue", "pval", "qval"]):
        return _format_scientific(v)

    if any(tok in col_norm for tok in ["logfoldchange", "logfc", "log2fc"]):
        sign = "+" if v > 0 else ("\u2212" if v < 0 else "")
        return _h(f"{sign}{abs(v):.3f}")

    if any(tok in col_norm for tok in ["tstat", "tstatistic"]):
        sign = "\u2212" if v < 0 else ""
        return _h(f"{sign}{abs(v):.1f}")

    if any(tok in col_norm for tok in ["effectsize", "cohensd", "cohens"]):
        sign = "\u2212" if v < 0 else ""
        return _h(f"{sign}{abs(v):.3f}")

    if abs(v) < 0.001 and v != 0:
        return _format_scientific(v)
    return _h(f"{v:.3f}")


def _yesno(x: Any) -> str:
    if x is True:
        return "Yes"
    if x is False:
        return "No"
    return ""


def _direction_from_delta(delta: Any) -> str:
    try:
        d = float(delta)
    except Exception:
        return ""
    if d > 0:
        return "gain of phosphosite"
    if d < 0:
        return "loss of phosphosite"
    return "neutral"


def _clean_site_label(label: str, center: str, center_pos: int) -> str:
    """Normalize site labels, avoiding placeholders like 'unknown_site'."""
    raw = (label or "").strip()
    low = raw.lower()
    if not raw or low in ("unknown_site", "unknown", "site"):
        if center and center_pos:
            return f"{center}{center_pos}"
        if center:
            return center
        if center_pos:
            return str(center_pos)
        return "Site"
    return raw


def _extract_pvalue(result: Dict[str, Any]) -> Optional[float]:
    """Find a p-value field in a result dict, if present."""
    if not isinstance(result, dict):
        return None
    for k, v in result.items():
        key = str(k).lower().replace("-", "").replace("_", "").replace(" ", "")
        if "pvalue" in key or key in ("pval", "p"):
            try:
                return float(v)
            except Exception:
                return None
    return None


def _render_table(headers: List[str], rows: List[List[str]]) -> str:
    th = "".join(f"<th>{_h(h)}</th>" for h in headers)
    trs = []
    for r in rows:
        tds = "".join(f"<td>{c}</td>" for c in r)  # cells may already be escaped
        trs.append(f"<tr>{tds}</tr>")
    tbody = "".join(trs)
    return (
        '<div class="tablewrap">'
        '<table class="table">'
        f"<thead><tr>{th}</tr></thead>"
        f"<tbody>{tbody}</tbody>"
        "</table>"
        "</div>"
    )


def _build_positions_ruler(overall_start: int, mutation_pos: int, overall_end: int) -> str:
    slots_count = overall_end - overall_start + 1
    width = slots_count * 5
    buf = [" "] * width

    def place_num(n: int, slot_index: int) -> None:
        s = str(n)
        off = max(0, min(width - 1, slot_index * 5))
        for i, ch in enumerate(s):
            if off + i < width:
                buf[off + i] = ch

    place_num(overall_start, 0)
    place_num(mutation_pos, mutation_pos - overall_start)
    place_num(overall_end, slots_count - 1)
    return "".join(buf)


def _slotify_sequence(
    seq: str,
    seq_start: int,
    overall_start: int,
    overall_end: int,
) -> List[str]:
    slots_count = overall_end - overall_start + 1
    slots = [" " * 5 for _ in range(slots_count)]
    for i, res in enumerate(seq):
        pos = seq_start + i
        if overall_start <= pos <= overall_end:
            slots[pos - overall_start] = f"{res}{' ' * 4}"
    return slots


def _highlight_slot(slots: List[str], mutation_pos: int, overall_start: int) -> None:
    idx = mutation_pos - overall_start
    if idx < 0 or idx >= len(slots):
        return
    raw = slots[idx]
    res = raw[0] if raw else " "
    slots[idx] = f'<span class="mut">{_h(res)}</span>' + (" " * 4)


def _render_window_codeblock(
    label: str,
    center_pos: int,
    mutation_pos: int,
    wt_15mer: str,
    mut_15mer: str,
) -> str:
    wt_start = center_pos - 7
    wt_end = center_pos + 7
    mut_start = center_pos - 7
    mut_end = center_pos + 7

    overall_start = min(wt_start, mut_start, mutation_pos)
    overall_end = max(wt_end, mut_end, mutation_pos)

    ruler = _build_positions_ruler(overall_start, mutation_pos, overall_end)

    wt_slots = _slotify_sequence(wt_15mer, wt_start, overall_start, overall_end)
    mut_slots = _slotify_sequence(mut_15mer, mut_start, overall_start, overall_end)

    _highlight_slot(wt_slots, mutation_pos, overall_start)
    _highlight_slot(mut_slots, mutation_pos, overall_start)

    wt_line = "".join(wt_slots)
    mut_line = "".join(mut_slots)

    pre = (
        f"Positions:  {ruler}\n"
        f"Wildtype:   {wt_line}\n"
        f"Mutant:     {mut_line}"
    )

    return (
        '<figure class="card">'
        f'  <figcaption><strong>{_h(label)}</strong></figcaption>'
        f'  <pre class="codeblock">\n{pre}\n  </pre>'
        "</figure>"
    )


def _merge_references(raw: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Merge references from multiple places; de-duplicate primarily by ref_id."""
    seen = set()
    out: List[Dict[str, Any]] = []

    for src in _as_list(raw.get("collector_output", {}).get("literature_sources")):
        if not isinstance(src, dict):
            continue
        key = _ref_key(src)
        if key in seen:
            continue
        seen.add(key)
        out.append(src)

    for hyp in _as_list(raw.get("impact_hypotheses_output", {}).get("hypotheses")):
        if not isinstance(hyp, dict):
            continue
        for src in _as_list(hyp.get("supporting_references")):
            if not isinstance(src, dict):
                continue
            key = _ref_key(src)
            if key in seen:
                continue
            seen.add(key)
            out.append(src)

    return out


def _create_reference_map(references: List[Dict[str, Any]]) -> Dict[tuple, int]:
    """Create a mapping from a stable reference key to reference numbers (1-indexed)."""
    ref_map: Dict[tuple, int] = {}
    for idx, ref in enumerate(references, start=1):
        if not isinstance(ref, dict):
            continue
        key = _ref_key(ref)
        if key not in ref_map:
            ref_map[key] = idx
    return ref_map


def _get_reference_number(ref: Dict[str, Any], ref_map: Dict[tuple, int]) -> Optional[int]:
    """Get the reference number for a reference dict."""
    if not isinstance(ref, dict):
        return None
    return ref_map.get(_ref_key(ref))


def _resolve_uniprot_from_gene(gene_name: str, organism: str = "9606") -> str:
    """Query UniProt REST API to resolve a gene name to its primary accession."""
    if not gene_name:
        return ""
    query = f"gene:{gene_name}+AND+organism_id:{organism}+AND+reviewed:true"
    url = (
        "https://rest.uniprot.org/uniprotkb/search"
        f"?query={query}&fields=accession&format=json"
    )
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            results = data.get("results")
            if results:
                acc = results[0].get("primaryAccession")
                if acc:
                    return acc
    except Exception:
        pass
    return ""


def _infer_uniprot_accession(raw: Dict[str, Any]) -> str:
    """Best-effort UniProt accession extraction.

    Checks (in order):
      1. Direct fields in inputs / collector_output (uniprot_id, uniprot_accession)
      2. Reference ref_ids and URLs containing a UniProt accession
      3. Collector notes mentioning a UniProt ID
      4. UniProt REST API lookup using the gene name (always available in inputs)
    """
    inputs = raw.get("inputs", {}) or {}
    collector = raw.get("collector_output", {}) or {}

    for src in (inputs, collector):
        for key in ("uniprot_id", "uniprot_accession", "uniprot"):
            v = src.get(key)
            if isinstance(v, str) and re.fullmatch(r"[A-Z0-9]{6,10}", v.strip()):
                return v.strip()

    refs = _merge_references(raw)
    for r in refs:
        if not isinstance(r, dict):
            continue
        ref_id = str(r.get("ref_id", "") or "")
        m = re.search(r"UniProt([A-Z0-9]{6,10})", ref_id)
        if m:
            return m.group(1)
        url = str(r.get("url", "") or "")
        m = re.search(r"/uniprot(?:kb)?/([A-Z0-9]{6,10})(?:/|$)", url)
        if m:
            return m.group(1)

    notes = str(collector.get("notes", "") or "")
    m = re.search(r"UniProt\s+ID\s+([A-Z0-9]{6,10})", notes, re.I)
    if m:
        return m.group(1)
    m = re.search(r"UniProt[:\s]+([A-Z0-9]{6,10})", notes, re.I)
    if m:
        return m.group(1)

    gene_name = (inputs.get("gene_name") or "").strip()
    if gene_name:
        acc = _resolve_uniprot_from_gene(gene_name)
        if acc:
            return acc

    return ""


def _try_fallback_urls(uniprot_id: str) -> Optional[str]:
    """
    Fallback function to try alternative methods for finding AlphaFold structure URLs.
    Only called if the primary API endpoint fails.
    """
    af_model_id = f"AF-{uniprot_id}-F1"

    model_versions = ["v6", "v5", "v4", "v3", "v2"]
    for version in model_versions:
        direct_url = f"https://alphafold.ebi.ac.uk/files/{af_model_id}-model_{version}.pdb"
        try:
            print(f"[DEBUG] Trying direct URL check: {direct_url}")
            test_resp = requests.head(direct_url, timeout=10, allow_redirects=True)
            if test_resp.status_code == 200:
                return direct_url
        except Exception:
            continue

    entry_page_url = f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
    try:
        print(f"[DEBUG] Trying entry page scrape: {entry_page_url}")
        resp = requests.get(entry_page_url, timeout=30)
        if resp.status_code == 200:
            content = resp.text
            pdb_patterns = [
                rf'https://alphafold\.ebi\.ac\.uk/files/{af_model_id}[^"\'<>\s]*\.pdb',
                rf'https://alphafold\.ebi\.ac\.uk/files/AF-{uniprot_id}-F1[^"\'<>\s]*\.pdb',
                rf'/files/{af_model_id}[^"\'<>\s]*\.pdb',
            ]
            for pattern in pdb_patterns:
                matches = re.findall(pattern, content)
                if matches:
                    url = matches[0]
                    if url.startswith('/'):
                        url = 'https://alphafold.ebi.ac.uk' + url
                    return url
    except Exception:
        pass

    return None


def _get_pdb_url_from_uniprot(uniprot_id: str) -> Optional[str]:
    """
    Get AlphaFold PDB URL for a UniProt ID using the API, with fallback mechanisms.
    Returns the PDB URL string or None if not found.
    """
    print(f"[DEBUG] Fetching PDB URL for UniProt ID: {uniprot_id}")
    pdb_url = None

    meta_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        print(f"[DEBUG] Making API request to: {meta_url}")
        r = requests.get(meta_url, timeout=30)
        if r.status_code == 200:
            info = r.json()
            if info and len(info) > 0:
                pdb_url = info[0].get("pdbUrl")
    except Exception:
        pass

    if not pdb_url:
        af_model_id = f"AF-{uniprot_id}-F1"
        meta_url = f"https://alphafold.ebi.ac.uk/api/prediction/{af_model_id}"
        try:
            print(f"[DEBUG] Trying fallback API request to: {meta_url}")
            r = requests.get(meta_url, timeout=30)
            if r.status_code == 200:
                info = r.json()
                if info and len(info) > 0:
                    pdb_url = info[0].get("pdbUrl")
        except Exception:
            pass

    if not pdb_url:
        print(f"[DEBUG] Trying fallback URL methods for: {uniprot_id}")
        pdb_url = _try_fallback_urls(uniprot_id)

    print(f"[DEBUG] PDB URL fetch completed for {uniprot_id}: {'Found' if pdb_url else 'Not found'}")
    return pdb_url


def _infer_pdb_url(raw: Dict[str, Any], gene: str) -> str:
    """Choose a PDB URL for structure usage. Uses raw overrides, then UniProt→AlphaFold API."""
    inputs = raw.get("inputs", {}) or {}
    collector = raw.get("collector_output", {}) or {}

    for key in [
        "protein_pdb_url",
        "pdb_url",
        "alphafold_pdb_url",
        "structure_pdb_url",
    ]:
        v = inputs.get(key) or collector.get(key)
        if isinstance(v, str) and v.strip():
            return v.strip()

    uniprot = _infer_uniprot_accession(raw)
    if uniprot:
        pdb_url = _get_pdb_url_from_uniprot(uniprot)
        if pdb_url:
            return pdb_url

    return ""


def _norm_site_key(wt_15: str, mut_15: str, center_pos: int) -> str:
    return f"{(wt_15 or '').strip()}|{(mut_15 or '').strip()}|{int(center_pos or 0)}"


def _is_better_label(label1: str, label2: str) -> bool:
    if not label1 or label1.strip().lower() in ("", "unknown_site", "site"):
        return False
    if not label2 or label2.strip().lower() in ("", "unknown_site", "site"):
        return True
    return False


def _collect_impacted_sites(raw: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Build impacted sites from:
      1) phosphosite windows that include kinase_specificity_changes
      2) hypotheses' implicated_impacted_sites

    Deduplicates based on sequence (wt_15mer, mut_15mer, center_pos).
    """
    out: List[Dict[str, Any]] = []
    seen: Dict[str, Dict[str, Any]] = {}

    windows = _as_list(raw.get("collector_output", {}).get("phosphosite_windows"))
    for w in windows:
        if not isinstance(w, dict):
            continue
        ksc = _as_list(w.get("kinase_specificity_changes"))
        if not ksc:
            continue

        wt_15 = (w.get("wildtype_15mer") or "").strip()
        mut_15 = (w.get("mutant_15mer") or "").strip()
        center_pos = _extract_int(w.get("evidence_info")) or _extract_int(w.get("site_label")) or 0
        center = (w.get("center") or "").strip()
        site_label = _clean_site_label(w.get("site_label") or "", center, center_pos)

        key = _norm_site_key(wt_15, mut_15, center_pos)

        events = []
        for e in ksc:
            if not isinstance(e, dict):
                continue
            events.append(
                {
                    "kinase": e.get("kinase", ""),
                    "wildtype_confidence": e.get("wildtype_confidence", None),
                    "mutant_confidence": e.get("mutant_confidence", None),
                    "delta_metric": e.get("delta_metric", None),
                    "direction": _direction_from_delta(e.get("delta_metric", None)),
                }
            )

        site_dict = {
            "site_label": site_label,
            "center_pos": center_pos,
            "wildtype_15mer": wt_15,
            "mutant_15mer": mut_15,
            "evidence_info": w.get("evidence_info", ""),
            "site_notes": w.get("window_notes", ""),
            "kinase_events": events,
        }

        if key in seen:
            existing = seen[key]
            if _is_better_label(site_label, existing.get("site_label", "")):
                existing["site_label"] = site_label
            if w.get("evidence_info") and not existing.get("evidence_info"):
                existing["evidence_info"] = w.get("evidence_info", "")
            if w.get("window_notes") and not existing.get("site_notes"):
                existing["site_notes"] = w.get("window_notes", "")
            existing_kinases = {e.get("kinase", ""): e for e in existing["kinase_events"]}
            for evt in events:
                kinase_name = evt.get("kinase", "")
                if kinase_name not in existing_kinases:
                    existing["kinase_events"].append(evt)
                    existing_kinases[kinase_name] = evt
        else:
            seen[key] = site_dict
            out.append(site_dict)

    hyps = _as_list(raw.get("impact_hypotheses_output", {}).get("hypotheses"))
    for hyp in hyps:
        if not isinstance(hyp, dict):
            continue
        for site in _as_list(hyp.get("implicated_impacted_sites")):
            if not isinstance(site, dict):
                continue
            wt_15 = (site.get("wildtype_15mer") or "").strip()
            mut_15 = (site.get("mutant_15mer") or "").strip()

            center = (site.get("center") or "").strip()
            center_pos = (
                _extract_int(site.get("center_pos"))
                or _extract_int(site.get("site_label"))
                or _extract_int(site.get("evidence_info"))
                or 0
            )
            raw_label = (site.get("site_label") or "").strip()
            label = _clean_site_label(raw_label, center, center_pos)

            key = _norm_site_key(wt_15, mut_15, center_pos)

            events = []
            for e in _as_list(site.get("kinase_events")):
                if not isinstance(e, dict):
                    continue
                events.append(
                    {
                        "kinase": e.get("kinase", ""),
                        "wildtype_confidence": e.get("wildtype_confidence", None),
                        "mutant_confidence": e.get("mutant_confidence", None),
                        "delta_metric": e.get("delta_metric", None),
                        "direction": e.get("direction", _direction_from_delta(e.get("delta_metric", None))),
                    }
                )

            if key in seen:
                existing = seen[key]
                if _is_better_label(label, existing.get("site_label", "")):
                    existing["site_label"] = label
                if site.get("evidence_info") and not existing.get("evidence_info"):
                    existing["evidence_info"] = site.get("evidence_info", "")
                if site.get("site_notes") and not existing.get("site_notes"):
                    existing["site_notes"] = site.get("site_notes", "")
                existing_kinases = {e.get("kinase", ""): e for e in existing["kinase_events"]}
                for evt in events:
                    kinase_name = evt.get("kinase", "")
                    if kinase_name not in existing_kinases:
                        existing["kinase_events"].append(evt)
                        existing_kinases[kinase_name] = evt
            else:
                site_dict = {
                    "site_label": label,
                    "center_pos": center_pos,
                    "wildtype_15mer": wt_15,
                    "mutant_15mer": mut_15,
                    "evidence_info": site.get("evidence_info", ""),
                    "site_notes": site.get("site_notes", ""),
                    "kinase_events": events,
                }
                seen[key] = site_dict
                out.append(site_dict)

    return out


def _render_kinase_events_table(events: List[Dict[str, Any]]) -> str:
    sortable: List[Dict[str, Any]] = []
    for e in events:
        if not isinstance(e, dict):
            continue
        sortable.append(e)

    def sort_key(evt: Dict[str, Any]) -> float:
        try:
            return abs(float(evt.get("delta_metric", 0.0)))
        except Exception:
            return 0.0

    sortable.sort(key=sort_key, reverse=True)

    rows = []
    for e in sortable:
        direction = (e.get("direction") or "").strip()
        dir_lower = direction.lower()
        if "gain" in dir_lower:
            dir_cell = f"<span style='color:#16a34a;font-weight:700'>{_h(direction)}</span>"
        elif "loss" in dir_lower:
            dir_cell = f"<span style='color:#dc2626;font-weight:700'>{_h(direction)}</span>"
        else:
            dir_cell = _h(direction)

        rows.append(
            [
                _h(e.get("kinase", "")),
                _h(_pct(e.get("wildtype_confidence"))),
                _h(_pct(e.get("mutant_confidence"))),
                _signed_pct(e.get("delta_metric")),
                dir_cell,
            ]
        )
    return _render_table(
        ["Kinase", "WT phosphorylation conf", "Mut phosphorylation conf", "\u0394 conf", "Direction"],
        rows or [[_h("\u2014"), _h(""), _h(""), _h(""), _h("")]],
    )


# ---------------------------
# NEW: Pipe-delimited interpretation renderer (unchanged from your version)
# ---------------------------

def _render_pipe_delimited_interpretation(
    text: str,
    cite_patterns: Optional[List[Dict[str, Any]]] = None,
    num_to_title: Optional[Dict[int, str]] = None,
) -> str:
    if not text:
        return ""

    sections = [s.strip() for s in text.split(" | ") if s.strip()]
    out: List[str] = []
    compact_headers = {"verdict", "confidence"}

    for sec in sections:
        if ":" in sec:
            header, body = sec.split(":", 1)
            header = header.strip()
            body = body.strip()

            if header.lower() in compact_headers:
                out.append(
                    f"<p><strong>{_h(header)}:</strong> "
                    f"{_render_text_with_cites(body, cite_patterns or [], num_to_title)}</p>"
                )
                continue

            clean_header = _strip_citation_patterns(header, cite_patterns or [])
            out.append(f"<h5>{_h(clean_header)}</h5>")

            if "•" in body:
                items = [x.strip() for x in body.split("•") if x.strip()]
                if items:
                    out.append(
                        "<ul>"
                        + "".join(
                            f"<li>{_render_text_with_cites(it, cite_patterns or [], num_to_title)}</li>"
                            for it in items
                        )
                        + "</ul>"
                    )
                else:
                    out.append(f"<p>{_render_text_with_cites(body, cite_patterns or [], num_to_title)}</p>")
            else:
                out.append(f"<p>{_render_text_with_cites(body, cite_patterns or [], num_to_title)}</p>")
        else:
            out.append(f"<p>{_render_text_with_cites(sec, cite_patterns or [], num_to_title)}</p>")

    return "".join(out)


def _render_markdownish(
    text: str,
    cite_patterns: Optional[List[Dict[str, Any]]] = None,
    num_to_title: Optional[Dict[int, str]] = None,
) -> str:
    if not text:
        return ""

    s = str(text).strip()
    if not s:
        return ""

    if ("\n" not in s) and (" | " in s):
        return _render_pipe_delimited_interpretation(s, cite_patterns=cite_patterns, num_to_title=num_to_title)

    lines = s.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    out: List[str] = []
    list_level = 0

    def close_lists(target_level: int = 0) -> None:
        nonlocal list_level
        while list_level > target_level:
            out.append("</ul>")
            list_level -= 1

    def open_lists(target_level: int) -> None:
        nonlocal list_level
        while list_level < target_level:
            out.append("<ul>")
            list_level += 1

    def looks_like_heading(line: str) -> bool:
        t = line.strip()
        if not t:
            return False
        if t.endswith(":"):
            return True
        if len(t) > 90:
            return False
        if any(tok in t for tok in ["≈", "->"]):
            return False
        if "." in t:
            return False
        return bool(re.match(r"^[A-Z][A-Za-z0-9 /–\-\(\),]+$", t))

    def looks_like_group_label(bullet_text: str) -> bool:
        t = bullet_text.strip()
        if not t:
            return False
        low = t.lower()

        if t.endswith(":"):
            return True
        if any(tok in t for tok in ["≈", "->"]):
            return False
        if re.search(r"\b(logfc|effect size|q\s*=?|p[_ ]value|q[_ ]value|t[_ ]stat|fdr)\b", low):
            return False
        if low.startswith("interpretation:"):
            return False

        if len(t) <= 70 and (t.count(" ") <= 6):
            return True
        return False

    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if not line.strip():
            close_lists(0)
            continue

        if looks_like_heading(line) and not re.match(r"^\s*(?:•|\-|\*)\s+", line):
            close_lists(0)
            t = line.strip().rstrip(":")
            t = _strip_citation_patterns(t, cite_patterns or [])
            out.append(f"<h4>{_h(t)}</h4>")
            continue

        m = re.match(r"^(\s*)(?:•|\-|\*)\s+(.*)$", line)
        if m:
            indent_spaces = len(m.group(1) or "")
            content = (m.group(2) or "").strip()

            target_level = 1 + max(0, indent_spaces // 2)

            if looks_like_group_label(content):
                close_lists(0)
                t = content.strip().rstrip(":")
                t = _strip_citation_patterns(t, cite_patterns or [])
                out.append(f"<h5>{_h(t)}</h5>")
                continue

            close_lists(target_level - 1)
            open_lists(target_level)

            if content.lower().startswith("interpretation:"):
                rest = content.split(":", 1)[1].strip() if ":" in content else ""
                out.append(
                    f"<li><strong>Interpretation:</strong> "
                    f"{_render_text_with_cites(rest, cite_patterns or [], num_to_title)}</li>"
                )
            else:
                out.append(f"<li>{_render_text_with_cites(content, cite_patterns or [], num_to_title)}</li>")
            continue

        close_lists(0)
        out.append(f"<p>{_render_text_with_cites(line.strip(), cite_patterns or [], num_to_title)}</p>")

    close_lists(0)
    return "".join(out)


# ---------------------------
# HTML builder
# ---------------------------

def _build_html(raw: Dict[str, Any]) -> str:
    inputs = raw.get("inputs", {}) or {}
    gene = inputs.get("gene_name", "") or ""
    m_cds = inputs.get("mutation_cds", "") or ""
    m_aa = inputs.get("mutation_aa", "") or ""

    mpos = inputs.get("aa_mut_start", None)
    try:
        mutation_pos = int(mpos)
    except Exception:
        mutation_pos = _extract_int(m_aa) or 0

    # Structure inference (kept)
    pdb_url = _infer_pdb_url(raw, gene)
    uniprot_id = _infer_uniprot_accession(raw)

    # NEW: Mol* viewer URL (AlphaFold-style)
    if uniprot_id:
        structure_viewer_url = f"https://molstar.org/viewer/?afdb={quote(uniprot_id)}"
        structure_link_url = f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
    elif pdb_url:
        structure_viewer_url = f"https://molstar.org/viewer/?url={quote(pdb_url, safe='')}&format=pdb"
        structure_link_url = pdb_url
    else:
        structure_viewer_url = ""
        structure_link_url = ""

    variant_summary = (raw.get("collector_output", {}) or {}).get("variant_summary", "") or ""

    impacted_sites = _collect_impacted_sites(raw)
    hyps = _as_list(raw.get("impact_hypotheses_output", {}).get("hypotheses"))
    ge_evals = _as_list(raw.get("gene_expression_evaluations"))
    references = _merge_references(raw)

    toc_links = [("#overview", "Overview")]

    ref_map = _create_reference_map(references)
    cite_patterns, num_to_title = _build_citation_index(references, ref_map) if references else ([], {})
    # Ensure cite_patterns is always a list (never None)
    if cite_patterns is None:
        cite_patterns = []
    if num_to_title is None:
        num_to_title = {}

    if hyps:
        toc_links.append(("#hypotheses", "Hypotheses"))
    if ge_evals:
        toc_links.append(("#gene-expression", "Differential gene expression"))

    # Collect per-hypothesis confidence levels
    hyp_confidences: List[Dict[str, str]] = []
    for hyp in (hyps or []):
        if not isinstance(hyp, dict):
            continue
        cl = (hyp.get("confidence_level") or "").strip()
        if cl:
            t = (hyp.get("title") or "Hypothesis").strip()
            hyp_confidences.append({"title": t, "level": cl})

    if hyp_confidences:
        toc_links.append(("#confidence", "Confidence"))
    if references:
        toc_links.append(("#references", "References"))

    toc_html = "".join(
        f'<a href="{_h(href)}" target="_self" aria-label="Jump to {_h(label)}">{_h(label)}</a>'
        for href, label in toc_links
    )

    input_rows = [
        ["Gene", _h(gene)],
        ["Mutation (CDS)", _h(m_cds)],
        ["Mutation (AA)", _h(m_aa)],
    ]
    inputs_block = _render_table(["Field", "Value"], input_rows)

    # NEW: structure block (static-ish figure, no custom animation)
    if structure_viewer_url:
        fallback_link = f"<p class='muted' style='margin-top:8px'>If the structure doesn't load, <a href='{_h(structure_link_url)}' target='_blank' rel='noopener'>view it here</a>.</p>" if structure_link_url else ""
        structure_block = (
            "<figure class='card'>"
            "<figcaption><strong>Predicted structure (AlphaFold DB / Mol*)</strong></figcaption>"
            f"<iframe class='molstar-frame' src='{_h(structure_viewer_url)}' "
            "allow='webgl; fullscreen' allowfullscreen "
            "loading='lazy' referrerpolicy='no-referrer' "
            "style='border: none;'></iframe>"
            f"{fallback_link}"
            "</figure>"
        )
    else:
        structure_block = (
            "<figure class='card'>"
            "<figcaption><strong>Predicted structure</strong></figcaption>"
            "<p class='muted'>No structure could be resolved for this report.</p>"
            "</figure>"
        )

    overview_parts: List[str] = []
    overview_parts.append(f"<div class='card'>{inputs_block}</div>")
    overview_parts.append(structure_block)  # NEW
    if variant_summary:
        overview_parts.append(
            "<div class='card'>"
            f"<p>{_render_text_with_cites(variant_summary, cite_patterns, num_to_title)}</p>"
            "</div>"
        )
    overview_html = "".join(overview_parts)

    # Kinase specificity
    kinase_html = ""
    if impacted_sites:
        cards: List[str] = []
        for site in impacted_sites:
            if not isinstance(site, dict):
                continue
            center_pos = int(site.get("center_pos") or 0)
            label = _clean_site_label(site.get("site_label") or "", "", center_pos)
            ev = site.get("evidence_info", "") or ""
            wt_15 = (site.get("wildtype_15mer") or "").strip()
            mut_15 = (site.get("mutant_15mer") or "").strip()
            events = _as_list(site.get("kinase_events"))

            align = ""
            if wt_15 and mut_15 and center_pos > 0 and mutation_pos > 0:
                align = _render_window_codeblock(label, center_pos, mutation_pos, wt_15, mut_15)

            position_html = f"<p class='muted'>Position: <strong>{_h(str(center_pos))}</strong></p>" if center_pos else ""

            cards.append(
                "<article class='card'>"
                f"<h3><strong>{_h(label)}</strong></h3>"
                + (f"<p class='muted'>{_h(ev)}</p>" if ev else "")
                + position_html
                + (align if align else "")
                + _render_kinase_events_table(events)
                + "</article>"
            )
        kinase_html = "".join(cards)

    # Hypotheses
    hypotheses_html = ""
    if hyps:
        parts: List[str] = []
        for hyp in hyps:
            if not isinstance(hyp, dict):
                continue

            title = hyp.get("title", "") or ""
            mech = hyp.get("mechanism_summary", "") or ""

            disease_ctx = _as_list(hyp.get("disease_context"))
            disease_html = ""
            if disease_ctx:
                disease_html = (
                    "<h4>Disease context</h4>"
                    "<ul>"
                    + "".join(
                        f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in disease_ctx
                    )
                    + "</ul>"
                )

            hyp_sites = _as_list(hyp.get("implicated_impacted_sites"))
            hyp_sites_html = ""
            if hyp_sites:
                site_bits: List[str] = []
                for s in hyp_sites:
                    if not isinstance(s, dict):
                        continue
                    raw_label = (s.get("site_label") or "").strip()
                    sev = s.get("evidence_info", "") or ""
                    wt_15 = (s.get("wildtype_15mer") or "").strip()
                    mut_15 = (s.get("mutant_15mer") or "").strip()
                    center_pos = _extract_int(raw_label) or _extract_int(sev) or 0
                    slabel = _clean_site_label(raw_label, "", center_pos)
                    events = _as_list(s.get("kinase_events"))

                    align = ""
                    if wt_15 and mut_15 and center_pos > 0 and mutation_pos > 0:
                        align = _render_window_codeblock(slabel, center_pos, mutation_pos, wt_15, mut_15)

                    position_html = f"<p class='muted'>Position: <strong>{_h(str(center_pos))}</strong></p>" if center_pos else ""

                    site_bits.append(
                        "<div class='card' style='margin:10px 0'>"
                        f"<h4><strong>{_h(slabel)}</strong></h4>"
                        + (f"<p class='muted'>{_render_text_with_cites(sev, cite_patterns, num_to_title)}</p>" if sev else "")
                        + position_html
                        + (align if align else "")
                        + _render_kinase_events_table(events)
                        + "</div>"
                    )
                hyp_sites_html = "<h4>Impacted sites</h4>" + "".join(site_bits)

            preds = _as_list(hyp.get("predicted_observations"))
            preds_html = ""
            if preds:
                preds_html = (
                    "<h4>Predicted observations</h4>"
                    "<ul>"
                    + "".join(
                        f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in preds
                    )
                    + "</ul>"
                )

            lab = hyp.get("lab_validation", {}) or {}
            lab_models = _as_list(lab.get("suggested_models"))
            lab_assays = _as_list(lab.get("assays"))
            lab_ctrls = _as_list(lab.get("controls_and_confounders"))
            lab_notes = lab.get("feasibility_notes", "") or ""

            lab_bits: List[str] = []
            if lab_models:
                lab_bits.append(
                    "<h4>Suggested models</h4>"
                    "<ul>" + "".join(f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in lab_models) + "</ul>"
                )
            if lab_assays:
                lab_bits.append(
                    "<h4>Assays</h4>"
                    "<ul>" + "".join(f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in lab_assays) + "</ul>"
                )
            if lab_ctrls:
                lab_bits.append(
                    "<h4>Controls and confounders</h4>"
                    "<ul>" + "".join(f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in lab_ctrls) + "</ul>"
                )
            if lab_notes:
                lab_bits.append(f"<h4>Feasibility notes</h4><p class='muted'>{_render_text_with_cites(lab_notes, cite_patterns, num_to_title)}</p>")
            lab_html = (
                "<details class='details' open>"
                "<summary>Laboratory validation</summary>"
                + "".join(lab_bits)
                + "</details>"
            ) if lab_bits else ""

            gev = hyp.get("gene_expression_validation", {}) or {}
            gev_html = ""
            if isinstance(gev, dict) and gev:
                is_app = gev.get("is_applicable", None)
                rationale = gev.get("rationale", "") or ""
                comparisons = _as_list(gev.get("proposed_comparisons"))
                markers = _as_list(gev.get("gene_sets_or_markers"))
                caveats = _as_list(gev.get("caveats"))

                bits: List[str] = []
                yn = _yesno(is_app) if is_app is not None else ""
                if yn:
                    bits.append(f"<p><strong>Applicable:</strong> {_h(yn)}</p>")
                if rationale:
                    bits.append(f"<p class='muted'>{_render_text_with_cites(rationale, cite_patterns, num_to_title)}</p>")
                if comparisons:
                    bits.append(
                        "<h4>Proposed comparisons</h4>"
                        "<ul>" + "".join(f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in comparisons) + "</ul>"
                    )
                if markers:
                    bits.append(
                        "<h4>Gene sets and markers</h4>"
                        "<ul>" + "".join(f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in markers) + "</ul>"
                    )
                if caveats:
                    bits.append(
                        "<h4>Caveats</h4>"
                        "<ul>" + "".join(f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in caveats) + "</ul>"
                    )

                if bits:
                    gev_html = (
                        "<details class='details' open>"
                        "<summary>Gene-expression validation</summary>"
                        + "".join(bits)
                        + "</details>"
                    )

            limitations = _as_list(hyp.get("limitations"))
            lim_html = ""
            if limitations:
                lim_html = (
                    "<details class='details' open>"
                    "<summary>Limitations</summary>"
                    "<ul>" + "".join(f"<li>{_render_text_with_cites(x, cite_patterns, num_to_title)}</li>" for x in limitations) + "</ul>"
                    "</details>"
                )

            clean_title = _strip_citation_patterns(title, cite_patterns) if title else ""
            parts.append(
                "<article class='card'>"
                + (f"<h3>{_h(clean_title)}</h3>" if clean_title else "<h3>Hypothesis</h3>")
                + (f"<p>{_render_text_with_cites(mech, cite_patterns, num_to_title)}</p>" if mech else "")
                + disease_html
                + (hyp_sites_html if hyp_sites_html else "")
                + (preds_html if preds_html else "")
                + (lab_html if lab_html else "")
                + (gev_html if gev_html else "")
                + (lim_html if lim_html else "")
                + "</article>"
            )

        hypotheses_html = "".join(parts)

    # Differential gene expression
    dge_html = ""
    if ge_evals:
        cards: List[str] = []
        for ge in ge_evals:
            if not isinstance(ge, dict):
                continue
            out = ge.get("gene_expression_output", {}) or {}

            genes = _as_list(out.get("genes_analyzed"))
            results = _as_list(out.get("results"))
            interpretation = out.get("interpretation", "") or ""
            group_filters = out.get("group_filters", {}) or {}

            genes_html = ""
            if genes:
                genes_html = "<div class='chips'>" + "".join(f"<span class='chip'>{_h(g)}</span>" for g in genes) + "</div>"

            gf_html = ""
            if group_filters:
                gf_items = []
                for k in ["group1_filter", "group2_filter"]:
                    if k in group_filters and str(group_filters.get(k, "")).strip():
                        label = "Group 1" if k == "group1_filter" else "Group 2"
                        filter_text = str(group_filters.get(k, "")).strip()
                        filter_parts = []
                        for line in filter_text.replace(";", "\n").split("\n"):
                            line = line.strip()
                            if not line:
                                continue
                            line = re.sub(r"^\s*==\s*", "", line)
                            if "==" in line:
                                parts2 = line.split("==", 1)
                                key = parts2[0].strip()
                                value = parts2[1].strip() if len(parts2) > 1 else ""
                                if key:
                                    filter_parts.append((key, value))
                            else:
                                filter_parts.append(("", line))

                        if filter_parts:
                            gf_items.append(f"<div style='margin-bottom:12px'><strong>{_h(label)}</strong>")
                            if filter_parts[0][0]:
                                gf_items.append("<ul style='margin-top:4px;margin-left:20px;list-style:none'>")
                                for key, value in filter_parts:
                                    if key and value:
                                        gf_items.append(f"<li style='font-family:monospace'><strong>{_h(key)}:</strong> {_h(value)}</li>")
                                    elif key:
                                        gf_items.append(f"<li style='font-family:monospace'><strong>{_h(key)}</strong></li>")
                                    elif value:
                                        gf_items.append(f"<li style='font-family:monospace'>{_h(value)}</li>")
                                gf_items.append("</ul>")
                            else:
                                gf_items.append("<ul style='margin-top:4px;margin-left:20px'>")
                                for _, value in filter_parts:
                                    gf_items.append(f"<li style='font-family:monospace'>{_h(value)}</li>")
                                gf_items.append("</ul>")
                            gf_items.append("</div>")

                if gf_items:
                    gf_html = "<div>" + "".join(gf_items) + "</div>"

            res_table = ""
            if results and all(isinstance(r, dict) for r in results):
                cols = []
                for r in results:
                    for k in r.keys():
                        if k not in cols:
                            cols.append(k)

                sig_results: List[Dict[str, Any]] = []
                nonsig_results: List[Dict[str, Any]] = []
                for r in results:
                    pval = _extract_pvalue(r)
                    if pval is not None and pval < 0.05:
                        sig_results.append(r)
                    else:
                        nonsig_results.append(r)

                sig_results.sort(key=lambda r: _extract_pvalue(r) if _extract_pvalue(r) is not None else 0.0)

                sig_rows = [[_format_ge_cell(c, r.get(c, "")) for c in cols] for r in sig_results]
                nonsig_rows = [[_format_ge_cell(c, r.get(c, "")) for c in cols] for r in nonsig_results]

                sig_table = _render_table(cols, sig_rows) if sig_rows else ""
                nonsig_table = _render_table(cols, nonsig_rows) if nonsig_rows else ""

                if sig_table and nonsig_table:
                    res_table = sig_table + "<div class='barrier'>Non-significant (p ≥ 0.05)</div>" + nonsig_table
                else:
                    res_table = sig_table or nonsig_table

            interp_html = _render_markdownish(interpretation, cite_patterns=cite_patterns, num_to_title=num_to_title) if interpretation else ""

            cards.append(
                "<article class='card'>"
                "<h3>Differential expression</h3>"
                + (genes_html if genes_html else "")
                + (
                    (
                        "<details class='details' open>"
                        "<summary>Group definitions</summary>"
                        + gf_html
                        + "</details>"
                    ) if gf_html else ""
                )
                + (res_table if res_table else "<p class='muted'>No differential-expression table was returned for this run.</p>")
                + (
                    (
                        "<details class='details' open>"
                        "<summary>Interpretation</summary>"
                        + interp_html
                        + "</details>"
                    ) if interp_html else ""
                )
                + "</article>"
            )

        dge_html = "".join(cards)

    # References (numbered)
    refs_html = ""
    if references:
        items: List[str] = []
        for idx, s in enumerate(references, start=1):
            if not isinstance(s, dict):
                continue
            t = s.get("title", "") or ""
            u = s.get("url", "") or ""
            summ = s.get("summary", "") or ""
            rid = s.get("ref_id", "") or ""

            if u:
                items.append(
                    f"<article class='card' id='ref-{idx}'>"
                    f"<p><strong>[{idx}]</strong> <a href='{_h(u)}' target='_blank' rel='noopener'>{_h(t or u)}</a></p>"
                    + (f"<p class='muted'>{_render_text_with_cites(summ, cite_patterns, num_to_title)}</p>" if summ else "")
                    + (f"<p class='muted' style='margin-top:6px'>ID: <code>{_h(rid)}</code></p>" if rid else "")
                    + "</article>"
                )
            else:
                items.append(
                    f"<article class='card' id='ref-{idx}'>"
                    f"<p><strong>[{idx}]</strong> {_h(t)}</p>"
                    + (f"<p class='muted'>{_render_text_with_cites(summ, cite_patterns, num_to_title)}</p>" if summ else "")
                    + (f"<p class='muted' style='margin-top:6px'>ID: <code>{_h(rid)}</code></p>" if rid else "")
                    + "</article>"
                )
        refs_html = "".join(items)

    title = f"{gene} {m_aa}".strip() or "Variant Report"

    html_doc = f"""<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{_h(title)}</title>
  <style>
    :root{{
      --bg:#f8fafc; --text:#0f172a; --muted:#475569; --card:#ffffff;
      --accent:#2563eb; --accent-2:#0ea5e9; --line:#e2e8f0; --chip-bg:#f1f5f9;
      --link:#2563eb; --link-visited:#6366f1;
    }}
    *{{box-sizing:border-box}}
    body{{margin:0;background:var(--bg);color:var(--text);font-family:ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,"Apple Color Emoji","Segoe UI Emoji";line-height:1.55}}
    a{{color:var(--link);text-decoration:none}}
    a:visited{{color:var(--link-visited)}}
    a:hover{{text-decoration:underline}}
    .wrap{{max-width:1100px;margin:0 auto;padding:24px}}
    header{{background:linear-gradient(135deg,var(--accent) 0%,var(--accent-2) 100%);color:white;padding:22px 0;position:relative;overflow:hidden}}
    header .wrap{{padding:0 24px;position:relative;z-index:2}}
    h1{{margin:0;font-size:clamp(1.35rem, 2.2vw + 1rem, 2.1rem);letter-spacing:-.02em}}
    nav{{position:sticky;top:0;background:rgba(248,250,252,.92);backdrop-filter:blur(10px);border-bottom:1px solid var(--line);z-index:10}}
    nav .wrap{{display:flex;gap:10px;flex-wrap:wrap;align-items:center;padding:12px 24px}}
    nav a{{display:inline-flex;align-items:center;padding:8px 12px;border-radius:999px;background:var(--chip-bg);border:1px solid var(--line);font-weight:700;font-size:.92rem}}
    main{{padding:22px 0}}
    section{{margin:18px 0 28px 0}}
    .section-h{{display:flex;align-items:flex-end;justify-content:space-between;gap:12px;margin-bottom:8px}}
    .section-h h2{{margin:0;font-size:clamp(1.5rem, 1.8vw + 1rem, 2rem);font-weight:800;letter-spacing:-0.02em}}
    .section-desc{{margin-top:0;margin-bottom:16px;font-size:.95rem}}
    h3{{font-size:1.15rem;font-weight:700}}
    h4{{font-size:1rem;font-weight:600}}
    h5{{font-size:0.95rem;font-weight:600}}
    .card{{background:var(--card);border:1px solid var(--line);border-radius:18px;padding:14px 14px;box-shadow:0 1px 0 rgba(15,23,42,.03), 0 10px 25px rgba(15,23,42,.05);margin:12px 0}}
    .muted{{color:var(--muted)}}
    .chips{{display:flex;gap:8px;flex-wrap:wrap;margin:8px 0 12px 0}}
    .chip{{display:inline-flex;align-items:center;padding:6px 10px;border-radius:999px;background:var(--chip-bg);border:1px solid var(--line);font-weight:800;font-size:.82rem}}
    figure.card{{margin:12px 0}}
    figcaption{{margin-bottom:8px}}
    .tablewrap{{width:100%;overflow:auto;border:1px solid var(--line);border-radius:14px}}
    table.table{{border-collapse:separate;border-spacing:0;width:100%;min-width:640px;background:var(--card)}}
    .table th, .table td{{padding:10px 12px;border-bottom:1px solid var(--line);vertical-align:top}}
    .table thead th{{position:sticky;top:0;background:var(--bg);z-index:1;text-align:left}}
    .table tbody tr:nth-child(odd) td{{background:rgba(241,245,249,.65)}}
    .details summary{{cursor:pointer;font-weight:900}}
    .ref-link{{color:var(--accent);text-decoration:none;font-weight:700}}
    .ref-link:hover{{text-decoration:underline}}
    .details h4{{margin:.9rem 0 .4rem 0}}
    .details h5{{margin:.7rem 0 .35rem 0; font-size: 1.0rem;}}
    .details p{{margin:.5rem 0}}
    .details ul {{
      margin: .4rem 0 .8rem 1.2rem;
      padding-left: 1.1rem;
      list-style: disc;
    }}
    .details ul ul {{
      margin-top: .3rem;
      list-style: circle;
    }}
    .barrier{{margin:10px 0;padding:8px 12px;border-radius:10px;background:#0f172a;color:#e2e8f0;font-weight:700;text-transform:uppercase;letter-spacing:.04em;font-size:.8rem}}
    .codeblock{{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono","Courier New",monospace;white-space:pre;tab-size:5;line-height:1.35}}
    .codeblock .mut{{background:rgba(37,99,235,.15);border-bottom:2px solid var(--accent);font-weight:800}}

    /* NEW: Mol* embed */
    .molstar-frame {{
      width: 100%;
      height: 520px;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: var(--card);
    }}

    @media print {{
      nav{{display:none}}
      header{{print-color-adjust:exact}}
      a[href]::after{{content:" (" attr(href) ")";font-size:.85em;color:var(--muted)}}
      .card{{break-inside:avoid}}
      section{{break-inside:avoid-page}}
      .molstar-frame{{display:none}}
    }}
  </style>
</head>
<body>
  <header>
    <div class="wrap">
      <h1>{_h(title)}</h1>
    </div>
  </header>

  <nav aria-label="Table of Contents">
    <div class="wrap">
      {toc_html}
    </div>
  </nav>

  <main class="wrap">

    <section id="overview" aria-label="Overview">
      <div class="section-h">
        <h2>Overview</h2>
      </div>
      <p class="section-desc muted">What the variant is and why it's being analyzed</p>
      {overview_html}
    </section>

    {""
      if not hyps else
    f'''
    <section id="hypotheses" aria-label="Hypotheses">
      <div class="section-h">
        <h2>Hypotheses</h2>
      </div>
      <p class="section-desc muted">Mechanistic explanations generated from the overview and literature.</p>
      {hypotheses_html}
    </section>
    '''
    }

    {""
      if not ge_evals else
    f'''
    <section id="gene-expression" aria-label="Differential gene expression">
      <div class="section-h">
        <h2>Differential gene expression</h2>
      </div>
      <p class="section-desc muted">Transcriptomic evidence used to support, deprioritize, or refine the hypotheses above (when applicable).</p>
      {dge_html}
    </section>
    '''
    }

    {""
      if not hyp_confidences else
    f'''
    <section id="confidence" aria-label="Confidence">
      <div class="section-h">
        <h2>Confidence</h2>
      </div>
      <p class="section-desc muted">Confidence level assigned to each hypothesis in this report.</p>
      {"".join(
          "<article class='card'>"
          f"<p><strong>{_h(c['title'])}</strong></p>"
          f"<p>Confidence: <strong>{_h(c['level'])}</strong></p>"
          "</article>"
          for c in hyp_confidences
      )}
    </section>
    '''
    }

    {""
      if not references else
    f'''
    <section id="references" aria-label="References">
      <div class="section-h">
        <h2>References</h2>
      </div>
      <p class="section-desc muted">Full bibliography for all inline citations in this report.</p>
      {refs_html}
    </section>
    '''
    }

  </main>
</body>
</html>"""
    return html_doc


# ---------------------------
# MCP tool
# ---------------------------

@mcp.tool()
def generate_and_store_html_report_from_raw_report(
    file_path: str,
    dir_path: str,
    filename: str,
) -> Dict[str, Any]:
    """
    Generate an HTML report from a JSON raw_report file and store it to a file.

    If the raw report does not contain any hypotheses, HTML generation is skipped
    and the function returns a success status with an empty file path.

    Args:
        file_path: Path to the raw report JSON file
        dir_path: Directory where the HTML report should be saved
        filename: Filename for the HTML report (without extension)

    Output rules:
      - Top block shows only: Gene, Mutation (CDS), Mutation (AA)
      - No internal/system terminology or status fields are displayed
      - No notes are displayed
      - Phosphosite windows are not displayed
      - Hypothesis IDs are not displayed
      - Hypothesis confidence levels are shown in a dedicated section before References
      - Boolean values are rendered as Yes/No where applicable
      - Gene-expression interpretation is rendered with headings and nested bullets
      - Hypothesis gene-expression analysis plan is omitted
    """
    if not file_path:
        return {"status": "error", "message": "file_path is empty"}

    try:
        with open(file_path, "r", encoding="utf-8") as f:
            raw_report_content = f.read()
    except FileNotFoundError:
        return {"status": "error", "message": f"File not found: {file_path}"}
    except Exception as e:
        return {"status": "error", "message": f"Failed to read file {file_path}: {str(e)}"}

    if not raw_report_content or len(raw_report_content) == 0:
        return {"status": "error", "message": "raw_report file is empty"}

    try:
        data = json.loads(raw_report_content)
    except Exception as e:
        return {"status": "error", "message": f"raw_report is not valid JSON: {str(e)}"}

    hyps = _as_list((data.get("impact_hypotheses_output") or {}).get("hypotheses"))
    if not hyps:
        return {"status": "success", "file_path": ""}

    raw_reports_path = Path(RAW_REPORTS_PATH)
    file_path_obj = Path(file_path)

    try:
        relative_path = file_path_obj.relative_to(raw_reports_path)
        subfolder = relative_path.parent if relative_path.parent != Path('.') else None
    except ValueError:
        subfolder = None

    output_dir = (Path(dir_path) / subfolder) if subfolder else Path(dir_path)

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        return {"status": "error", "message": f"Failed to create directory {output_dir}: {str(e)}"}

    try:
        print(f"[DEBUG] Building HTML for: {filename}")
        html_output = _build_html(data)
        print(f"[DEBUG] HTML build completed for: {filename}")
    except Exception as e:
        return {"status": "error", "message": f"Failed to build HTML: {str(e)}"}

    output_file = output_dir / filename

    try:
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(html_output)
    except Exception as e:
        return {"status": "error", "message": f"Failed to write files: {str(e)}"}

    return {"status": "success", "file_path": str(output_file)}


def rerun(delete: bool = False) -> None:
    """
    Generate HTML reports for all raw report files.

    Args:
        delete: If True, delete all existing HTML reports before generating new ones.
    """
    raw_reports_dir = Path(RAW_REPORTS_PATH)
    html_reports_dir = Path(HTML_REPORTS_PATH)

    if not raw_reports_dir.exists():
        print(f"Error: Raw reports directory does not exist: {RAW_REPORTS_PATH}")
        return

    html_reports_dir.mkdir(parents=True, exist_ok=True)

    if delete:
        for html_file in html_reports_dir.rglob("*.html"):
            try:
                html_file.unlink()
                print(f"Deleted: {html_file}")
            except Exception as e:
                print(f"Warning: Failed to delete {html_file}: {str(e)}")

    raw_json_files = list(raw_reports_dir.rglob("*.json"))

    if not raw_json_files:
        print(f"No JSON files found in {RAW_REPORTS_PATH}")
        return

    print(f"Processing {len(raw_json_files)} raw report file(s)...")

    success_count = 0
    error_count = 0

    for raw_file in raw_json_files:
        filename = raw_file.stem
        print(f"\n[DEBUG] Starting to process raw report: {raw_file.name}")

        result = generate_and_store_html_report_from_raw_report(
            file_path=str(raw_file),
            dir_path=str(html_reports_dir),
            filename=f"{filename}.html"
        )

        if result.get("status") == "success":
            print(f"Generated: {result.get('file_path')}")
            success_count += 1
        else:
            print(f"Error processing {raw_file.name}: {result.get('message')}")
            error_count += 1

    print(f"\nCompleted: {success_count} successful, {error_count} errors")


if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) > 1:
        if sys.argv[1] == "rerun":
            delete_flag = "--delete" in sys.argv
            rerun(delete=delete_flag)
        elif sys.argv[1] == "figure":
            script_dir = Path(__file__).resolve().parent
            dummy_json = script_dir / "data" / "dummy.json"
            output_html = script_dir / "data" / "dummy_report.html"
            result = generate_and_store_html_report_from_raw_report(
                file_path=str(dummy_json),
                dir_path=str(script_dir / "data"),
                filename="dummy_report.html",
            )
            if result.get("status") == "success":
                print(f"Figure HTML written to: {result.get('file_path')}")
            else:
                print(f"Error: {result.get('message')}")
                sys.exit(1)
        else:
            assert False, "Invalid argument"
    else:
        mcp.run(transport="stdio")
