from __future__ import annotations

from dataclasses import dataclass, field
import re
import time
from typing import List, Optional, Sequence, Tuple
from urllib.parse import parse_qs, urljoin, urlparse

import requests
from bs4 import BeautifulSoup

__all__ = [
    "GorillaError",
    "GeneHit",
    "GORecord",
    "RunMetadata",
    "GorillaResult",
    "submit_ranked_genes",
    "parse_results_html",
    "GorillaClient",
]

DEFAULT_BASE_URL = "https://cbl-gorilla.cs.technion.ac.il"
RUN_PATH_SEGMENT = "GOrilla"


class GorillaError(RuntimeError):
    """Raised when a submission to GOrilla fails."""


@dataclass(frozen=True)
class GeneHit:
    """A gene reported for a GO term."""

    symbol: str
    description: str


@dataclass(frozen=True)
class GORecord:
    """A single GO enrichment entry returned by GOrilla."""

    go_id: str
    description: str
    p_value: Optional[float]
    fdr_q_value: Optional[float]
    enrichment: Optional[float]
    counts: Tuple[int, int, int, int]
    genes: List[GeneHit]


@dataclass
class RunMetadata:
    """Summary metadata reported by GOrilla."""

    species: Optional[str] = None
    submitted_gene_count: Optional[int] = None
    recognized_gene_count: Optional[int] = None
    recognized_by_symbol: Optional[int] = None
    recognized_by_other_ids: Optional[int] = None
    duplicates_removed: Optional[int] = None
    total_after_duplicates: Optional[int] = None
    associated_with_go_terms: Optional[int] = None
    warnings: List[str] = field(default_factory=list)


@dataclass
class GorillaResult:
    """Structured representation of a GOrilla enrichment run."""

    run_id: str
    result_url: str
    go_terms: List[GORecord]
    metadata: RunMetadata
    raw_html: str


class GorillaClient:
    """Thin convenience wrapper that reuses an HTTP session for multiple runs."""

    def __init__(
        self,
        *,
        base_url: str = DEFAULT_BASE_URL,
        session: Optional[requests.Session] = None,
        timeout: float = 120.0,
        processing_timeout: Optional[float] = None,
        poll_interval: float = 2.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.session = session or requests.Session()
        self.timeout = timeout
        self.processing_timeout = processing_timeout
        self.poll_interval = poll_interval

    def submit_ranked_genes(
        self,
        genes: Sequence[str],
        *,
        species: str = "HOMO_SAPIENS",
        ontology: str = "proc",
        pvalue_threshold: float = 1e-3,
        fast_mode: bool = True,
        analysis_name: Optional[str] = None,
        user_email: Optional[str] = None,
        output_excel: bool = False,
        output_unresolved: bool = False,
        output_revigo: bool = False,
        background_genes: Optional[Sequence[str]] = None,
        run_mode: Optional[str] = None,
        processing_timeout: Optional[float] = None,
        poll_interval: Optional[float] = None,
    ) -> GorillaResult:
        return submit_ranked_genes(
            genes,
            species=species,
            ontology=ontology,
            pvalue_threshold=pvalue_threshold,
            fast_mode=fast_mode,
            analysis_name=analysis_name,
            user_email=user_email,
            output_excel=output_excel,
            output_unresolved=output_unresolved,
            output_revigo=output_revigo,
            background_genes=background_genes,
            run_mode=run_mode,
            processing_timeout=(
                processing_timeout if processing_timeout is not None else self.processing_timeout
            ),
            poll_interval=poll_interval if poll_interval is not None else self.poll_interval,
            base_url=self.base_url,
            session=self.session,
            timeout=self.timeout,
        )


def submit_ranked_genes(
    genes: Sequence[str],
    *,
    species: str = "HOMO_SAPIENS",
    ontology: str = "proc",
    pvalue_threshold: float = 1e-3,
    fast_mode: bool = True,
    analysis_name: Optional[str] = None,
    user_email: Optional[str] = None,
    output_excel: bool = False,
    output_unresolved: bool = False,
    output_revigo: bool = False,
    background_genes: Optional[Sequence[str]] = None,
    run_mode: Optional[str] = None,
    processing_timeout: Optional[float] = None,
    poll_interval: float = 2.0,
    base_url: str = DEFAULT_BASE_URL,
    session: Optional[requests.Session] = None,
    timeout: float = 120.0,
) -> GorillaResult:
    """Submit a ranked gene list to GOrilla and return the parsed results.

    Parameters
    ----------
    genes:
        Ordered gene identifiers (top-ranked genes first).
    species:
        A supported GOrilla species identifier (e.g., ``HOMO_SAPIENS``).
    ontology:
        Which ontology to query: ``proc`` (biological process), ``func``,
        ``comp``, or ``all``.
    pvalue_threshold:
        P-value cutoff reported by GOrilla.
    fast_mode:
        Whether to speed up the mHG algorithm (matches the web form default).
    analysis_name, user_email:
        Optional metadata fields accepted by GOrilla.
    output_excel, output_unresolved, output_revigo:
        Optional toggles corresponding to the advanced checkboxes on the site.
    background_genes:
        Optional background gene list. When provided (or when ``run_mode='hg'``),
        GOrilla runs the target-versus-background hypergeometric test.
    run_mode:
        Explicitly select ``'mhg'`` (ranked list) or ``'hg'`` (target/background).
        Defaults to ``'hg'`` when a background list is supplied, otherwise ``'mhg'``.
    processing_timeout:
        Maximum wall-clock time (seconds) to wait for the web results to become
        available. Defaults to ``timeout`` if omitted.
    poll_interval:
        Delay between polling attempts while waiting for the results page.
    base_url:
        Base URL of the GOrilla deployment.
    session:
        Optional ``requests.Session`` to reuse connections.
    timeout:
        Request timeout (seconds).
    """

    if not genes:
        raise ValueError("At least one gene is required.")

    if poll_interval <= 0:
        raise ValueError("poll_interval must be greater than zero.")

    created_session = False
    if session is None:
        session = requests.Session()
        created_session = True

    target_set = "\n".join(_normalize_gene(g) for g in genes if g)
    if not target_set.strip():
        raise ValueError("Gene list cannot be empty after normalization.")

    if background_genes is not None:
        background_set = "\n".join(_normalize_gene(g) for g in background_genes if g)
        if not background_set.strip():
            raise ValueError("Background gene list cannot be empty after normalization.")
    else:
        background_set = None

    resolved_run_mode = _resolve_run_mode(run_mode, background_set is not None)
    if resolved_run_mode == "hg" and background_set is None:
        raise ValueError("background_genes must be supplied when run_mode='hg'.")
    if resolved_run_mode == "mhg" and background_set is not None:
        raise ValueError("background_genes provided but run_mode set to 'mhg'.")

    submit_url = urljoin(base_url.rstrip("/") + "/", "servlet/GOrilla")
    files = _build_multipart_payload(
        run_mode=resolved_run_mode,
        target_set=target_set,
        background_set=background_set,
        species=species,
        ontology=ontology,
        pvalue_threshold=pvalue_threshold,
        fast_mode=fast_mode,
        analysis_name=analysis_name,
        user_email=user_email,
        output_excel=output_excel,
        output_unresolved=output_unresolved,
        output_revigo=output_revigo,
    )

    total_timeout = processing_timeout if processing_timeout is not None else timeout
    if total_timeout <= 0:
        raise ValueError("processing_timeout/timeout must be greater than zero.")
    deadline = time.monotonic() + total_timeout
    try:
        response = session.post(submit_url, files=files, timeout=timeout)
        response.raise_for_status()
        result_url, html = _resolve_results_page(
            session=session,
            base_url=base_url,
            initial_url=response.url,
            initial_html=response.text,
            deadline=deadline,
            poll_interval=poll_interval,
        )
    except requests.RequestException as exc:
        raise GorillaError(f"Submission to GOrilla failed: {exc}") from exc
    finally:
        if created_session:
            session.close()

    return parse_results_html(html, result_url=result_url)


def parse_results_html(html: str, *, result_url: Optional[str] = None) -> GorillaResult:
    """Parse a GOrilla results HTML page into structured objects."""

    soup = BeautifulSoup(html, "html.parser")
    go_terms = _parse_go_table(soup)
    metadata = _parse_metadata(soup)
    run_id = _extract_run_id(result_url)
    return GorillaResult(
        run_id=run_id,
        result_url=result_url or "",
        go_terms=go_terms,
        metadata=metadata,
        raw_html=html,
    )


def _build_multipart_payload(
    *,
    run_mode: str,
    target_set: str,
    background_set: Optional[str],
    species: str,
    ontology: str,
    pvalue_threshold: float,
    fast_mode: bool,
    analysis_name: Optional[str],
    user_email: Optional[str],
    output_excel: bool,
    output_unresolved: bool,
    output_revigo: bool,
) -> List[Tuple[str, Tuple[Optional[str], str]]]:
    """Construct the multipart fields that mirror the web form submission."""

    files: List[Tuple[str, Tuple[Optional[str], str]]] = []

    def add_field(name: str, value: str) -> None:
        files.append((name, (None, value)))

    add_field("application", "gorilla")
    add_field("species", species)
    add_field("run_mode", run_mode)
    add_field("target_set", target_set)
    add_field("db", ontology)
    add_field("pvalue_thresh", f"{pvalue_threshold:.12g}")
    add_field("run_gogo_button", "Search Enriched GO terms")

    if run_mode == "hg" and background_set:
        add_field("background_set", background_set)

    if run_mode == "mhg" and fast_mode:
        add_field("fast_mode", "on")
    if analysis_name:
        add_field("analysis_name", analysis_name)
    if user_email:
        add_field("user_email", user_email)
    if output_excel:
        add_field("output_excel", "on")
    if output_unresolved:
        add_field("output_unresolved", "on")
    if output_revigo:
        add_field("output_revigo", "on")

    return files


def _parse_go_table(soup: BeautifulSoup) -> List[GORecord]:
    table = soup.find("table", id="table1")
    if table is None:
        text = soup.get_text(" ", strip=True)
        lowered = text.lower()
        if "no go enrichment found" in lowered or "no go terms" in lowered:
            return []
        snippet = text[:200]
        raise GorillaError(
            "Could not locate the GO results table in the HTML response; "
            f"got the following page snippet instead: {snippet!r}"
        )

    records: List[GORecord] = []
    for row in table.find_all("tr")[1:]:
        cells = row.find_all("td")
        if len(cells) < 6:
            continue

        go_identifier = cells[0].get_text(strip=True)
        description = cells[1].get_text(" ", strip=True)
        p_value = _parse_float(cells[2].get_text(strip=True))
        fdr_q_value = _parse_float(cells[3].get_text(strip=True))

        enrichment_text = cells[4].get_text(strip=True)
        enrichment, counts = _parse_enrichment(enrichment_text)

        genes_div = cells[5].find("div")
        genes = _parse_gene_block(genes_div.get_text("\n", strip=True)) if genes_div else []

        records.append(
            GORecord(
                go_id=go_identifier,
                description=description,
                p_value=p_value,
                fdr_q_value=fdr_q_value,
                enrichment=enrichment,
                counts=counts,
                genes=genes,
            )
        )

    return records


def _parse_metadata(soup: BeautifulSoup) -> RunMetadata:
    metadata = RunMetadata()
    metadata.warnings = [
        header.get_text(" ", strip=True)
        for header in soup.find_all("h1")
        if header.get_text(strip=True)
    ]

    for paragraph in soup.find_all("p"):
        text = paragraph.get_text(" ", strip=True)
        if not text:
            continue
        if text.startswith("Species used:"):
            metadata.species = text.split(":", 1)[1].strip()
            continue
        if text.startswith("The system has recognized"):
            _populate_counts_from_summary(metadata, text)

    return metadata


def _populate_counts_from_summary(metadata: RunMetadata, text: str) -> None:
    patterns = {
        "recognized_gene_count": r"recognized (\d+) genes out of (\d+)",
        "recognized_by_symbol": r"(\d+) genes were recognized by gene symbol and (\d+) genes by other gene IDs",
        "duplicates_removed": r"(\d+) duplicate genes were removed .* total of (\d+)",
        "associated_with_go_terms": r"Only (\d+) of these genes are associated with a GO term",
    }

    m = re.search(patterns["recognized_gene_count"], text)
    if m:
        metadata.recognized_gene_count = int(m.group(1))
        metadata.submitted_gene_count = int(m.group(2))

    m = re.search(patterns["recognized_by_symbol"], text)
    if m:
        metadata.recognized_by_symbol = int(m.group(1))
        metadata.recognized_by_other_ids = int(m.group(2))

    m = re.search(patterns["duplicates_removed"], text)
    if m:
        metadata.duplicates_removed = int(m.group(1))
        metadata.total_after_duplicates = int(m.group(2))

    m = re.search(patterns["associated_with_go_terms"], text)
    if m:
        metadata.associated_with_go_terms = int(m.group(1))


def _parse_gene_block(raw_text: str) -> List[GeneHit]:
    genes: List[GeneHit] = []
    for line in raw_text.splitlines():
        clean = line.strip()
        if not clean:
            continue
        symbol, description = _split_gene_line(clean)
        genes.append(GeneHit(symbol=symbol, description=description))
    return genes


def _split_gene_line(line: str) -> Tuple[str, str]:
    parts = re.split(r"\s*-\s*", line, maxsplit=1)
    if len(parts) == 2:
        return parts[0].strip(), parts[1].strip()
    return parts[0].strip(), ""


def _parse_enrichment(value: str) -> Tuple[Optional[float], Tuple[int, int, int, int]]:
    counts = (0, 0, 0, 0)
    value = value.strip()
    if not value:
        return None, counts

    match = re.match(r"([0-9.+-Ee]+)\s*\(([^)]+)\)", value)
    if match:
        enrichment = _parse_float(match.group(1))
        counts = _parse_counts(match.group(2))
        return enrichment, counts

    # No explicit enrichment factor, just counts.
    return None, _parse_counts(value)


def _parse_counts(group: str) -> Tuple[int, int, int, int]:
    parts = [part.strip() for part in group.replace("(", "").replace(")", "").split(",") if part.strip()]
    numbers = []
    for part in parts[:4]:
        try:
            numbers.append(int(part))
        except ValueError:
            numbers.append(0)
    while len(numbers) < 4:
        numbers.append(0)
    return tuple(numbers)  # type: ignore[return-value]


def _parse_float(value: str) -> Optional[float]:
    clean = value.strip()
    if not clean or clean.upper() in {"NA", "N/A"}:
        return None
    try:
        return float(clean)
    except ValueError:
        return None


def _normalize_gene(gene: str) -> str:
    return gene.strip()


def _extract_run_id(result_url: Optional[str]) -> str:
    if not result_url:
        return "unknown"

    parsed = urlparse(result_url)
    query_ids = parse_qs(parsed.query).get("id")
    if query_ids:
        return query_ids[0]

    segments = [seg for seg in parsed.path.split("/") if seg]
    try:
        idx = segments.index(RUN_PATH_SEGMENT)
        run_id = segments[idx + 1]
        return run_id
    except (ValueError, IndexError):
        return "unknown"


def _resolve_run_mode(run_mode: Optional[str], has_background: bool) -> str:
    if run_mode is None:
        return "hg" if has_background else "mhg"
    normalized = run_mode.lower()
    if normalized not in {"mhg", "hg"}:
        raise ValueError("run_mode must be either 'mhg' or 'hg'.")
    return normalized
    parsed = urlparse(result_url)
    query_ids = parse_qs(parsed.query).get("id")
    if query_ids:
        return query_ids[0]
    segments = [seg for seg in parsed.path.split("/") if seg]
    try:
        idx = segments.index(RUN_PATH_SEGMENT)
        run_id = segments[idx + 1]
        return run_id
    except (ValueError, IndexError):
        return "unknown"


def _resolve_results_page(
    *,
    session: requests.Session,
    base_url: str,
    initial_url: str,
    initial_html: str,
    deadline: float,
    poll_interval: float = 2.0,
) -> Tuple[str, str]:
    """Return a tuple of (final_result_url, html) waiting for processing to finish."""

    if _page_signals_completion(initial_html):
        return initial_url, initial_html

    run_id = _extract_run_id(initial_url)
    candidate_urls = []
    if run_id != "unknown":
        candidate_urls.append(_build_result_html_url(base_url, run_id))
    candidate_urls.append(initial_url)

    last_error: Optional[str] = None
    while True:
        remaining_time = deadline - time.monotonic()
        if remaining_time <= 0:
            break
        remaining = max(0.5, remaining_time)

        for url in candidate_urls:
            try:
                response = session.get(url, timeout=remaining)
            except requests.RequestException as exc:
                last_error = str(exc)
                continue

            if response.status_code == 200:
                html = response.text
                if _page_signals_completion(html):
                    return response.url, html
                last_error = "processing not finished"
            elif response.status_code == 404:
                last_error = "results not ready (404)"
            else:
                last_error = f"unexpected status {response.status_code}"

        sleep_time = deadline - time.monotonic()
        if sleep_time <= 0:
            break
        sleep_for = min(poll_interval, max(0.1, sleep_time))
        if sleep_for <= 0:
            break
        time.sleep(sleep_for)

    raise GorillaError(
        "Timed out waiting for GOrilla results to become available"
        + (f" ({last_error})" if last_error else "")
    )


def _build_result_html_url(base_url: str, run_id: str) -> str:
    return urljoin(base_url.rstrip("/") + "/", f"{RUN_PATH_SEGMENT}/{run_id}/GOResults.html")


def _page_signals_completion(html: str) -> bool:
    if 'id="table1"' in html:
        return True
    lowered = html.lower()
    completion_markers = [
        "no go enrichment found",
        "no go terms with an enrichment",
        "the system has recognized",
        "results page will be available",
    ]
    return any(marker in lowered for marker in completion_markers)
