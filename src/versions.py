"""Upstream source version fetcher for clinvar-ingest.

clinvar-ingest pulls from three logical sources:
  - infores:clinvar       — the VCF, .tbi, and submission_summary
  - infores:medgen        — MedGenIDMappings
  - infores:mondo         — mondo.sssom.tsv (used here for label lookup,
                             not as a true mapping rewire — track at primary)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from kozahub_metadata_schema import (
    now_iso,
    urls_from_download_yaml,
    version_from_http_last_modified,
)


INGEST_DIR = Path(__file__).resolve().parents[1]
DOWNLOAD_YAML = INGEST_DIR / "download.yaml"


def get_source_versions() -> list[dict[str, Any]]:
    clinvar_urls = urls_from_download_yaml(DOWNLOAD_YAML, contains=["pub/clinvar"])
    medgen_urls = urls_from_download_yaml(DOWNLOAD_YAML, contains=["pub/medgen"])
    mondo_urls = urls_from_download_yaml(DOWNLOAD_YAML, contains=["data.monarchinitiative.org/mappings"])
    now = now_iso()

    sources: list[dict[str, Any]] = []

    if clinvar_urls:
        ver, method = version_from_http_last_modified(clinvar_urls[0])
        sources.append({
            "id": "infores:clinvar",
            "name": "ClinVar",
            "urls": clinvar_urls,
            "version": ver,
            "version_method": method,
            "retrieved_at": now,
        })

    if medgen_urls:
        ver, method = version_from_http_last_modified(medgen_urls[0])
        sources.append({
            "id": "infores:medgen",
            "name": "MedGen",
            "urls": medgen_urls,
            "version": ver,
            "version_method": method,
            "retrieved_at": now,
        })

    if mondo_urls:
        ver, method = version_from_http_last_modified(mondo_urls[0])
        sources.append({
            "id": "infores:mondo",
            "name": "Mondo Disease Ontology (SSSOM)",
            "urls": mondo_urls,
            "version": ver,
            "version_method": method,
            "retrieved_at": now,
        })

    return sources
