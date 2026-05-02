"""Emit output/release-metadata.yaml for clinvar-ingest.

Standard boilerplate — content is in src/versions.py and the schema package.
"""

from __future__ import annotations

import sys
from pathlib import Path

INGEST_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(INGEST_DIR / "src"))

from versions import get_source_versions  # noqa: E402
from kozahub_metadata_schema.writer import write_metadata  # noqa: E402


if __name__ == "__main__":
    src = INGEST_DIR / "src"
    transform_paths = list(src.rglob("*.py")) + list(src.rglob("*.yaml"))

    output_dir = INGEST_DIR / "output"
    # Default to globbing every TSV / NT file in output/ as artifacts.
    # Override the artifacts list explicitly if your ingest produces a fixed set.
    artifacts = sorted(
        p.name
        for p in output_dir.glob("*")
        if p.is_file() and p.suffix in {".tsv", ".gz", ".jsonl", ".nt"}
    )

    metadata = write_metadata(
        ingest_name="clinvar-ingest",
        source_versions=get_source_versions(),
        transform_paths=transform_paths,
        artifacts=artifacts,
        output_dir=output_dir,
    )
    print(f"Wrote {output_dir / 'release-metadata.yaml'}")
    print(f"  build_version: {metadata['build_version']}")
    for s in metadata["sources"]:
        print(
            f"  source {s['id']}: version={s['version']} via {s['version_method']} "
            f"({len(s.get('urls') or [])} url(s))"
        )
