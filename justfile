# clinvar-ingest justfile

# Package directory
PKG := "src"

# Explicitly enumerate transforms
TRANSFORMS := "transform"

# List all commands
_default:
    @just --list

# Install dependencies
[group('project management')]
install:
    uv sync --group dev

# Full pipeline: test -> download -> preprocess -> transform -> postprocess -> metadata
[group('ingest')]
run: test download preprocess transform-all postprocess metadata
    @echo "Done!"

# Download source data
[group('ingest')]
download: install
    uv run downloader download.yaml

# Preprocess: convert VCF to TSV
[group('ingest')]
preprocess:
    uv run python scripts/vcf_to_tsv.py data/clinvar.vcf.gz data/clinvar.tsv

# Run all transforms
[group('ingest')]
transform-all: download
    #!/usr/bin/env bash
    set -euo pipefail
    # PYTHONPATH=src so koza's spec-loaded transform module can import
    # sibling files (clinvar_helpers, etc.) — same effect as pytest's
    # [tool.pytest.ini_options] pythonpath = ["src"].
    export PYTHONPATH=src
    for t in {{TRANSFORMS}}; do
        if [ -n "$t" ]; then
            echo "Transforming $t..."
            uv run koza transform {{PKG}}/$t.yaml
        fi
    done

# Emit output/release-metadata.yaml describing this build's upstream sources and artifacts
[group('ingest')]
metadata:
    uv run python scripts/write_metadata.py

# Run specific transform
[group('ingest')]
transform NAME:
    PYTHONPATH=src uv run koza transform {{PKG}}/{{NAME}}.yaml

# Postprocess (no-op for clinvar)
[group('ingest')]
postprocess:
    @echo "No postprocessing required"

# Run tests
[group('development')]
test: install
    uv run pytest

# Run tests with coverage
[group('development')]
test-cov: install
    uv run pytest --cov=. --cov-report=term-missing

# Lint code
[group('development')]
lint:
    uv run ruff check .

# Format code
[group('development')]
format:
    uv run ruff format .

# Clean output directory
[group('ingest')]
clean:
    rm -rf output/
