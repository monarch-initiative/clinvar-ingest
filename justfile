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

# Full pipeline: download -> preprocess -> transform -> postprocess
[group('ingest')]
run: download preprocess transform-all postprocess
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
    for t in {{TRANSFORMS}}; do
        if [ -n "$t" ]; then
            echo "Transforming $t..."
            uv run koza transform {{PKG}}/$t.yaml
        fi
    done

# Run specific transform
[group('ingest')]
transform NAME:
    uv run koza transform {{PKG}}/{{NAME}}.yaml

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
