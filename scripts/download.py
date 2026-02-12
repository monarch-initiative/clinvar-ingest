#!/usr/bin/env python3
"""Download source data files specified in download.yaml."""

import urllib.request
from pathlib import Path

import yaml


def download_files():
    """Download all files specified in download.yaml."""
    with open("download.yaml") as f:
        config = yaml.safe_load(f)

    downloads = config.get("downloads", [])
    if not downloads:
        print("No downloads configured in download.yaml")
        return

    for item in downloads:
        url = item["url"]
        local_name = item["local_name"]
        local_path = Path(local_name)

        # Create parent directories if needed
        local_path.parent.mkdir(parents=True, exist_ok=True)

        if local_path.exists():
            print(f"Skipping {local_name} (already exists)")
            continue

        print(f"Downloading {url} -> {local_name}")
        urllib.request.urlretrieve(url, local_path)
        print(f"  Downloaded {local_path.stat().st_size:,} bytes")


if __name__ == "__main__":
    download_files()
