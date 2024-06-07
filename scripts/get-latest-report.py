import json

import requests


def main():
    url = "https://api.github.com/repos/monarch-initiative/clinvar-ingest/releases/latest"

    # Get the latest release from the GitHub API
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(
            f"\n\tFailed to get latest release from {url}\n\tStatus: {response.status_code} - {response.text}"
        )
    data = json.loads(response.text)

    # Get the download URLs for the reports
    reports = {}
    for asset in data["assets"]:
        report_name = asset["name"]
        if "report.tsv" in asset["name"].split("_"):
            file_url = asset["browser_download_url"]
            reports[report_name] = file_url

    if not reports:
        raise Exception("No reports found in the latest release")

    # Download the reports
    for fn, url in reports.items():
        response = requests.get(url)
        output_fn = "_".join(fn.split("_")[-2:])
        with open(f"docs/{output_fn}", "wb") as f:
            f.write(response.content)


if __name__ == "__main__":
    main()

