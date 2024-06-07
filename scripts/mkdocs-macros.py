import csv
from pathlib import Path

nodes_report_file = Path("docs/nodes_report.tsv")
edges_report_file = Path("docs/edges_report.tsv")


def define_env(env):
    """Define mkdocs macros."""

    @env.macro
    def get_nodes_report() -> str:
        if not nodes_report_file.exists():
            return ""

        report = "## Nodes Report\n\n"
        with open(nodes_report_file, newline="\n") as csvfile:
            reader = csv.reader(csvfile, delimiter="\t")

            # turn into markdown table
            headers = next(reader)
            report += "|" + "|".join(headers) + "|\n"
            report += "|" + "|".join(["---" for _ in range(0, len(headers))]) + "|\n"
            for row in reader:
                report += "|" + "|".join(row) + "|\n"

        return report

    @env.macro
    def get_edges_report() -> str:
        if not edges_report_file.exists():
            return ""

        report = "## Edges Report\n\n"
        with open(edges_report_file, newline="\n") as csvfile:
            reader = csv.reader(csvfile, delimiter="\t")

            # turn into markdown table
            headers = next(reader)
            report += "|" + "|".join(headers) + "|\n"
            report += "|" + "|".join(["---" for _ in range(0, len(headers))]) + "|\n"
            for row in reader:
                report += "|" + "|".join(row) + "|\n"

        return report
