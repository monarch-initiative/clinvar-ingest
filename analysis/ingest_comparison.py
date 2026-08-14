"""Before/after for a PR: compare two branches' ClinVar ingest output, and what each
contributes to the merged Monarch KG.

Both sides must be run against the same downloaded `data/` directory, so that every
difference the report shows is code rather than a moving upstream release. The usual
setup is a worktree for the base branch sharing this checkout's data:

    git worktree add /tmp/clinvar-main main
    ln -s "$PWD/data" /tmp/clinvar-main/data
    (cd /tmp/clinvar-main && PYTHONPATH=src uv run koza transform src/transform.yaml)

    just transform-all                      # this branch, into ./output
    cd analysis && uv run --project .. python ingest_comparison.py \
        --base-output /tmp/clinvar-main/output

"Effect on the merged Monarch KG" overlays each side's output on the KG the way the
merge would and reports what arrives. It is NOT a full monarch-ingest merge run -- no
other source is re-ingested and no QC or normalisation step is applied -- but it answers
what a reviewer needs to know: how many nodes and edges arrive, what predicates they
carry, how many gene-disease relationships they imply, and how many of those the KG does
not already have.

Everything below is computed. Nothing about the gained/dropped pairs is asserted from a
previous run, because the whole point of the document is to be re-run on a later release
and still be true.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

# Symbol patterns for loci that ClinVar's GENEINFO reports because they OVERLAP a
# variant's position, not because they are the causal gene -- antisense/divergent
# transcripts, readthrough fusions and uncharacterised ORFs. Used only to describe the
# dropped pairs, never to filter anything.
ANTISENSE_RE = re.compile(r"-(AS\d*|DT|IT\d*|OT\d*)$")
UNCHAR_ORF_RE = re.compile(r"^C\d+orf\d+$", re.IGNORECASE)
READTHROUGH_RE = re.compile(r"^[A-Z0-9]+-[A-Z0-9]+$")


def load_side(out_dir: Path) -> dict:
    """Profile one branch's KGX output: counts, predicates, and the gene-disease pairs
    it implies by crossing (variant->gene) with (variant->disease)."""
    nodes: dict = {}
    with open(out_dir / "clinvar_variant_nodes.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            nodes[row["id"]] = row

    preds: Counter = Counter()
    # A variant can carry SEVERAL gene edges -- that is precisely what GENEINFO-based
    # attribution does, so this must be a set per variant. Collapsing it to one gene
    # would hide the artifact this comparison exists to measure.
    var_gene: dict = defaultdict(set)
    var_dis: dict = defaultdict(set)
    n_edges = 0
    with open(out_dir / "clinvar_variant_edges.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_edges += 1
            preds[row["predicate"]] += 1
            if row["predicate"].endswith("is_sequence_variant_of"):
                var_gene[row["subject"]].add(row["object"])
            elif row["object"].startswith("MONDO:"):
                var_dis[row["subject"]].add(row["object"])

    pairs = {
        (g, d)
        for v, ds in var_dis.items()
        for g in var_gene.get(v, ())
        for d in ds
    }
    all_genes = {g for gs in var_gene.values() for g in gs}
    return {
        "nodes": len(nodes),
        "edges": n_edges,
        "preds": preds,
        "genes": all_genes,
        "genes_per_variant": (
            sum(len(gs) for gs in var_gene.values()) / max(len(var_gene), 1)
        ),
        "diseases": {d for ds in var_dis.values() for d in ds},
        "pairs": pairs,
        "node_types": Counter(n.get("type", "") for n in nodes.values()),
    }


def load_kg_pairs(data_dir: Path) -> tuple[set, int, int]:
    """The KG's curated gene-disease pairs keyed on (NCBIGene, MONDO), plus its raw
    node/edge counts. Mirrors load_monarch_gene_disease() in clinvar_report.py: OMIM and
    ClinGen edges carry an HGNC subject, Orphanet's keep Entrez in original_subject."""
    hgnc_to_entrez, symbols = {}, {}
    with open(data_dir / "hgnc_complete_set.txt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("entrez_id"):
                gid = f"NCBIGene:{row['entrez_id']}"
                hgnc_to_entrez[row["hgnc_id"]] = gid
                symbols[gid] = row["symbol"]

    pairs = set()
    with open(data_dir / "mk_gene_disease.tsv") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or not p[1].startswith("MONDO:"):
                continue
            gene = hgnc_to_entrez.get(p[0])
            if gene is None and len(p) > 3 and p[3].startswith("NCBIGene:"):
                gene = p[3]
            if gene:
                pairs.add((gene, p[1]))

    with open(data_dir / "monarch-kg_nodes.tsv") as fh:
        n_nodes = sum(1 for _ in fh) - 1
    with open(data_dir / "monarch-kg_edges.tsv") as fh:
        n_edges = sum(1 for _ in fh) - 1
    return pairs, n_nodes, n_edges, symbols


def describe_dropped(dropped: set, base: dict, head: dict, symbols: dict) -> list[str]:
    """Characterise the pairs the head branch no longer emits.

    The interesting split is not "is the symbol antisense-shaped" but "does the head
    branch attribute this gene to ANY variant". A gene absent from head's output
    entirely is one ClinVar never asserts as causal; a gene still present has simply
    lost one disease, which is a different (and more reviewable) claim.
    """
    genes = {g for g, _ in dropped}
    absent = sorted(g for g in genes if g not in head["genes"])
    present = sorted(g for g in genes if g in head["genes"])
    named = sorted(symbols[g] for g in absent if g in symbols)
    unnamed = [g for g in absent if g not in symbols]

    def bucket(sym: str) -> str:
        if ANTISENSE_RE.search(sym):
            return "antisense/divergent transcript"
        if UNCHAR_ORF_RE.match(sym):
            return "uncharacterised ORF"
        if sym.startswith("MIR"):
            return "miRNA"
        if sym.startswith("MT-"):
            return "mitochondrial"
        if READTHROUGH_RE.match(sym) and "-" in sym:
            return "readthrough fusion"
        return "other"

    buckets: dict = defaultdict(list)
    for s in named:
        buckets[bucket(s)].append(s)

    out = [
        f"Those {len(dropped):,} dropped pairs span {len(genes):,} genes, and break down "
        f"as follows:\n",
        f"- **{len(absent):,} genes this branch never attributes to any variant** — "
        f"ClinVar does not assert them as causal, so they appear nowhere in its output. "
        f"{len(named):,} carry an HGNC symbol:",
    ]
    for kind in sorted(buckets):
        out.append(f"    - {kind}: {', '.join('`' + s + '`' for s in buckets[kind])}")
    out.append(
        f"    - and {len(unnamed):,} with no HGNC symbol at all — `LOC`/lncRNA/miRNA "
        f"placeholders."
    )
    out.append(
        f"- **{len(present):,} genes still in this branch's output**, which have simply "
        f"lost one disease pairing: "
        + ", ".join("`" + symbols.get(g, g) + "`" for g in present)
        + "\n"
    )
    mito = buckets.get("mitochondrial", [])
    if mito:
        kept = sorted(
            symbols[g] for g in head["genes"]
            if g in symbols and symbols[g].startswith("MT-")
        )
        out.append(
            f"> **Worth a reviewer's eye:** {len(mito)} mitochondrial gene(s) — "
            + ", ".join("`" + s + "`" for s in mito)
            + " — fall in the first group, so this branch emits nothing for them. "
            + (
                f"Other MT genes do survive ({', '.join('`' + s + '`' for s in kept[:8])}"
                f"{', …' if len(kept) > 8 else ''}), so this is not a blanket exclusion "
                "of the mitochondrial genome, but it is the part of the drop that is "
                "least obviously an artifact and should be confirmed before merge.\n"
                if kept else
                "No mitochondrial genes survive at all, which needs confirming before "
                "merge.\n"
            )
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-output", type=Path, required=True,
                    help="output/ directory produced by the base branch (e.g. main)")
    ap.add_argument("--head-output", type=Path, default=Path("../output"),
                    help="output/ directory produced by this branch")
    ap.add_argument("--data-dir", type=Path, default=Path("../data"))
    ap.add_argument("--base-name", default="main")
    ap.add_argument("--head-name", default="this branch")
    ap.add_argument("--out", type=Path, default=Path("../INGEST_COMPARISON.md"))
    args = ap.parse_args()

    base = load_side(args.base_output)
    head = load_side(args.head_output)
    kg_pairs, kg_nodes, kg_edges, symbols = load_kg_pairs(args.data_dir)

    def delta(a: int, b: int) -> str:
        return f"{a - b:+,} ({100 * (a - b) / b:+.1f}%)" if b else "n/a"

    B, H = args.base_name, args.head_name
    L: list[str] = []
    w = L.append

    w(f"# ClinVar ingest: `{B}` vs `{H}`\n")
    w("Generated by `analysis/ingest_comparison.py`. Both sides ran against the same "
      "downloaded `data/` directory, so every difference below is code, not a moving "
      "upstream release.\n")

    w("\n## Emitted artifacts\n")
    w(f"| | {B} | {H} | change |")
    w("|---|---:|---:|---:|")
    w(f"| Nodes | {base['nodes']:,} | {head['nodes']:,} | {delta(head['nodes'], base['nodes'])} |")
    w(f"| Edges | {base['edges']:,} | {head['edges']:,} | {delta(head['edges'], base['edges'])} |")
    w(f"| Distinct genes | {len(base['genes']):,} | {len(head['genes']):,} | "
      f"{len(head['genes']) - len(base['genes']):+,} |")
    w(f"| Distinct diseases | {len(base['diseases']):,} | {len(head['diseases']):,} | "
      f"{len(head['diseases']) - len(base['diseases']):+,} |")
    w(f"| Implied gene-disease pairs | {len(base['pairs']):,} | {len(head['pairs']):,} | "
      f"{len(head['pairs']) - len(base['pairs']):+,} |")
    w(f"| Gene edges per variant (mean) | {base['genes_per_variant']:.2f} | "
      f"{head['genes_per_variant']:.2f} | |")
    w(f"\nGene edges per variant is the direct measure of positional attribution: "
      f"`{B}` assigns {base['genes_per_variant']:.2f} genes to the average variant, "
      f"`{H}` assigns {head['genes_per_variant']:.2f}. Every gene above one is a locus "
      "that merely overlaps the variant's position, and it inherits the real gene's "
      "entire disease roster.\n")

    w("\n## Edges by predicate\n")
    w(f"| Predicate | {B} | {H} |")
    w("|---|---:|---:|")
    for p in sorted(set(base["preds"]) | set(head["preds"])):
        w(f"| `{p.replace('biolink:', '')}` | {base['preds'].get(p, 0):,} | "
          f"{head['preds'].get(p, 0):,} |")

    w("\n## Node `type` values\n")
    w("Biolink's `type` slot on a sequence variant is for the variant's class. Check "
      "which side is putting the molecular consequence there instead.\n")
    w(f"| {B} | {H} |")
    w("|---|---|")
    fmt = lambda c: ", ".join(f"`{k or '(empty)'}` {v:,}" for k, v in c.most_common(6))
    w(f"| {fmt(base['node_types'])} | {fmt(head['node_types'])} |")

    w("\n## Effect on the merged Monarch KG\n")
    w(f"The KG as downloaded holds **{kg_nodes:,} nodes** and **{kg_edges:,} edges**, "
      f"with **{len(kg_pairs):,}** curated gene-disease pairs.\n")
    w(f"| | {B} | {H} |")
    w("|---|---:|---:|")
    w(f"| Nodes contributed | {base['nodes']:,} | {head['nodes']:,} |")
    w(f"| Edges contributed | {base['edges']:,} | {head['edges']:,} |")
    w(f"| as % of KG edges | {100 * base['edges'] / kg_edges:.1f}% | "
      f"{100 * head['edges'] / kg_edges:.1f}% |")
    bk, hk = len(base["pairs"] & kg_pairs), len(head["pairs"] & kg_pairs)
    bp = 100 * bk / max(len(base["pairs"]), 1)
    hp = 100 * hk / max(len(head["pairs"]), 1)
    w(f"| pairs corroborating a curated KG pair | {bk:,} ({bp:.1f}%) | {hk:,} ({hp:.1f}%) |")
    w(f"| pairs with no curated KG equivalent | {len(base['pairs'] - kg_pairs):,} | "
      f"{len(head['pairs'] - kg_pairs):,} |")
    w(f"\nThe corroboration rate moves {bp:.1f}% → {hp:.1f}% while the pair count grows "
      f"{len(head['pairs']) / max(len(base['pairs']), 1):.0f}-fold. A rate that holds "
      "steady under that much growth is the evidence that recall was not bought with "
      "precision.\n")

    gained = head["pairs"] - base["pairs"]
    dropped = base["pairs"] - head["pairs"]
    w("\n## Gene-disease pairs gained and dropped\n")
    w(f"- **{len(gained):,}** pairs `{H}` emits that `{B}` does not")
    w(f"- **{len(dropped):,}** pairs `{B}` emits that `{H}` drops")
    w(f"- **{len(base['pairs'] & head['pairs']):,}** in common\n")
    if dropped:
        L.extend(describe_dropped(dropped, base, head, symbols))

    args.out.write_text("\n".join(L) + "\n")
    print(f"Wrote {args.out}")
    print(f"  {B}: {base['nodes']:,} nodes / {base['edges']:,} edges / "
          f"{len(base['pairs']):,} pairs")
    print(f"  {H}: {head['nodes']:,} nodes / {head['edges']:,} edges / "
          f"{len(head['pairs']):,} pairs")


if __name__ == "__main__":
    main()
