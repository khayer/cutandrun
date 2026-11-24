#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# species_composition.py
# Authors: Katharina E. Hayer, co-created with Claude
# Version: 1.0.0
# Date:    2026-11-23
# ---------------------------------------------------------------------------
# Aggregates samtools idxstats across samples into per-species read counts and
# percentages, and draws them.
#
# WHY TWO PANELS. A 100% stacked bar answers "how much is host vs virus", which
# is the normalisation question. It CANNOT answer "is anything contaminating
# this library", because 0.4% mycoplasma and 0.0% mycoplasma look identical at
# that scale. The second panel plots non-host counts on a log axis, where a
# trace is visible. Both are needed; neither alone is enough.
#
# Colour follows the ROLE, not an arbitrary identity: host is the lab purple,
# viral contigs take the green ramp, and decoys/contaminants are grey so a
# non-zero grey bar reads as "wrong" at a glance.
# ---------------------------------------------------------------------------

import argparse, os, re, sys
from collections import OrderedDict

HOST_PURPLE = "#6b35a0"
VIRAL_GREENS = ["#2e7d4f", "#4ca06f", "#7cc49a", "#a9dcc0"]
DECOY_GREYS = ["#8a867e", "#b0aca4", "#c9c6bf", "#ddd9d2"]
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
HEADING = "#5b2a86"


def parse_groups(spec):
    """'human=^(chr|GL|KI);mycoplasma=^Myco_' -> [(label, compiled regex), ...]"""
    out = []
    for part in spec.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            sys.exit(f"bad --groups entry (need label=regex): {part!r}")
        label, pattern = part.split("=", 1)
        out.append((label.strip(), re.compile(pattern.strip())))
    return out


def species_of(contig, groups):
    for label, rx in groups:
        if rx.search(contig):
            return label
    return contig


def sample_name(path):
    base = os.path.basename(path)
    for suffix in (".idxstats", ".target.filtered", ".target.sorted", ".target",
                   ".spikein.filtered", ".spikein.sorted", ".spikein", ".filtered",
                   ".sorted"):
        base = base.replace(suffix, "")
    return base


def read_idxstats(path, groups):
    counts = OrderedDict()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] == "*":
                continue
            counts[species_of(f[0], groups)] = counts.get(species_of(f[0], groups), 0) + int(f[2])
    return counts


def classify(label, host_label, viral_labels):
    if label == host_label:
        return "host"
    return "viral" if label in viral_labels else "decoy"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("idxstats", nargs="+")
    ap.add_argument("--groups", default="human=^(chr|GL|KI);mycoplasma=^Myco_")
    ap.add_argument("--host", default="human")
    ap.add_argument("--viral", default="Ad5,HSV17_NR",
                    help="comma-separated species treated as viral (green)")
    ap.add_argument("--prefix", default="species_composition")
    ap.add_argument("--title", default="Reads per species")
    a = ap.parse_args()

    groups = parse_groups(a.groups)
    viral = {v.strip() for v in a.viral.split(",") if v.strip()}

    per_sample = OrderedDict()
    for p in sorted(a.idxstats):
        per_sample[sample_name(p)] = read_idxstats(p, groups)
    if not per_sample:
        sys.exit("no idxstats files read")

    # stable species order: host, viral, then decoys, each by total descending
    totals = {}
    for counts in per_sample.values():
        for sp, n in counts.items():
            totals[sp] = totals.get(sp, 0) + n
    order = sorted(totals, key=lambda s: (
        {"host": 0, "viral": 1, "decoy": 2}[classify(s, a.host, viral)], -totals[s]))

    # ---- tables ----------------------------------------------------------
    with open(f"{a.prefix}.tsv", "w") as fh:
        fh.write("sample\t" + "\t".join(order) + "\ttotal\n")
        for s, c in per_sample.items():
            row = [c.get(sp, 0) for sp in order]
            fh.write(f"{s}\t" + "\t".join(str(v) for v in row) + f"\t{sum(row)}\n")

    # MultiQC custom content: it renders this as a stacked bargraph on its own
    with open(f"{a.prefix}_mqc.tsv", "w") as fh:
        fh.write("# id: 'species_composition'\n")
        fh.write("# section_name: 'Species composition'\n")
        fh.write("# description: 'Reads assigned to each species in the combined "
                 "reference, from samtools idxstats. Decoy contigs should be ~0.'\n")
        fh.write("# plot_type: 'bargraph'\n")
        fh.write("# pconfig:\n")
        fh.write("#    id: 'species_composition_plot'\n")
        fh.write("#    title: 'Reads per species'\n")
        fh.write("#    ylab: 'Reads'\n")
        fh.write("Sample\t" + "\t".join(order) + "\n")
        for s, c in per_sample.items():
            fh.write(f"{s}\t" + "\t".join(str(c.get(sp, 0)) for sp in order) + "\n")

    # ---- figure ----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "font.size": 9, "text.color": INK,
        "axes.labelcolor": INK2, "axes.edgecolor": "#c3c2b7",
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": 150,
    })

    colour, ng, nd = {}, 0, 0
    for sp in order:
        kind = classify(sp, a.host, viral)
        if kind == "host":
            colour[sp] = HOST_PURPLE
        elif kind == "viral":
            colour[sp] = VIRAL_GREENS[min(ng, len(VIRAL_GREENS) - 1)]; ng += 1
        else:
            colour[sp] = DECOY_GREYS[min(nd, len(DECOY_GREYS) - 1)]; nd += 1

    samples = list(per_sample)
    y = np.arange(len(samples))
    h = max(2.6, 0.42 * len(samples) + 1.5)
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.2, h),
                                   gridspec_kw={"width_ratios": [1.55, 1]})

    left = np.zeros(len(samples))
    for sp in order:
        vals = np.array([per_sample[s].get(sp, 0) for s in samples], dtype=float)
        tot = np.array([max(1, sum(per_sample[s].values())) for s in samples], dtype=float)
        pct = 100 * vals / tot
        axl.barh(y, pct, left=left, height=0.62, color=colour[sp],
                 edgecolor=SURFACE, linewidth=1.4, label=sp)
        for i, (v, l) in enumerate(zip(pct, left)):
            if v >= 6:
                axl.text(l + v / 2, y[i], f"{v:.1f}", ha="center", va="center",
                         fontsize=7.5, color="#ffffff")
        left += pct
    axl.set_yticks(y, samples, color=INK2, fontsize=8)
    axl.invert_yaxis(); axl.set_xlim(0, 100)
    axl.set_xlabel("percent of mapped reads")
    axl.set_title("Composition", color=HEADING, loc="left", fontsize=10.5, pad=8)
    axl.grid(axis="y", visible=False)

    # right panel: everything that is not host, on a log axis
    nonhost = [sp for sp in order if sp != a.host]
    width = 0.8 / max(1, len(nonhost))
    for k, sp in enumerate(nonhost):
        vals = np.array([max(per_sample[s].get(sp, 0), 0.5) for s in samples], dtype=float)
        axr.barh(y + (k - (len(nonhost) - 1) / 2) * width, vals, height=width,
                 color=colour[sp], edgecolor=SURFACE, linewidth=0.5, label=sp)
    axr.set_xscale("log")
    axr.set_xlim(0.5, None)
    axr.set_yticks(y, ["" for _ in samples])
    axr.invert_yaxis()
    axr.set_xlabel("non-host reads (log scale; 0 drawn at 0.5)")
    axr.set_title("Non-host detail", color=HEADING, loc="left", fontsize=10.5, pad=8)
    axr.grid(axis="y", visible=False)

    handles, labels = axl.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=min(len(order), 6),
               bbox_to_anchor=(0.5, -0.14), labelcolor=INK2, fontsize=8,
               handlelength=1.1, columnspacing=1.4)
    fig.suptitle(a.title, x=0.045, ha="left", fontsize=12, color=HEADING, y=1.0)
    fig.subplots_adjust(bottom=0.26)
    fig.savefig(f"{a.prefix}.png", bbox_inches="tight")
    print(f"wrote {a.prefix}.tsv, {a.prefix}_mqc.tsv, {a.prefix}.png "
          f"({len(samples)} samples, {len(order)} species)")


if __name__ == "__main__":
    main()
