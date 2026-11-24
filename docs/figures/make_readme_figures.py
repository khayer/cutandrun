#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_readme_figures.py
# Authors: Katharina E. Hayer, co-created with Claude
# Version: 1.5.0
# Date:    2026-11-23
# ---------------------------------------------------------------------------
# Renders the example figures used in readme_of_changes.md.
#
# EVERY NUMBER HERE IS SYNTHETIC. The figures exist to show the SHAPE of each
# output so a reader knows what to expect before running the pipeline; they are
# not results and must never be quoted as such. Each figure carries a
# "synthetic example data" stamp so a screenshot cannot be mistaken for a result.
#
# Each panel is drawn to the definition the pipeline actually implements, so the
# figure and the prose cannot drift apart:
#   - dual normalisation  : scale = (C / spikein_reads) * (mean_target / target)
#   - Top-N selection     : rank within loss/gain by raw p asc, |log2FC| desc,
#                           baseMean desc; keep first N
#   - random unchanged    : sample without replacement from status == unchanged,
#                           fixed seed 42
#   - Top 5 labels        : per class, rank by volcano metric asc then
#                           |log2FC| desc
#
# v1.1.0  Volcano fixes found by looking at the rendered PNG: the p-value floor
#         produced a fake flat ceiling of points, and all ten gene labels landed
#         on it on top of each other. p is now drawn from a continuous model with
#         no hard clip, and labels are staggered outward per class.
# v1.2.0  Added the pyGenomeTracks top-N snapshot figure. That panel deliberately
#         uses the MODULE'S OWN track colours (#1f78b4 own rep, #33a02c other
#         rep, #808080 control) rather than the validated palette, because its
#         job is to show what the real output looks like; a figure in different
#         colours would mislead. Those two signal colours were validated anyway
#         (deutan dE 23.4, normal dE 24.9, both >= 3:1 contrast) and every track
#         carries a title label, which covers the tritan pair at dE 6.7.
#
# Palette: the Weitzman Lab purple and green.
#   series purple #6b35a0, series green #2e7d4f, headings #5b2a86
# The house purple #5b2a86 is excellent as heading TEXT but sits just below the
# chart lightness band as a MARK (OKLab L 0.40 against a 0.43 floor), so marks
# use #6b35a0, the nearest passing step. That pair was validated as a set:
# CVD dE 19.8, normal-vision dE 27.4, both >= 3:1 contrast on this surface.
#
# The six genomic-feature classes cannot be told apart in two hues, so that one
# panel uses a purple SEQUENTIAL ramp ordered outward from the TSS (promoter,
# 5'UTR, exon, intron, 3'UTR, intergenic). The ramp is monotonic in lightness,
# and identity never rests on colour alone: every readable segment carries its
# own number and the legend is always present.
#
# v1.3.0  Restyled from the default categorical palette to the lab purple/green.
#         The pyGenomeTracks panel is the deliberate exception - it keeps the
#         module's own hardcoded track colours, because its job is to look like
#         the real output.
# v1.4.0  De-identified the sample labels. Earlier versions used the lab's real
#         construct and timepoint naming, which is unpublished experimental
#         design and this repository is public. The labels are now generic and
#         carry the same teaching load: two targets, a matched control per
#         timepoint, and one positive-control histone mark.
# v1.5.0  Replaced two matplotlib conveniences newer than the deeptools
#         container ships (3.2/3.3 on Python 3.8):
#           set_xticks(ticks, labels, ...) positionally   needs 3.5+
#           legend(..., labelcolor=...)                   needs 3.3+
#         Both now go through the setter API, stable since 2.x, so this runs on
#         the cluster as well as a modern desktop. Output is pixel-identical.
#
# Usage:  python3 make_readme_figures.py [outdir]
# ---------------------------------------------------------------------------

import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SURFACE   = "#fcfcfb"
INK       = "#0b0b0b"
INK2      = "#52514e"
MUTED     = "#898781"
GRID      = "#e1e0d9"
BASELINE  = "#c3c2b7"
PURPLE    = "#6b35a0"               # series purple (passing step of #5b2a86)
GREEN     = "#2e7d4f"               # series green
HEADING   = "#5b2a86"               # house purple, used as text only
# purple sequential ramp, light -> dark, for the six feature classes
CAT       = ["#ddd0ea", "#c6b0dc", "#ae8fcd", "#9570bd", "#7c50ab", "#6b35a0"]
DARK_TEXT_ON = {"#ddd0ea", "#c6b0dc", "#ae8fcd"}   # light steps need dark text
DIV_LOW, DIV_HIGH, NEUTRAL = PURPLE, GREEN, "#898781"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.size": 9,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": BASELINE,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "figure.dpi": 150,
})


def ticklabels(ax, axis, ticks, labels, **kw):
    """set_?ticks + set_?ticklabels, instead of the 3.5+ positional form."""
    getattr(ax, f"set_{axis}ticks")(list(ticks))
    getattr(ax, f"set_{axis}ticklabels")(list(labels), **kw)


def colour_legend(leg, colour):
    """Instead of legend(labelcolor=...), which needs matplotlib 3.3+."""
    for text in leg.get_texts():
        text.set_color(colour)
    return leg


def stamp(fig):
    fig.text(0.995, 0.005, "synthetic example data - not a result",
             ha="right", va="bottom", fontsize=7, color=MUTED, style="italic")


def finish(fig, path):
    stamp(fig)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", path)


# --- 1. peak feature annotation -------------------------------------------
def fig_feature_annotation(outdir):
    rng = np.random.default_rng(7)
    samples = ["Target_A_R1", "Target_A_R2", "Target_B_R1", "Target_B_R2", "ActiveMark"]
    feats = ["promoter", "5'UTR", "exon", "intron", "3'UTR", "intergenic"]
    base = np.array([
        [34, 5, 11, 26, 4, 20],
        [31, 6, 12, 27, 4, 20],
        [18, 4, 14, 33, 6, 25],
        [20, 4, 13, 32, 6, 25],
        [71, 4,  6, 11, 2,  6],      # an active mark is promoter-dominated, as expected
    ], dtype=float)
    base += rng.normal(0, 0.7, base.shape)
    pct = 100 * base / base.sum(axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    left = np.zeros(len(samples))
    y = np.arange(len(samples))
    for j, f in enumerate(feats):
        w = pct[:, j]
        ax.barh(y, w, left=left, height=0.62, color=CAT[j],
                edgecolor=SURFACE, linewidth=1.4, label=f)   # 2 px surface gap
        for i, (wi, li) in enumerate(zip(w, left)):
            if wi >= 7:                                       # required relief
                ax.text(li + wi / 2, y[i], f"{wi:.0f}",
                        ha="center", va="center", fontsize=7.5,
                        color=INK if CAT[j] in DARK_TEXT_ON else "#ffffff")
        left += w
    ticklabels(ax, 'y', y, samples, color=INK2)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("percent of peaks")
    ax.set_title("Peak feature annotation, per sample", color=HEADING, loc="left",
                 fontsize=11, pad=10)
    ax.grid(axis="y", visible=False)
    colour_legend(
        ax.legend(ncol=6, loc="upper left", bbox_to_anchor=(0, -0.22),
                  handlelength=1.1, columnspacing=1.2, fontsize=8), INK2)
    finish(fig, os.path.join(outdir, "fig_peak_feature_annotation.png"))


# --- 2. ChIPseeker distance to TSS ----------------------------------------
def fig_tss_distance(outdir):
    rng = np.random.default_rng(11)
    near = rng.normal(0, 900, 5200)
    far = rng.normal(0, 18000, 3000)
    d = np.concatenate([near, far])
    d = d[np.abs(d) <= 20000]

    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    ax.hist(d, bins=80, color=PURPLE, edgecolor=SURFACE, linewidth=0.4)
    ax.axvline(0, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.text(0, ax.get_ylim()[1] * 0.96, " TSS", ha="left", va="top",
            fontsize=8.5, color=INK2)
    for x in (-3000, 3000):
        ax.axvline(x, color=MUTED, linewidth=0.9, linestyle=(0, (1, 2)))
    ax.annotate("promoter window,\n--chipseeker_tss_dist 3000",
                xy=(3000, ax.get_ylim()[1] * 0.55), xytext=(9000, ax.get_ylim()[1] * 0.72),
                fontsize=8, color=INK2,
                arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.9))
    ax.set_xlabel("distance to nearest TSS (bp)")
    ax.set_ylabel("peaks")
    ax.set_title("Distance to TSS, all peaks", color=HEADING, loc="left",
                 fontsize=11, pad=10)
    ax.grid(axis="x", visible=False)
    finish(fig, os.path.join(outdir, "fig_chipseeker_tss_distance.png"))


# --- 3. dual normalisation -------------------------------------------------
def fig_dual_normalisation(outdir):
    samples = ["Control_T1", "Target_A_T1", "Control_T2", "Target_A_T2", "Target_B_T2"]
    spikein = np.array([42000, 38500, 51000, 47000, 29500], dtype=float)
    target = np.array([9.1e6, 8.7e6, 9.4e6, 9.0e6, 26.8e6])     # last: high viral load
    C = 10000.0
    single = C / spikein
    dual = single * (target.mean() / target)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.3))
    x = np.arange(len(samples))

    ax = axes[0]
    ax.bar(x, target / 1e6, width=0.58, color=PURPLE, edgecolor=SURFACE, linewidth=1.4)
    ax.axhline(target.mean() / 1e6, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.text(len(samples) - 0.5, target.mean() / 1e6, " mean", ha="left", va="center",
            fontsize=8, color=INK2)
    ticklabels(ax, 'x', x, samples, rotation=30, ha="right", color=INK2, fontsize=8)
    ax.set_ylabel("target-genome reads (millions)")
    ax.set_title("Input: target abundance varies", color=HEADING, loc="left", fontsize=10.5, pad=8)
    ax.grid(axis="x", visible=False)

    ax = axes[1]
    w = 0.33
    ax.bar(x - w / 2, single, width=w, color=PURPLE, edgecolor=SURFACE,
           linewidth=1.4, label="spike-in only")
    ax.bar(x + w / 2, dual, width=w, color=GREEN, edgecolor=SURFACE,
           linewidth=1.4, label="dual")
    ticklabels(ax, 'x', x, samples, rotation=30, ha="right", color=INK2, fontsize=8)
    ax.set_ylabel("scale factor")
    ax.set_title("Effect: dual down-weights the high-load sample",
                 color=HEADING, loc="left", fontsize=10.5, pad=8)
    ax.grid(axis="x", visible=False)
    ax.set_ylim(0, max(single.max(), dual.max()) * 1.32)   # headroom for legend
    colour_legend(
        ax.legend(loc="upper left", ncol=2, fontsize=8.5,
                  handlelength=1.1, columnspacing=1.4), INK2)

    fig.text(0.5, -0.13,
             "spike-in only:  scale = C / spikein_reads        "
             "dual:  scale = (C / spikein_reads) * (mean_target_reads / target_reads)",
             ha="center", fontsize=8, color=INK2, family="monospace")
    fig.tight_layout()
    finish(fig, os.path.join(outdir, "fig_dual_normalisation.png"))


# --- 4. promoter GC volcano, Top-N raw p-value mode ------------------------
def fig_volcano_topn(outdir):
    rng = np.random.default_rng(3)
    n = 4000
    log2fc = rng.normal(0, 0.75, n)
    log2fc[:250] -= rng.gamma(2.2, 0.8, 250)        # real losses
    log2fc[250:480] += rng.gamma(2.2, 0.8, 230)     # real gains
    # significance grows with effect size, with spread; no hard floor, so there
    # is no artificial ceiling of points at the clip value
    y = np.abs(log2fc) * rng.gamma(1.9, 0.85, n) + rng.exponential(0.35, n)
    pval = 10.0 ** (-y)
    base_mean = rng.lognormal(5, 1, n)
    rank = np.argsort(np.argsort(pval)) + 1
    padj = np.clip(pval * n / rank, 0, 1)
    status = np.where((padj < 0.05) & (np.abs(log2fc) > 1),
                      np.where(log2fc < 0, "loss", "gain"), "unchanged")

    TOP_N = 200                                     # 1000 by default; 200 reads better
    # rank: raw p asc, then |log2FC| desc, then baseMean desc  (lexsort: last key
    # is primary), keep the first N of each class
    order = np.lexsort((-base_mean, -np.abs(log2fc), pval))
    loss_rank = [i for i in order if log2fc[i] < 0][:TOP_N]
    gain_rank = [i for i in order if log2fc[i] > 0][:TOP_N]
    unch = np.flatnonzero(status == "unchanged")
    pick = np.random.default_rng(42).choice(
        unch, size=min(TOP_N, unch.size), replace=False)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.scatter(log2fc, y, s=6, color="#dcdbd5", linewidths=0, zorder=1,
               label="all promoter peaks")
    ax.scatter(log2fc[pick], y[pick], s=18, facecolors="none",
               edgecolors=NEUTRAL, linewidths=0.9, zorder=2,
               label=f"random unchanged (n={len(pick)}, seed 42)")
    ax.scatter(log2fc[loss_rank], y[loss_rank], s=20, color=DIV_LOW,
               edgecolors=SURFACE, linewidths=0.5, zorder=3,
               label=f"Top {TOP_N} loss")
    ax.scatter(log2fc[gain_rank], y[gain_rank], s=20, color=DIV_HIGH,
               edgecolors=SURFACE, linewidths=0.5, zorder=3,
               label=f"Top {TOP_N} gain")

    genes = {"loss": ["ADNP", "CHD4", "HP1BP3", "POGZ", "SMARCA4"],
             "gain": ["TRIM28", "ZNF143", "CTCF", "RAD51", "BRCA2"]}
    for cls, idxs, colr in (("loss", loss_rank, DIV_LOW), ("gain", gain_rank, DIV_HIGH)):
        lab = sorted(idxs, key=lambda i: (pval[i], -abs(log2fc[i])))[:5]
        lab = sorted(lab, key=lambda i: -y[i])        # label top-down
        ha = "right" if cls == "loss" else "left"
        dx = -7 if cls == "loss" else 7
        for k, i in enumerate(lab):
            ax.annotate(genes[cls][k], (log2fc[i], y[i]),
                        textcoords="offset points",
                        xytext=(dx, 13 - 6 * k),      # stagger so they cannot stack
                        ha=ha, va="center", fontsize=7.5, color=colr,
                        arrowprops=dict(arrowstyle="-", color=colr,
                                        linewidth=0.5, shrinkA=0, shrinkB=2))
    ax.axvline(0, color=BASELINE, linewidth=0.9)
    xpad = 0.20 * (log2fc.max() - log2fc.min())
    ax.set_xlim(log2fc.min() - xpad, log2fc.max() + xpad)
    ax.set_ylim(bottom=-0.3)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(raw p-value)")
    ax.set_title("Promoter GC differential binding, Top-N raw p-value mode",
                 color=HEADING, loc="left", fontsize=11, pad=10)
    colour_legend(
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2,
                  fontsize=8, handlelength=1.1, scatterpoints=1), INK2)
    finish(fig, os.path.join(outdir, "fig_promoter_gc_volcano_topn.png"))


# --- 5. pyGenomeTracks top-N snapshot --------------------------------------
PYGT_OWN, PYGT_OTHER, PYGT_CTRL = "#1f78b4", "#33a02c", "#808080"


def fig_pygenometracks(outdir):
    """Mimics the track stack the module writes: x-axis on top, a GTF gene
    track, the consensus-peak BED, then one bigWig signal track per sample."""
    rng = np.random.default_rng(19)
    start, end = 41_200, 49_800                       # a 2 kb peak +/- flank
    x = np.linspace(start, end, 1400)
    peak_c, peak_w = 45_500, 520

    def signal(amp, noise, extra=0.0):
        y = amp * np.exp(-0.5 * ((x - peak_c) / peak_w) ** 2)
        y += extra * np.exp(-0.5 * ((x - 43_300) / 700) ** 2)
        y += np.abs(rng.normal(0, noise, x.size))
        return np.convolve(y, np.ones(9) / 9, mode="same")

    tracks = [
        ("Target_A_R1", signal(26, 0.9), PYGT_OWN),
        ("Target_A_R2", signal(23, 0.9), PYGT_OWN),
        ("Target_A_R3 (other rep)", signal(19, 1.0, 3), PYGT_OTHER),
        ("Control_R1 (control rep)", signal(2.2, 0.8), PYGT_CTRL),
    ]
    heights = [1.0, 3, 1.2] + [1.5] * len(tracks)     # module's own heights
    fig, axes = plt.subplots(len(heights), 1, figsize=(7.4, 5.4),
                             gridspec_kw={"height_ratios": heights, "hspace": 0.22})
    for ax in axes:
        ax.set_xlim(start, end)
        ax.grid(False)
        for sp in ("top", "right", "left", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])

    # x-axis, drawn at the top as the module sets it (where = top), with the
    # labels ABOVE the spine so the panel below cannot clip them
    axt = axes[0]
    axt.spines["top"].set_visible(True)
    axt.spines["top"].set_color(BASELINE)
    axt.xaxis.set_ticks_position("top")
    ticks = np.linspace(start, end, 5)
    ticklabels(axt, 'x', ticks, [f"{t/1000:.1f} kb" for t in ticks], fontsize=7.5)
    axt.tick_params(axis="x", length=3, pad=2, labeltop=True, labelbottom=False)
    for lab, ha in ((axt.get_xticklabels()[0], "left"),
                    (axt.get_xticklabels()[-1], "right")):
        lab.set_horizontalalignment(ha)

    # gene track, flybase style: thick exons on a thin intron line with arrows
    axg = axes[1]
    axg.set_ylim(0, 1)
    for gstart, gend, name, strand in [(41_600, 44_900, "UL54", "+"),
                                       (46_100, 49_400, "UL55", "-")]:
        axg.plot([gstart, gend], [0.55, 0.55], color=INK, linewidth=0.9, zorder=1)
        exons = np.linspace(gstart, gend, 7)
        for a, b in zip(exons[::2], exons[1::2]):
            axg.add_patch(plt.Rectangle((a, 0.38), b - a, 0.34, color=INK, zorder=2))
        for xa in np.linspace(gstart, gend, 9)[1:-1]:
            axg.annotate("", xy=(xa + (160 if strand == "+" else -160), 0.55),
                         xytext=(xa, 0.55),
                         arrowprops=dict(arrowstyle="->", color=INK, linewidth=0.7))
        axg.text((gstart + gend) / 2, 0.14, name, ha="center", fontsize=8, color=INK)
    axg.set_ylabel("Genes", rotation=0, ha="right", va="center",
                   fontsize=8.5, color=INK2, labelpad=8)

    # consensus peak track
    axp = axes[2]
    axp.set_ylim(0, 1)
    axp.add_patch(plt.Rectangle((peak_c - 1000, 0.3), 2000, 0.4, color=INK))
    axp.set_ylabel("Consensus\npeaks", rotation=0, ha="right", va="center",
                   fontsize=8.5, color=INK2, labelpad=8)

    ymax = max(t[1].max() for t in tracks) * 1.12
    for ax, (name, y, colr) in zip(axes[3:], tracks):
        ax.fill_between(x, 0, y, color=colr, linewidth=0)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(name, rotation=0, ha="right", va="center",
                      fontsize=8, color=INK2, labelpad=8)
        ax.text(start + 120, ymax * 0.72, f"[0 - {ymax:.0f}]",
                fontsize=6.5, color=MUTED)

    fig.suptitle("Top consensus peak 1 of 10  |  nearest gene UL54, 310 bp from TSS",
                 x=0.125, ha="left", fontsize=10.5, color=HEADING, y=0.955)
    finish(fig, os.path.join(outdir, "fig_pygenometracks_top10.png"))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)
    print("rendering synthetic example figures into", out)
    fig_feature_annotation(out)
    fig_tss_distance(out)
    fig_dual_normalisation(out)
    fig_volcano_topn(out)
    fig_pygenometracks(out)
    print("done")
