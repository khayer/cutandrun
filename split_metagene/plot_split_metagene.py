#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# plot_split_metagene.py
# Authors : Katharina (Nina) Hayer, co-created with Claude Fable 5.1 (Anthropic)
# Version : 1.2.0
# Date    : 2026-10-06
# Purpose : Draw a split metagene figure from two deepTools computeMatrix
#           outputs (TSS-anchored and TES-anchored): mean profiles on top,
#           row-grouped heatmaps below, one shared row order across panels.
#           Layout mirrors the "intron-containing vs intronless" panel style.
# Fonts   : Calibri (falls back to Carlito/DejaVu); nothing below 12 pt.
# -----------------------------------------------------------------------------
import argparse
import gzip
import json

import logging

import numpy as np
import matplotlib
matplotlib.use("Agg")
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, FancyArrow

# ---- house style -------------------------------------------------------------
PURPLE_DARK = "#2E2A7A"   # intron-containing line
GREEN_LIGHT = "#A7D86B"   # intronless line
HEAT_CMAP = LinearSegmentedColormap.from_list(
    "whitepurple", ["#F7F6FB", "#C9BEE6", "#7A5BC2", "#3B1E7A"])
plt.rcParams.update({
    "font.family": ["Calibri", "Carlito", "DejaVu Sans"],
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 14,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    "axes.linewidth": 1.6, "pdf.fonttype": 42,
})


def read_matrix(path):
    """Return (header, values[region x bin] averaged over samples, names, n_samples).

    computeMatrix lays samples side by side in each row (sample_boundaries);
    replicates are averaged per region/bin so the figure shows the condition.
    """
    with gzip.open(path, "rt") as fh:
        first = fh.readline()
        assert first.startswith("@"), "not a computeMatrix file"
        hdr = json.loads(first[1:])
        rows, names = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            names.append(f[3])
            rows.append(np.array([np.nan if v in ("nan", "") else float(v)
                                  for v in f[6:]]))
    m = np.vstack(rows)
    sb = hdr["sample_boundaries"]
    blocks = [m[:, sb[i]:sb[i + 1]] for i in range(len(sb) - 1)]
    blocks = np.stack(blocks, axis=0)          # rep x region x bin
    avg = np.nanmean(blocks, axis=0)
    return hdr, avg, names, blocks


def groups_from_header(hdr, n_rows):
    b = hdr["group_boundaries"]
    labels = hdr.get("group_labels", [f"group{i}" for i in range(len(b) - 1)])
    return [(labels[i], b[i], b[i + 1]) for i in range(len(b) - 1)]


def nice_label(lbl):
    l = lbl.lower()
    if "intron_containing" in l or "intron-containing" in l:
        return "intron-\ncontaining"
    if "intronless" in l:
        return "intronless"
    return lbl.replace("_", " ")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tss", required=True, help="computeMatrix TSS gz")
    ap.add_argument("--tes", required=True, help="computeMatrix TES gz")
    ap.add_argument("--title", default="")
    ap.add_argument("--out", default="split_metagene")
    ap.add_argument("--vmax", type=float, default=None,
                    help="heatmap max; default 98th percentile")
    ap.add_argument("--ymax", type=float, default=None, help="profile y max")
    ap.add_argument("--ylabel", default=None, help="y label; default by --signal")
    ap.add_argument("--sort-by", choices=["tss", "tes", "mean"], default="mean",
                    help="order rows within a group by signal in this panel")
    ap.add_argument("--ci", choices=["95", "sem", "none"], default="95",
                    help="band around profile lines when >1 replicate: 95%% "
                         "t-interval across replicates, +/-1 SEM, or none")
    ap.add_argument("--alpha", type=float, default=0.8,
                    help="opacity of the profile lines (0-1)")
    ap.add_argument("--signal", choices=["subtract", "log2ratio"], default="subtract",
                    help="subtract: >=0 scale, white->purple; "
                         "log2ratio: symmetric scale, blue<-0->purple")
    args = ap.parse_args()

    hT, mT, nT, rT = read_matrix(args.tss)
    hE, mE, nE, rE = read_matrix(args.tes)
    assert nT == nE, "TSS and TES matrices must contain the same regions in order"
    nrep = rT.shape[0]
    if args.signal == "subtract":
        mT, mE = np.clip(mT, 0, None), np.clip(mE, 0, None)
        rT, rE = np.clip(rT, 0, None), np.clip(rE, 0, None)
    groups = groups_from_header(hT, len(nT))

    upT, dnT = hT["upstream"][0], hT["downstream"][0]
    upE, dnE = hE["upstream"][0], hE["downstream"][0]
    binsz = hT["bin size"][0]

    # ---- shared row order: within each group, sort descending by signal ------
    order = []
    for _, a, b in groups:
        if args.sort_by == "tss":
            key = np.nanmean(mT[a:b], axis=1)
        elif args.sort_by == "tes":
            key = np.nanmean(mE[a:b], axis=1)
        else:
            key = np.nanmean(np.hstack([mT[a:b], mE[a:b]]), axis=1)
        order.extend(a + np.argsort(-np.nan_to_num(key)))
    order = np.array(order)
    mT, mE = mT[order], mE[order]
    rT, rE = rT[:, order], rE[:, order]

    allv = np.concatenate([mT.ravel(), mE.ravel()])
    if args.signal == "log2ratio":
        vmax = args.vmax if args.vmax else np.nanpercentile(np.abs(allv), 98)
        vmin, cmap = -vmax, LinearSegmentedColormap.from_list(
            "bluewhitepurple", ["#1B5E9E", "#F7F6FB", "#3B1E7A"])
    else:
        vmax = args.vmax if args.vmax else np.nanpercentile(allv, 98)
        vmin, cmap = 0.0, HEAT_CMAP

    # ---- figure ---------------------------------------------------------------
    n_total = len(order)
    heights = [b - a for _, a, b in groups]
    fig = plt.figure(figsize=(7.2, 7.6))
    outer = gridspec.GridSpec(
        3, 3, figure=fig, height_ratios=[1.4, 3.1, 0.55],
        width_ratios=[1, 1, 0.07], hspace=0.12, wspace=0.12,
        left=0.19, right=0.90, top=0.90, bottom=0.07)

    axPT = fig.add_subplot(outer[0, 0])
    axPE = fig.add_subplot(outer[0, 1], sharey=axPT)
    axHT = fig.add_subplot(outer[1, 0])
    axHE = fig.add_subplot(outer[1, 1], sharey=axHT)
    axCB = fig.add_subplot(outer[1, 2])
    axGT = fig.add_subplot(outer[2, 0], sharex=axHT)
    axGE = fig.add_subplot(outer[2, 1], sharex=axHE)

    # title bar
    if args.title:
        fig.text(0.545, 0.955, args.title, ha="center", va="center",
                 fontsize=15, color="#6A2C91", fontweight="bold")
        fig.add_artist(plt.Line2D([0.19, 0.90], [0.93, 0.93],
                                  color="black", lw=2))

    # ---- profiles -------------------------------------------------------------
    xT = np.arange(-upT, dnT, binsz) + binsz / 2
    xE = np.arange(-upE, dnE, binsz) + binsz / 2
    colors = [PURPLE_DARK, GREEN_LIGHT]

    def band(reps_block):
        """mean profile and half-width of the band across replicates."""
        per_rep = np.nanmean(reps_block, axis=1)           # rep x bin
        mean = per_rep.mean(axis=0)
        n = per_rep.shape[0]
        if n < 2 or args.ci == "none":
            return mean, None
        sem = per_rep.std(axis=0, ddof=1) / np.sqrt(n)
        if args.ci == "sem":
            return mean, sem
        from scipy import stats as _st
        return mean, _st.t.ppf(0.975, n - 1) * sem

    for i, (lbl, a, b) in enumerate(groups):
        c = colors[i % len(colors)]
        for ax, x, R in ((axPT, xT, rT), (axPE, xE, rE)):
            mean, hw = band(R[:, a:b, :])
            if hw is not None:
                ax.fill_between(x, mean - hw, mean + hw, color=c, alpha=0.22,
                                linewidth=0)
            ax.plot(x, mean, color=c, lw=3.2, alpha=args.alpha,
                    label=nice_label(lbl).replace("\n", "") if ax is axPE else None)
    axPT.set_ylabel(args.ylabel or ("log2(target / IgG)" if args.signal == "log2ratio"
                                    else "IgG-subtracted\ncoverage"))
    if args.signal == "log2ratio":
        axPT.axhline(0, color="#888888", lw=1)
        if args.ymax:
            axPT.set_ylim(-args.ymax, args.ymax)
    elif args.ymax:
        axPT.set_ylim(0, args.ymax)
    else:
        axPT.set_ylim(0, None)
    leg = axPE.legend(frameon=False, loc="center right", handlelength=1.4)
    for lh in leg.get_lines():
        lh.set_alpha(1.0)
    if nrep > 1 and args.ci != "none":
        fig.text(0.90, 0.912, ("band: 95% CI" if args.ci == "95" else "band: +/-1 SEM")
                 + f" across n={nrep} replicates",
                 ha="right", va="top", fontsize=12, color="#555555")
    for ax in (axPT, axPE):
        ax.set_xlim(-upT if ax is axPT else -upE, dnT if ax is axPT else dnE)
        ax.tick_params(labelbottom=False, length=4, width=1.4)
        ax.spines[["top", "right"]].set_visible(False)
    axPE.tick_params(labelleft=False)

    # ---- heatmaps -------------------------------------------------------------
    imT = axHT.imshow(mT, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                      extent=[-upT, dnT, n_total, 0], interpolation="nearest")
    axHE.imshow(mE, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                extent=[-upE, dnE, n_total, 0], interpolation="nearest")
    for ax, x0 in ((axHT, 0), (axHE, 0)):
        ax.axvline(x0, color="black", ls=(0, (3, 2)), lw=1.3)
        ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
        for sp in ax.spines.values():
            sp.set_linewidth(1.8)
    # group separators + side bars
    y = 0
    for i, (lbl, a, b) in enumerate(groups):
        h = b - a
        if i > 0:
            for ax in (axHT, axHE):
                ax.axhline(y, color="black", lw=3)
        # axes coords run bottom-up; heatmap rows run top-down -> flip
        bar = Rectangle((-0.06, 1 - (y + h) / n_total), 0.025, h / n_total,
                        transform=axHT.transAxes, clip_on=False,
                        facecolor=colors[i % len(colors)], edgecolor="none")
        axHT.add_patch(bar)
        axHT.text(-0.10, 1 - (y + h / 2) / n_total, nice_label(lbl),
                  transform=axHT.transAxes, ha="right", va="center", fontsize=13)
        y += h
    cb = fig.colorbar(imT, cax=axCB)
    cb.ax.tick_params(labelsize=12, length=4)
    cb.outline.set_linewidth(1.5)

    # ---- gene cartoons ------------------------------------------------------
    for ax, up, dn, anchor in ((axGT, upT, dnT, "TSS"), (axGE, upE, dnE, "TES")):
        ax.set_xlim(-up, dn)
        ax.set_ylim(0, 1)
        ax.axis("off")
        if anchor == "TSS":
            ax.add_patch(Rectangle((0, 0.25), dn, 0.38, facecolor="#4D4D4D"))
            ax.plot([-up, 0], [0.44, 0.44], color="black", lw=1.5)
            ax.add_patch(FancyArrow(0, 0.75, dn * 0.18, 0, width=0.05,
                                    head_width=0.22, head_length=dn * 0.05,
                                    color="black", length_includes_head=True))
            ax.plot([0, 0], [0.44, 0.75], color="black", lw=1.5)
            ax.text(dn / 2, 0.44, "gene", ha="center", va="center",
                    color="white", fontsize=13)
            ax.text(-up, 1.0, f"{up//1000}kb", ha="left", va="bottom", fontsize=12)
            ax.text(dn, 1.0, f"{dn//1000}kb", ha="right", va="bottom", fontsize=12)
            ax.text(0, -0.12, "TSS", ha="center", va="top", fontsize=13)
        else:
            ax.add_patch(Rectangle((-up, 0.25), up, 0.38, facecolor="#4D4D4D"))
            ax.plot([0, dn], [0.44, 0.44], color="black", lw=1.5)
            ax.text(-up / 2, 0.44, "gene", ha="center", va="center",
                    color="white", fontsize=13)
            ax.text(-up, 1.0, f"{up//1000}kb", ha="left", va="bottom", fontsize=12)
            ax.text(dn, 1.0, f"{dn//1000}kb", ha="right", va="bottom", fontsize=12)
            ax.text(0, -0.12, "TES", ha="center", va="top", fontsize=13)

    for ext in ("pdf", "png"):
        fig.savefig(f"{args.out}.{ext}", dpi=200 if ext == "png" else None)
    print(f"wrote {args.out}.pdf / .png  (n={n_total}; replicates averaged={nrep}; "
          + ", ".join(f"{nice_label(l).replace(chr(10),' ')}={b-a}" for l, a, b in groups)
          + f"; vmax={vmax:.2f})")


if __name__ == "__main__":
    main()
