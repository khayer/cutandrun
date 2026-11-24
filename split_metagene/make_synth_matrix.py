#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_synth_matrix.py
# Authors: Katharina E. Hayer, co-created with Claude
# Version: 1.1.0
# Date:    2026-11-23
# ---------------------------------------------------------------------------
# Writes a pair of deepTools computeMatrix files (TSS- and TES-anchored) filled
# with SYNTHETIC signal, so plot_split_metagene.py can be demonstrated without
# any real data. Nothing here comes from an experiment.
#
# The output is the real computeMatrix format, not a mock: a gzipped file whose
# first line is '@' followed by the JSON header, then one tab-separated row per
# region (chrom, start, end, name, score, strand, then the values, with every
# sample's bins laid side by side). The header keys plot_split_metagene.py
# actually reads are sample_boundaries, group_boundaries, group_labels,
# upstream, downstream and 'bin size'.
#
# Group sizes deliberately mirror the real imbalance (many intron-containing,
# few intronless), because a demo that showed two equal bands would misrepresent
# what the figure looks like on real input.
# ---------------------------------------------------------------------------

import argparse, gzip, json
import numpy as np


def build(path, up, down, binsize, samples, n_ic, n_il, seed, anchored_peak):
    rng = np.random.default_rng(seed)
    nbins = (up + down) // binsize
    x = np.arange(nbins)
    anchor = up // binsize

    def profile_ic():
        y = np.zeros(nbins)
        if anchored_peak:                      # TSS: peak just downstream
            y += 0.33                                   # upstream baseline
            y += 2.4 * np.exp(-0.5 * ((x - (anchor + 3)) / 2.4) ** 2)
            body = x > anchor + 3
            # concave rise through the gene body, not a straight ramp
            t = np.linspace(0, 1, body.sum())
            y[body] += 0.5 + 1.7 * t ** 0.75
        else:                                  # TES: plateau then a drop
            y += 2.5
            tail = x > anchor
            y[tail] = np.linspace(2.4, 0.9, tail.sum())
        return y

    def profile_il():
        return np.full(nbins, 0.42)

    rows = []
    for grp, n, prof in (("ic", n_ic, profile_ic), ("il", n_il, profile_il)):
        for i in range(n):
            gene_scale = rng.lognormal(0.0, 0.55)       # gene-to-gene spread
            base = prof() * gene_scale
            per_sample = []
            for s in range(samples):
                rep = base * rng.normal(1.0, 0.07)      # replicate effect
                rep = rep + rng.normal(0, 0.18, nbins)  # bin noise
                per_sample.append(np.clip(rep, 0, None))
            rows.append((grp, i, np.concatenate(per_sample)))

    hdr = {
        "upstream": [up] * samples,
        "downstream": [down] * samples,
        "body": [0] * samples,
        "bin size": [binsize] * samples,
        "ref point": ["TSS" if anchored_peak else "TES"] * samples,
        "verbose": False,
        "bin avg type": "mean",
        "missing data as zero": True,
        "min threshold": None,
        "max threshold": None,
        "scale": 1,
        "skip zeros": False,
        "nan after end": False,
        "proc number": 1,
        "sort regions": "keep",
        "sort using": "mean",
        "unscaled 5 prime": [0] * samples,
        "unscaled 3 prime": [0] * samples,
        "group_labels": ["genes.intron_containing.bed", "genes.intronless.bed"],
        "group_boundaries": [0, n_ic, n_ic + n_il],
        "sample_labels": [f"Target_A_R{s+1}" for s in range(samples)],
        "sample_boundaries": [nbins * s for s in range(samples + 1)],
    }

    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(hdr) + "\n")
        for grp, i, vals in rows:
            start = 100000 + i * 20000
            name = f"{'GENE' if grp == 'ic' else 'SOLO'}{i+1:05d}"
            strand = "+" if i % 2 == 0 else "-"
            fh.write("\t".join(["chrS", str(start), str(start + 12000),
                                name, "0", strand]
                               + [f"{v:.4f}" for v in vals]) + "\n")
    print(f"  wrote {path}: {len(rows)} regions x {nbins} bins x {samples} samples")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="synthetic_demo")
    ap.add_argument("--samples", type=int, default=3)
    ap.add_argument("--n-intron-containing", type=int, default=900)
    ap.add_argument("--n-intronless", type=int, default=40)
    ap.add_argument("--bin-size", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20261123)
    a = ap.parse_args()

    build(f"{a.prefix}.TSS.mat.gz", 1000, 5000, a.bin_size, a.samples,
          a.n_intron_containing, a.n_intronless, a.seed, True)
    build(f"{a.prefix}.TES.mat.gz", 5000, 1000, a.bin_size, a.samples,
          a.n_intron_containing, a.n_intronless, a.seed + 1, False)
    print("  these matrices contain synthetic signal only")


if __name__ == "__main__":
    main()
