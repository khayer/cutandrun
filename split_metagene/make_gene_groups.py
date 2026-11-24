#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# make_gene_groups.py
# Authors : Katharina (Nina) Hayer, co-created with Claude Fable 5.1 (Anthropic)
# Version : 1.2.0
# Date    : 2026-11-24
# Purpose : From a GTF, write strand-aware BED6 files of protein-coding genes
#           split into intron-containing vs intronless, for computeMatrix.
#
# v1.1.0: the length filter is applied AFTER the intronless call, so the
#         dropped-for-length count can be reported per group, and a sweep over
#         candidate --min-length values is printed. Intronless protein-coding
#         genes are short as a class, so --min-length and the intronless group
#         pull against each other; the sweep makes that trade visible instead
#         of leaving it to be discovered from an n=22 curve.
#
# v1.2.0: --length-match. Intronless genes are short and intron-containing
#         genes are long, so the two groups differ by about an order of
#         magnitude in length before splicing status is considered at all, and
#         gene length drives body coverage on its own. Matching inside fixed
#         length bins makes the comparison about introns rather than length.
# -----------------------------------------------------------------------------
import argparse
import gzip
import random
import re
import sys
from collections import Counter, defaultdict

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_attrs(s):
    return dict(ATTR_RE.findall(s))


def drop_neighbours(kept, args):
    """Mark genes that share sequence with a neighbour as None, in place.

    BOTH members of an overlapping pair are dropped. Keeping one would mean
    choosing arbitrarily between two genes that are equally contaminated, and
    the kept one would still carry its neighbour's reads.
    """
    def window(g):
        chrom, s, e, _name, strand, _il = g
        if args.neighbour_filter == "any":
            return s, e
        return (s - args.tss_up, s + args.tss_down) if strand == "+" \
            else (e - args.tss_down, e + args.tss_up)

    by_chrom = defaultdict(list)
    for i, g in enumerate(kept):
        by_chrom[g[0]].append(i)

    doomed = set()
    for chrom, idxs in by_chrom.items():
        spans = sorted(((window(kept[i]), i) for i in idxs), key=lambda t: t[0][0])
        # sweep: compare each against those still open to its left
        active = []
        for (lo, hi), i in spans:
            active = [(a_hi, a_i) for a_hi, a_i in active if a_hi > lo]
            for _a_hi, a_i in active:
                doomed.add(i); doomed.add(a_i)
            active.append((hi, i))

    counts = {"total": len(doomed), "ic": 0, "il": 0}
    for i in doomed:
        counts["il" if kept[i][5] else "ic"] += 1
        kept[i] = None
    return counts


def median(values):
    v = sorted(values)
    return v[len(v) // 2] if v else 0


def length_match(kept, args):
    """Subsample intron-containing genes to the intronless length distribution.

    Matching is done inside fixed-width length bins, not on the mean, so the
    two groups share a length distribution bin for bin. A bin with too few
    intron-containing genes contributes what it has; the shortfall is reported
    rather than made up from a neighbouring bin, which would reintroduce the
    very skew the matching exists to remove.
    """
    rng = random.Random(args.seed)
    b = args.length_match_bin
    il = [g for g in kept if g[5]]
    ic_by_bin = defaultdict(list)
    for g in kept:
        if not g[5]:
            ic_by_bin[(g[2] - g[1]) // b].append(g)

    want = Counter((g[2] - g[1]) // b for g in il)
    picked = []
    deficit = 0
    for key in sorted(want):
        n = want[key] * args.length_match
        pool = sorted(ic_by_bin.get(key, []))      # deterministic input order
        if len(pool) <= n:
            picked.extend(pool)
            deficit += n - len(pool)
        else:
            picked.extend(rng.sample(pool, n))

    report = {
        "n_il": len(il),
        "n_ic": len(picked),
        "asked": len(il) * args.length_match,
        "deficit": deficit,
        "med_il": median([g[2] - g[1] for g in il]),
        "med_ic": median([g[2] - g[1] for g in picked]),
    }
    return il + picked, report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("gtf", help="GENCODE/Ensembl GTF (.gtf or .gtf.gz)")
    ap.add_argument("--out-prefix", default="genes", help="prefix for output BEDs")
    ap.add_argument("--min-length", type=int, default=6000,
                    help="min gene length (bp) so TSS and TES windows do not "
                         "overlap; default 6000 (window is 1+5 kb)")
    ap.add_argument("--biotype", default="protein_coding",
                    help="gene_type/gene_biotype to keep; 'all' for no filter")
    ap.add_argument("--neighbour-filter", choices=["none", "tss", "any"],
                    default="none",
                    help="drop genes that share sequence with a neighbour. "
                         "'tss': drop both genes of any pair whose TSS windows "
                         "overlap. 'any': drop both genes of any pair whose "
                         "spans overlap. Default 'none' keeps every gene, which "
                         "is the historical behaviour.")
    ap.add_argument("--tss-up", type=int, default=1000,
                    help="bp upstream of the TSS, for --neighbour-filter tss; "
                         "must match the matrix window to mean anything")
    ap.add_argument("--tss-down", type=int, default=5000,
                    help="bp downstream of the TSS, for --neighbour-filter tss")
    ap.add_argument("--length-match", type=int, default=0, metavar="N",
                    help="subsample intron-containing genes to N per intronless "
                         "gene within each length bin; 0 disables matching")
    ap.add_argument("--length-match-bin", type=int, default=500,
                    help="width in bp of the length bins used by --length-match")
    ap.add_argument("--seed", type=int, default=0,
                    help="seed for --length-match sampling, so reruns agree")
    ap.add_argument("--exclude-contigs", default="chrM,chrY",
                    help="comma list of contigs to drop")
    ap.add_argument("--keep-ad5", default="",
                    help="name of a viral contig to drop (host-only figure)")
    args = ap.parse_args()

    drop = set(c for c in args.exclude_contigs.split(",") if c)
    if args.keep_ad5:
        drop.add(args.keep_ad5)

    genes = {}                       # gene_id -> [chr, start, end, strand, name]
    exons = defaultdict(lambda: defaultdict(int))  # gene -> tx -> exon count

    with open_any(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, feat, start, end, strand = f[0], f[2], int(f[3]), int(f[4]), f[6]
            if chrom in drop:
                continue
            if feat not in ("gene", "exon"):
                continue
            a = parse_attrs(f[8])
            gid = a.get("gene_id")
            if not gid:
                continue
            if feat == "gene":
                bt = a.get("gene_type") or a.get("gene_biotype") or ""
                if args.biotype != "all" and bt != args.biotype:
                    continue
                genes[gid] = [chrom, start - 1, end, strand, a.get("gene_name", gid)]
            else:
                tx = a.get("transcript_id", "NA")
                exons[gid][tx] += 1

    # Classify first, filter on length second, so the length cost is known per
    # group. Classifying after the length filter hides which group it empties.
    classified = []                 # (chrom, start, end, name, strand, intronless)
    for gid, (chrom, s, e, strand, name) in genes.items():
        txs = exons.get(gid, {})
        if not txs:
            continue
        # intronless = every transcript is single-exon (exon features, not CDS:
        # a contiguous CDS with a spliced 5'UTR is still a spliced gene)
        classified.append((chrom, s, e, name, strand,
                           all(c == 1 for c in txs.values())))

    kept = [g for g in classified if g[2] - g[1] >= args.min_length]
    short = [g for g in classified if g[2] - g[1] < args.min_length]
    n_short = len(short)
    n_short_ic = sum(1 for g in short if not g[5])
    n_short_il = sum(1 for g in short if g[5])

    n_neighbour = drop_neighbours(kept, args) if args.neighbour_filter != "none" else 0
    if args.neighbour_filter != "none":
        kept = [g for g in kept if g is not None]

    match_report = None
    if args.length_match > 0:
        kept, match_report = length_match(kept, args)

    n_ic = n_il = 0
    with open(f"{args.out_prefix}.intron_containing.bed", "w") as ic, \
         open(f"{args.out_prefix}.intronless.bed", "w") as il:
        for chrom, s, e, name, strand, intronless in sorted(kept):
            row = f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}\n"
            if intronless:
                il.write(row); n_il += 1
            else:
                ic.write(row); n_ic += 1

    print(f"intron-containing: {n_ic}", file=sys.stderr)
    print(f"intronless       : {n_il}", file=sys.stderr)
    print(f"dropped (< {args.min_length} bp): {n_short} "
          f"({n_short_ic} intron-containing, {n_short_il} intronless)",
          file=sys.stderr)
    if args.neighbour_filter != "none":
        print(f"dropped (--neighbour-filter {args.neighbour_filter}): "
              f"{n_neighbour['total']} "
              f"({n_neighbour['ic']} intron-containing, "
              f"{n_neighbour['il']} intronless)", file=sys.stderr)

    # --min-length sweep, before the neighbour filter. A TSS window of
    # -tss_up/+tss_down and a TES window only separate when the gene is at
    # least tss_up + tss_down long, which is where the 6000 bp default comes
    # from; these counts say what that costs each group.
    if match_report is not None:
        r = match_report
        print(f"--length-match {args.length_match}: {r['n_ic']} intron-containing "
              f"kept for {r['n_il']} intronless "
              f"(asked {r['asked']}, short by {r['deficit']})", file=sys.stderr)
        print(f"  median length: {r['med_il']} intronless, "
              f"{r['med_ic']} intron-containing", file=sys.stderr)
        if r["deficit"] > 0:
            print("  NOTE: some length bins held fewer intron-containing genes "
                  "than asked; those bins are under-filled rather than topped "
                  "up from adjacent bins.", file=sys.stderr)

    floor = args.tss_up + args.tss_down
    sweep = sorted({1000, 2000, 3000, 5000, floor, args.min_length})
    print("min-length sweep (before --neighbour-filter):", file=sys.stderr)
    for cut in sweep:
        ic = sum(1 for g in classified if g[2] - g[1] >= cut and not g[5])
        il = sum(1 for g in classified if g[2] - g[1] >= cut and g[5])
        mark = "  <- current" if cut == args.min_length else ""
        print(f"  >= {cut:>6} bp: {ic:>6} intron-containing, "
              f"{il:>5} intronless{mark}", file=sys.stderr)


if __name__ == "__main__":
    main()
