#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# make_gene_groups.py
# Authors : Katharina (Nina) Hayer, co-created with Claude Fable 5.1 (Anthropic)
# Version : 1.0.0
# Date    : 2026-10-06
# Purpose : From a GTF, write strand-aware BED6 files of protein-coding genes
#           split into intron-containing vs intronless, for computeMatrix.
# -----------------------------------------------------------------------------
import argparse
import gzip
import re
import sys
from collections import defaultdict

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

    n_short = 0
    kept = []                       # (chrom, start, end, name, strand, intronless)
    for gid, (chrom, s, e, strand, name) in genes.items():
        if e - s < args.min_length:
            n_short += 1
            continue
        txs = exons.get(gid, {})
        if not txs:
            continue
        # intronless = every transcript is single-exon (exon features, not CDS:
        # a contiguous CDS with a spliced 5'UTR is still a spliced gene)
        kept.append((chrom, s, e, name, strand, all(c == 1 for c in txs.values())))

    n_neighbour = drop_neighbours(kept, args) if args.neighbour_filter != "none" else 0
    if args.neighbour_filter != "none":
        kept = [g for g in kept if g is not None]

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
    print(f"dropped (< {args.min_length} bp): {n_short}", file=sys.stderr)
    if args.neighbour_filter != "none":
        print(f"dropped (--neighbour-filter {args.neighbour_filter}): "
              f"{n_neighbour['total']} "
              f"({n_neighbour['ic']} intron-containing, "
              f"{n_neighbour['il']} intronless)", file=sys.stderr)


if __name__ == "__main__":
    main()
