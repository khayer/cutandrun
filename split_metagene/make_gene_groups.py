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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("gtf", help="GENCODE/Ensembl GTF (.gtf or .gtf.gz)")
    ap.add_argument("--out-prefix", default="genes", help="prefix for output BEDs")
    ap.add_argument("--min-length", type=int, default=6000,
                    help="min gene length (bp) so TSS and TES windows do not "
                         "overlap; default 6000 (window is 1+5 kb)")
    ap.add_argument("--biotype", default="protein_coding",
                    help="gene_type/gene_biotype to keep; 'all' for no filter")
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

    n_ic = n_il = n_short = 0
    with open(f"{args.out_prefix}.intron_containing.bed", "w") as ic, \
         open(f"{args.out_prefix}.intronless.bed", "w") as il:
        for gid, (chrom, s, e, strand, name) in genes.items():
            if e - s < args.min_length:
                n_short += 1
                continue
            txs = exons.get(gid, {})
            if not txs:
                continue
            # intronless = every transcript is single-exon
            intronless = all(c == 1 for c in txs.values())
            row = f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}\n"
            if intronless:
                il.write(row); n_il += 1
            else:
                ic.write(row); n_ic += 1

    print(f"intron-containing: {n_ic}", file=sys.stderr)
    print(f"intronless       : {n_il}", file=sys.stderr)
    print(f"dropped (< {args.min_length} bp): {n_short}", file=sys.stderr)


if __name__ == "__main__":
    main()
