#!/usr/bin/env python3

"""Select top consensus peaks per group and annotate nearest gene for plotting."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class Peak:
    chrom: str
    start: int
    end: int
    score: float
    rank: int = 0


@dataclass
class Gene:
    start: int
    end: int
    name: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare top consensus peak regions for pyGenomeTracks.")
    parser.add_argument("--consensus-bed", required=True, help="Consensus peak BED file.")
    parser.add_argument("--gtf", required=True, help="GTF file used for nearest gene labeling.")
    parser.add_argument("--group-id", required=True, help="Group identifier used in titles and filenames.")
    parser.add_argument("--top-n", type=int, default=10, help="Number of top peaks to select.")
    parser.add_argument("--flank", type=int, default=1000, help="Flank size in bp around peak midpoint.")
    parser.add_argument("--out-regions-tsv", required=True, help="Output TSV with selected plotting regions.")
    parser.add_argument("--out-top-bed", required=True, help="Output BED with selected top peaks.")
    return parser.parse_args()


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("_") or "na"


def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_field.split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, val = part.split(" ", 1)
        attrs[key.strip()] = val.strip().strip('"')
    return attrs


def load_genes(gtf_path: Path) -> Dict[str, List[Gene]]:
    genes_by_chrom: Dict[str, List[Gene]] = {}
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attrs = fields
            if feature not in {"gene", "transcript", "exon"}:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            parsed = parse_gtf_attributes(attrs)
            gene_name = (
                parsed.get("gene_name")
                or parsed.get("gene_id")
                or parsed.get("Name")
                or parsed.get("transcript_id")
                or "NA"
            )
            genes_by_chrom.setdefault(chrom, []).append(Gene(start=start_i, end=end_i, name=gene_name))

    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda g: (g.start, g.end, g.name))
    return genes_by_chrom


def nearest_gene(chrom: str, pos: int, genes_by_chrom: Dict[str, List[Gene]]) -> Tuple[str, int]:
    genes = genes_by_chrom.get(chrom, [])
    if not genes:
        return "NA", -1

    best_name = "NA"
    best_dist = None
    for gene in genes:
        if gene.start <= pos <= gene.end:
            return gene.name, 0
        dist = min(abs(pos - gene.start), abs(pos - gene.end))
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_name = gene.name
    return best_name, int(best_dist if best_dist is not None else -1)


def load_consensus_peaks(consensus_path: Path) -> List[Peak]:
    peaks: List[Peak] = []
    with consensus_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            score = 0.0
            if len(fields) >= 10:
                try:
                    score = float(fields[9])
                except ValueError:
                    score = 0.0
            elif len(fields) >= 5:
                try:
                    score = float(fields[4])
                except ValueError:
                    score = 0.0
            peaks.append(Peak(chrom=fields[0], start=start, end=end, score=score))

    peaks.sort(key=lambda p: (-p.score, -(p.end - p.start), p.chrom, p.start, p.end))
    return peaks


def write_outputs(
    peaks: List[Peak],
    genes_by_chrom: Dict[str, List[Gene]],
    group_id: str,
    flank: int,
    out_regions_tsv: Path,
    out_top_bed: Path,
) -> int:
    selected = peaks
    with out_regions_tsv.open("w", encoding="utf-8") as regions, out_top_bed.open("w", encoding="utf-8") as bed:
        regions.write(
            "rank\tchrom\tpeak_start\tpeak_end\twindow_start\twindow_end\tnearest_gene\tdistance_bp\ttitle\tfile_stub\n"
        )

        for idx, peak in enumerate(selected, start=1):
            peak.rank = idx
            center = (peak.start + peak.end) // 2
            window_start = max(0, center - flank)
            window_end = center + flank
            gene_name, distance = nearest_gene(peak.chrom, center, genes_by_chrom)

            title = (
                f"{group_id} | rank {idx} | {peak.chrom}:{peak.start}-{peak.end} | "
                f"nearest gene: {gene_name} ({distance} bp)"
            )
            stub = safe_name(f"{group_id}_top{idx:02d}_{gene_name}")

            regions.write(
                "\t".join(
                    [
                        str(idx),
                        peak.chrom,
                        str(peak.start),
                        str(peak.end),
                        str(window_start),
                        str(window_end),
                        gene_name,
                        str(distance),
                        title,
                        stub,
                    ]
                )
                + "\n"
            )
            bed.write(
                f"{peak.chrom}\t{peak.start}\t{peak.end}\t{group_id}_top{idx:02d}\t{int(peak.score)}\t.\n"
            )

    return len(selected)


def main() -> None:
    args = parse_args()

    consensus_path = Path(args.consensus_bed)
    gtf_path = Path(args.gtf)
    out_regions_tsv = Path(args.out_regions_tsv)
    out_top_bed = Path(args.out_top_bed)

    peaks = load_consensus_peaks(consensus_path)
    selected = peaks[: max(0, args.top_n)]
    genes_by_chrom = load_genes(gtf_path)

    count = write_outputs(
        peaks=selected,
        genes_by_chrom=genes_by_chrom,
        group_id=args.group_id,
        flank=max(0, args.flank),
        out_regions_tsv=out_regions_tsv,
        out_top_bed=out_top_bed,
    )

    print(f"Prepared {count} regions for group {args.group_id}")


if __name__ == "__main__":
    main()
