#!/usr/bin/env python3

"""Select top consensus peaks per group and annotate nearest gene for plotting."""

from __future__ import annotations

import argparse
import csv
import re
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class Peak:
    chrom: str
    start: int
    end: int
    score: float
    input_rank: int = 0
    rank: int = 0


@dataclass
class Gene:
    start: int
    end: int
    name: str


@dataclass
class IntervalIndex:
    starts_by_chrom: Dict[str, List[int]]
    intervals_by_chrom: Dict[str, List[Tuple[int, int]]]


@dataclass
class SelectedPeak:
    peak: Peak
    feature: str
    feature_rank: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare top consensus peak regions for pyGenomeTracks.")
    parser.add_argument("--consensus-bed", required=True, help="Consensus peak BED file.")
    parser.add_argument("--gtf", required=True, help="GTF file used for nearest gene labeling.")
    parser.add_argument("--chipseeker-annot", required=False, help="ChIPseeker annotation table for the same consensus BED.")
    parser.add_argument("--homer-annot", required=False, help="HOMER annotatePeaks table for the same consensus BED.")
    parser.add_argument("--group-id", required=True, help="Group identifier used in titles and filenames.")
    parser.add_argument("--top-n", type=int, default=10, help="Number of top peaks to select.")
    parser.add_argument("--flank", type=int, default=1000, help="Flank size in bp around peak midpoint.")
    parser.add_argument(
        "--window-bp",
        type=int,
        default=0,
        help="Exact window size in bp centered on peak midpoint. If > 0, overrides --flank.",
    )
    parser.add_argument(
        "--rank-mode",
        choices=["global", "by_feature"],
        default="global",
        help="Selection mode: global top-N peaks or top-N per genomic feature class.",
    )
    parser.add_argument(
        "--feature-types",
        default="promoter,TES,UTR,exon,intron,intergenic",
        help="Comma-separated feature classes when --rank-mode by_feature is enabled.",
    )
    parser.add_argument(
        "--anchor-window",
        type=int,
        default=3000,
        help="Window in bp around TSS/TES used to classify promoter and TES peaks.",
    )
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


def add_interval(store: Dict[str, List[Tuple[int, int]]], chrom: str, start: int, end: int) -> None:
    store.setdefault(chrom, []).append((start, end))


def build_interval_index(intervals_by_chrom: Dict[str, List[Tuple[int, int]]]) -> IntervalIndex:
    starts: Dict[str, List[int]] = {}
    norm_intervals: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, intervals in intervals_by_chrom.items():
        cleaned = [(s, e) for s, e in intervals if e >= s]
        cleaned.sort(key=lambda x: (x[0], x[1]))
        norm_intervals[chrom] = cleaned
        starts[chrom] = [s for s, _ in cleaned]
    return IntervalIndex(starts_by_chrom=starts, intervals_by_chrom=norm_intervals)


def point_overlaps(chrom: str, pos: int, index: IntervalIndex) -> bool:
    starts = index.starts_by_chrom.get(chrom)
    if not starts:
        return False
    intervals = index.intervals_by_chrom[chrom]
    i = bisect_right(starts, pos) - 1
    if i < 0:
        return False
    s, e = intervals[i]
    return s <= pos <= e


def range_overlaps(chrom: str, start: int, end: int, index: IntervalIndex) -> bool:
    starts = index.starts_by_chrom.get(chrom)
    if not starts:
        return False
    intervals = index.intervals_by_chrom[chrom]
    i = bisect_right(starts, end) - 1
    while i >= 0:
        s, e = intervals[i]
        if e < start:
            return False
        if s <= end and e >= start:
            return True
        i -= 1
    return False


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


def build_feature_indexes(gtf_path: Path, anchor_window: int) -> Dict[str, IntervalIndex]:
    transcript_intervals: Dict[str, List[Tuple[int, int]]] = {}
    exon_intervals: Dict[str, List[Tuple[int, int]]] = {}
    utr_intervals: Dict[str, List[Tuple[int, int]]] = {}
    promoter_intervals: Dict[str, List[Tuple[int, int]]] = {}
    tes_intervals: Dict[str, List[Tuple[int, int]]] = {}

    with gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, _ = fields
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            feature_l = feature.lower()
            if feature_l in {"transcript", "gene"}:
                add_interval(transcript_intervals, chrom, start_i, end_i)
                if strand == "-":
                    tss = end_i
                    tes = start_i
                else:
                    tss = start_i
                    tes = end_i
                add_interval(promoter_intervals, chrom, max(1, tss - anchor_window), tss + anchor_window)
                add_interval(tes_intervals, chrom, max(1, tes - anchor_window), tes + anchor_window)
            elif feature_l == "exon":
                add_interval(exon_intervals, chrom, start_i, end_i)
            elif feature_l in {"utr", "three_prime_utr", "five_prime_utr", "three_prime_utr", "five_prime_utr"}:
                add_interval(utr_intervals, chrom, start_i, end_i)

    return {
        "promoter": build_interval_index(promoter_intervals),
        "TES": build_interval_index(tes_intervals),
        "UTR": build_interval_index(utr_intervals),
        "exon": build_interval_index(exon_intervals),
        "transcript": build_interval_index(transcript_intervals),
    }


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
        for input_rank, line in enumerate(handle, start=1):
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
            peaks.append(Peak(chrom=fields[0], start=start, end=end, score=score, input_rank=input_rank))

    peaks.sort(key=lambda p: (-p.score, -(p.end - p.start), p.chrom, p.start, p.end))
    return peaks


def load_homer_annotations(homer_annot_path: Path) -> Dict[object, str]:
    annotations: Dict[object, str] = {}
    if not homer_annot_path.exists():
        return annotations

    with homer_annot_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return annotations

        annotation_col = None
        for candidate in ["Annotation", "Detailed Annotation"]:
            for column in reader.fieldnames:
                if column and column.strip().lower() == candidate.lower():
                    annotation_col = column
                    break
            if annotation_col:
                break
        if annotation_col is None:
            return annotations

        chr_col = None
        start_col = None
        end_col = None
        for column in reader.fieldnames:
            lowered = column.strip().lower() if column else ""
            if chr_col is None and lowered in {"chr", "chrom", "chromosome"}:
                chr_col = column
            elif start_col is None and lowered == "start":
                start_col = column
            elif end_col is None and lowered == "end":
                end_col = column

        if not chr_col or not start_col or not end_col:
            return annotations

        for row_index, row in enumerate(reader, start=1):
            try:
                chrom = str(row[chr_col]).strip()
                start = int(str(row[start_col]).strip())
                end = int(str(row[end_col]).strip())
            except Exception:
                continue
            label = canonical_homer_annotation(str(row.get(annotation_col, "")))
            annotations[(chrom, start, end)] = label

            peak_id = str(row.get("PeakID", "")).strip()
            match = re.match(r"^peak_(\d+)$", peak_id, flags=re.IGNORECASE)
            if match:
                annotations[int(match.group(1))] = label
            else:
                annotations[row_index] = label

    return annotations


def load_chipseeker_annotations(chipseeker_annot_path: Path) -> Dict[Tuple[str, int, int], str]:
    annotations: Dict[Tuple[str, int, int], str] = {}
    if not chipseeker_annot_path.exists():
        return annotations

    with chipseeker_annot_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return annotations

        annotation_col = None
        for candidate in ["annotation", "Annotation"]:
            for column in reader.fieldnames:
                if column and column.strip().lower() == candidate.lower():
                    annotation_col = column
                    break
            if annotation_col:
                break
        if annotation_col is None:
            return annotations

        chr_col = None
        start_col = None
        end_col = None
        for column in reader.fieldnames:
            lowered = column.strip().lower() if column else ""
            if chr_col is None and lowered in {"seqnames", "seqname", "chr", "chrom", "chromosome"}:
                chr_col = column
            elif start_col is None and lowered == "start":
                start_col = column
            elif end_col is None and lowered == "end":
                end_col = column

        if not chr_col or not start_col or not end_col:
            return annotations

        for row in reader:
            try:
                chrom = str(row[chr_col]).strip()
                start = int(str(row[start_col]).strip())
                end = int(str(row[end_col]).strip())
            except Exception:
                continue
            annotations[(chrom, start, end)] = canonical_homer_annotation(str(row.get(annotation_col, "")))

    return annotations


def canonical_homer_annotation(annotation: str) -> str:
    if not annotation:
        return "Unknown"

    ann = str(annotation).strip()
    ann_lower = ann.lower()

    if ann_lower.startswith("promoter") or "tss" in ann_lower:
        return "promoter"
    if ann_lower.startswith("tts"):
        return "TES"
    if "5' utr" in ann_lower or "5 utr" in ann_lower or "3' utr" in ann_lower or "3 utr" in ann_lower:
        return "UTR"
    if ann_lower.startswith("exon"):
        return "exon"
    if ann_lower.startswith("intron"):
        return "intron"
    if "intergenic" in ann_lower:
        return "intergenic"
    return ann


def classify_feature(
    peak: Peak,
    feature_indexes: Dict[str, IntervalIndex],
    chipseeker_annotations: Dict[Tuple[str, int, int], str],
    homer_annotations: Dict[object, str],
) -> str:
    chrom = peak.chrom
    start = peak.start
    end = peak.end

    chipseeker_label = chipseeker_annotations.get((chrom, start, end))
    if chipseeker_label:
        return chipseeker_label

    homer_label = homer_annotations.get((chrom, start, end))
    if homer_label:
        return homer_label

    homer_label = homer_annotations.get(peak.input_rank)
    if homer_label:
        return homer_label

    # Precedence is intentional: if a peak overlaps multiple feature classes,
    # prefer promoter/TES/UTR labels before exon/intron.
    if range_overlaps(chrom, start, end, feature_indexes["promoter"]):
        return "promoter"
    if range_overlaps(chrom, start, end, feature_indexes["TES"]):
        return "TES"
    if range_overlaps(chrom, start, end, feature_indexes["UTR"]):
        return "UTR"
    if range_overlaps(chrom, start, end, feature_indexes["exon"]):
        return "exon"
    if range_overlaps(chrom, start, end, feature_indexes["transcript"]):
        return "intron"
    return "intergenic"


def parse_feature_types(value: str) -> List[str]:
    allowed = ["promoter", "TES", "UTR", "exon", "intron", "intergenic"]
    wanted = [x.strip() for x in value.split(",") if x.strip()]
    if not wanted:
        return allowed
    normalized: List[str] = []
    for item in wanted:
        item_l = item.lower()
        if item_l == "tes":
            normalized.append("TES")
        elif item_l == "utr":
            normalized.append("UTR")
        elif item_l in {"promoter", "exon", "intron", "intergenic"}:
            normalized.append(item_l)
    return [x for x in allowed if x in normalized]


def select_peaks(
    peaks: List[Peak],
    top_n: int,
    rank_mode: str,
    feature_types: List[str],
    feature_indexes: Dict[str, IntervalIndex],
    chipseeker_annotations: Dict[Tuple[str, int, int], str],
    homer_annotations: Dict[object, str],
) -> List[SelectedPeak]:
    if rank_mode == "global":
        selected: List[SelectedPeak] = []
        for i, peak in enumerate(peaks[: max(0, top_n)], start=1):
            selected.append(
                SelectedPeak(
                    peak=peak,
                    feature=classify_feature(peak, feature_indexes, chipseeker_annotations, homer_annotations),
                    feature_rank=i,
                )
            )
        return selected

    limits = {feat: max(0, top_n) for feat in feature_types}
    buckets: Dict[str, List[Peak]] = {feat: [] for feat in feature_types}
    complete = set()
    for peak in peaks:
        feature = classify_feature(peak, feature_indexes, chipseeker_annotations, homer_annotations)
        if feature not in buckets:
            continue
        if len(buckets[feature]) >= limits[feature]:
            continue
        buckets[feature].append(peak)
        if len(buckets[feature]) >= limits[feature]:
            complete.add(feature)
        if len(complete) == len(limits):
            break

    selected = []
    for feature in feature_types:
        for idx, peak in enumerate(buckets[feature], start=1):
            selected.append(SelectedPeak(peak=peak, feature=feature, feature_rank=idx))
    return selected


def write_outputs(
    selected: List[SelectedPeak],
    genes_by_chrom: Dict[str, List[Gene]],
    group_id: str,
    flank: int,
    window_bp: int,
    out_regions_tsv: Path,
    out_top_bed: Path,
) -> int:
    with out_regions_tsv.open("w", encoding="utf-8") as regions, out_top_bed.open("w", encoding="utf-8") as bed:
        regions.write(
            "rank\tchrom\tpeak_start\tpeak_end\twindow_start\twindow_end\tnearest_gene\tdistance_bp\ttitle\tfile_stub\n"
        )

        for idx, item in enumerate(selected, start=1):
            peak = item.peak
            peak.rank = idx
            center = (peak.start + peak.end) // 2
            if window_bp > 0:
                half = window_bp // 2
                window_start = max(0, center - half)
                window_end = window_start + window_bp
            else:
                window_start = max(0, center - flank)
                window_end = center + flank
            gene_name, distance = nearest_gene(peak.chrom, center, genes_by_chrom)

            title = (
                f"{group_id} | {item.feature} rank {item.feature_rank} | {peak.chrom}:{peak.start}-{peak.end} | "
                f"nearest gene: {gene_name} ({distance} bp)"
            )
            stub = safe_name(f"{group_id}_{item.feature}_top{item.feature_rank:02d}_{gene_name}")

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
                f"{peak.chrom}\t{peak.start}\t{peak.end}\t{group_id}_{item.feature}_top{item.feature_rank:02d}\t{int(peak.score)}\t.\n"
            )

    return len(selected)


def main() -> None:
    args = parse_args()

    consensus_path = Path(args.consensus_bed)
    gtf_path = Path(args.gtf)
    out_regions_tsv = Path(args.out_regions_tsv)
    out_top_bed = Path(args.out_top_bed)

    peaks = load_consensus_peaks(consensus_path)
    genes_by_chrom = load_genes(gtf_path)
    feature_indexes = build_feature_indexes(gtf_path, max(0, args.anchor_window))
    feature_types = parse_feature_types(args.feature_types)
    chipseeker_annotations = load_chipseeker_annotations(Path(args.chipseeker_annot)) if args.chipseeker_annot else {}
    homer_annotations = load_homer_annotations(Path(args.homer_annot)) if args.homer_annot else {}

    selected = select_peaks(
        peaks=peaks,
        top_n=max(0, args.top_n),
        rank_mode=args.rank_mode,
        feature_types=feature_types,
        feature_indexes=feature_indexes,
        chipseeker_annotations=chipseeker_annotations,
        homer_annotations=homer_annotations,
    )

    count = write_outputs(
        selected=selected,
        genes_by_chrom=genes_by_chrom,
        group_id=args.group_id,
        flank=max(0, args.flank),
        window_bp=max(0, args.window_bp),
        out_regions_tsv=out_regions_tsv,
        out_top_bed=out_top_bed,
    )

    print(f"Prepared {count} regions for group {args.group_id}")


if __name__ == "__main__":
    main()
