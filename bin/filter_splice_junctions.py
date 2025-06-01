#!/usr/bin/env python3

"""
Custom junction filtering script for STAR two-pass alignment
Implements the original awk logic:
1. Keep non-chrM, annotated junctions (col 6 is 1)
2. Keep non-chrM, unannotated (col 6 is 0) junctions with >5 unique reads
"""

import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter and merge STAR splice junctions"
    )
    parser.add_argument("--input_files", nargs="+", help="Input SJ.out.tab files")
    parser.add_argument(
        "--output", required=True, help="Output filtered junctions file"
    )
    parser.add_argument(
        "--min_reads_unannotated",
        type=int,
        default=5,
        help="Min unique reads for unannotated junctions",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Store junction information: key -> {motif, annotation, total_reads,
    # max_overhang}
    junctions = defaultdict(
        lambda: {
            "motif": "1",
            "annotation": "0",
            "total_reads": 0,
            "max_overhang": "50",
        }
    )

    for sj_file in args.input_files:
        with open(sj_file, "r") as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split("\t")
                    if len(fields) >= 9:
                        (
                            chrom,
                            start,
                            end,
                            strand,
                            motif,
                            annotated,
                            unique_reads,
                            multimap_reads,
                            max_overhang,
                        ) = fields[:9]

                        # Skip chrM junctions
                        if chrom == "chrM":
                            continue

                        unique_reads = int(unique_reads)
                        annotated = int(annotated)

                        junction_key = f"{chrom}:{start}-{end}:{strand}"

                        # Aggregate data for this junction
                        junctions[junction_key]["motif"] = motif
                        junctions[junction_key]["annotation"] = str(annotated)
                        junctions[junction_key]["total_reads"] += unique_reads
                        junctions[junction_key]["max_overhang"] = max_overhang

    # Write filtered junctions following original awk logic
    with open(args.output, "w") as f:
        for junction_key, data in junctions.items():
            chrom, coords, strand = junction_key.split(":")
            start, end = coords.split("-")

            # Original awk logic:
            # 1. Keep annotated junctions (annotation == 1)
            # 2. Keep unannotated junctions (annotation == 0) with
            #    >min_reads_unannotated
            annotation = int(data["annotation"])
            total_reads = data["total_reads"]

            keep_junction = annotation == 1 or (
                annotation == 0 and total_reads > args.min_reads_unannotated
            )

            if keep_junction:
                # Write in STAR SJ format
                f.write(
                    f"{chrom}\t{start}\t{end}\t{strand}\t"
                    f"{data['motif']}\t{data['annotation']}\t"
                    f"{total_reads}\t0\t{data['max_overhang']}\n"
                )


if __name__ == "__main__":
    main()
