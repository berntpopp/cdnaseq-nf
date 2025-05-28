#!/usr/bin/env python3

"""
Custom junction filtering script for STAR two-pass alignment
"""

import sys
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Filter and merge STAR splice junctions')
    parser.add_argument('--input_files', nargs='+', help='Input SJ.out.tab files')
    parser.add_argument('--output', required=True, help='Output filtered junctions file')
    parser.add_argument('--min_unique_reads', type=int, default=3, help='Minimum unique reads supporting junction')
    parser.add_argument('--min_samples', type=int, default=1, help='Minimum samples supporting junction')
    return parser.parse_args()

def main():
    args = parse_args()
    
    junction_counts = defaultdict(lambda: {'samples': set(), 'total_reads': 0, 'max_reads': 0})
    
    for sj_file in args.input_files:
        sample_name = sj_file.split('/')[-1].replace('.SJ.out.tab', '')
        
        with open(sj_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    chrom, start, end, strand, motif, annotated, unique_reads, multimap_reads, max_overhang = fields
                    
                    unique_reads = int(unique_reads)
                    if unique_reads >= args.min_unique_reads:
                        junction_key = f"{chrom}:{start}-{end}:{strand}"
                        junction_counts[junction_key]['samples'].add(sample_name)
                        junction_counts[junction_key]['total_reads'] += unique_reads
                        junction_counts[junction_key]['max_reads'] = max(junction_counts[junction_key]['max_reads'], unique_reads)
    
    # Write filtered junctions
    with open(args.output, 'w') as f:
        for junction, data in junction_counts.items():
            if len(data['samples']) >= args.min_samples:
                chrom, coords, strand = junction.split(':')
                start, end = coords.split('-')
                # Write in STAR SJ format: chr start end strand motif annotated unique_reads multimap_reads max_overhang
                f.write(f"{chrom}\\t{start}\\t{end}\\t{strand}\\t1\\t0\\t{data['max_reads']}\\t0\\t50\\n")

if __name__ == '__main__':
    main()
