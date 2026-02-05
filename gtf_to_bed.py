#!/usr/bin/env python3
"""
Convert GTF file to BED12 format.
Groups exons by transcript/parent and creates proper BED12 entries with blocks.
"""

import re
import sys
from collections import defaultdict

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string and extract gene name, ID, and Parent."""
    result = {}
    
    # Look for gene "GENE_NAME" pattern
    gene_match = re.search(r'gene "([^"]+)"', attr_string)
    if gene_match:
        result['gene'] = gene_match.group(1)
    
    # Look for ID "..." pattern
    id_match = re.search(r'ID "([^"]+)"', attr_string)
    if id_match:
        result['id'] = id_match.group(1)
    
    # Look for Parent "..." pattern
    parent_match = re.search(r'Parent "([^"]+)"', attr_string)
    if parent_match:
        result['parent'] = parent_match.group(1)
    
    return result

def gtf_to_bed(gtf_file, bed_file):
    """
    Convert GTF to BED12 format with ID and Parent in columns 13-14.
    
    Parameters:
    - gtf_file: Input GTF file path
    - bed_file: Output BED file path
    """
    # Dictionary to group exons by parent (transcript/RNA)
    transcripts = defaultdict(list)
    
    with open(gtf_file, 'r') as infile:
        for line in infile:
            # Skip comments and empty lines
            if line.startswith('#') or line.strip() == '':
                continue
            
            fields = line.strip().split('\t')
            
            # We need at least 9 fields for a valid GTF
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])  # GTF is 1-based
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            attributes = fields[8]
            
            # Process both exon and CDS features (CDS when no exons exist)
            if feature_type not in ['exon', 'CDS']:
                continue
            
            # Extract attributes
            attrs = parse_gtf_attributes(attributes)
            if 'gene' not in attrs or 'parent' not in attrs:
                continue
            
            gene_name = attrs['gene']
            exon_id = attrs.get('id', '.')
            parent_id = attrs['parent']
            
            # Store exon/CDS information grouped by parent
            # Use feature_type to track what we're storing
            transcripts[parent_id].append({
                'chrom': chrom,
                'start': start - 1,  # Convert to 0-based
                'end': end,
                'gene': gene_name,
                'score': score if score != '.' else '0',
                'strand': strand,
                'exon_id': exon_id,
                'feature_type': feature_type
            })
    
    # Write BED12 format
    with open(bed_file, 'w') as outfile:
        for parent_id, exons in sorted(transcripts.items()):
            # Sort exons by start position
            exons.sort(key=lambda x: x['start'])
            
            # Get transcript bounds
            chrom = exons[0]['chrom']
            tx_start = min(e['start'] for e in exons)
            tx_end = max(e['end'] for e in exons)
            gene_name = exons[0]['gene']
            score = exons[0]['score']
            strand = exons[0]['strand']
            
            # Calculate block information
            block_count = len(exons)
            block_sizes = ','.join(str(e['end'] - e['start']) for e in exons) + ','
            block_starts = ','.join(str(e['start'] - tx_start) for e in exons) + ','
            
            # Collect all exon IDs for this transcript
            exon_ids = ','.join(e['exon_id'] for e in exons)
            
            # BED12 format: chrom, chromStart, chromEnd, name, score, strand, 
            #               thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
            # Plus columns 13-14: exon IDs, parent ID
            bed_line = (f'{chrom}\t{tx_start}\t{tx_end}\t{gene_name}\t{score}\t{strand}\t'
                       f'{tx_start}\t{tx_end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\t'
                       f'{exon_ids}\t{parent_id}\n')
            outfile.write(bed_line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python gtf_to_bed.py <input.gtf> <output.bed>")
        sys.exit(1)
    
    input_gtf = sys.argv[1]
    output_bed = sys.argv[2]
    
    print(f"Converting {input_gtf} to BED format...")
    gtf_to_bed(input_gtf, output_bed)
    print(f"Output written to {output_bed}")
