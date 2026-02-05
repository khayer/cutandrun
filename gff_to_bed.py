#!/usr/bin/env python3
"""
Convert GFF file to BED12 format.
Groups exons/CDS by transcript/parent and creates proper BED12 entries with blocks.
"""

import re
import sys
from collections import defaultdict

def parse_gff_attributes(attr_string):
    """Parse GFF attribute string and extract gene name, ID, and Parent."""
    result = {}
    
    # Look for ID=... pattern
    id_match = re.search(r'ID=([^;]+)', attr_string)
    if id_match:
        result['id'] = id_match.group(1)
    
    # Look for gene=... pattern
    gene_match = re.search(r'gene=([^;]+)', attr_string)
    if gene_match:
        result['gene'] = gene_match.group(1)
    
    # Look for Parent=... pattern
    parent_match = re.search(r'Parent=([^;]+)', attr_string)
    if parent_match:
        result['parent'] = parent_match.group(1)
    
    return result

def gff_to_bed(gff_file, bed_file):
    """
    Convert GFF to BED12 format with exon IDs and parent ID in columns 13-14.
    
    Parameters:
    - gff_file: Input GFF file path
    - bed_file: Output BED file path
    """
    # Dictionary to store gene/mRNA/transcript information
    genes = {}
    # Dictionary to group CDS/exons by parent (mRNA/transcript)
    transcripts = defaultdict(list)
    
    with open(gff_file, 'r') as infile:
        for line in infile:
            # Skip comments and empty lines
            if line.startswith('#') or line.strip() == '':
                continue
            
            fields = line.strip().split('\t')
            
            # We need at least 9 fields for a valid GFF
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            source = fields[1]
            feature_type = fields[2]
            start = int(fields[3])  # GFF is 1-based
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            frame = fields[7]
            attributes = fields[8]
            
            # Extract attributes
            attrs = parse_gff_attributes(attributes)
            
            # Store gene information
            if feature_type == 'gene':
                gene_id = attrs.get('id', '')
                genes[gene_id] = {
                    'chrom': chrom,
                    'start': start - 1,  # Convert to 0-based
                    'end': end,
                    'name': gene_id,
                    'score': score if score != '.' else '0',
                    'strand': strand
                }
            
            # Store mRNA information (will be used if no CDS exists)
            elif feature_type == 'mRNA':
                parent_id = attrs.get('id', '')
                gene_name = attrs.get('gene', parent_id)
                genes[parent_id] = {
                    'chrom': chrom,
                    'start': start - 1,
                    'end': end,
                    'name': gene_name,
                    'score': score if score != '.' else '0',
                    'strand': strand,
                    'parent_id': parent_id
                }
            
            # Process CDS features (primary feature for coding regions)
            elif feature_type == 'CDS':
                gene_name = attrs.get('gene', '')
                parent_id = attrs.get('gene', '')  # Use gene name as parent for simple genes
                cds_id = attrs.get('id', f'cds-{gene_name}')
                
                # If no gene name, skip
                if not gene_name:
                    continue
                
                # Store CDS information grouped by gene/parent
                transcripts[parent_id].append({
                    'chrom': chrom,
                    'start': start - 1,  # Convert to 0-based
                    'end': end,
                    'gene': gene_name,
                    'score': score if score != '.' else '0',
                    'strand': strand,
                    'cds_id': cds_id,
                    'feature_type': feature_type
                })
    
    # Write BED12 format
    with open(bed_file, 'w') as outfile:
        for parent_id, features in sorted(transcripts.items()):
            # Sort features by start position
            features.sort(key=lambda x: x['start'])
            
            # Get transcript bounds
            chrom = features[0]['chrom']
            tx_start = min(f['start'] for f in features)
            tx_end = max(f['end'] for f in features)
            gene_name = features[0]['gene']
            score = features[0]['score']
            strand = features[0]['strand']
            
            # Calculate block information
            block_count = len(features)
            block_sizes = ','.join(str(f['end'] - f['start']) for f in features) + ','
            block_starts = ','.join(str(f['start'] - tx_start) for f in features) + ','
            
            # Collect all CDS IDs for this transcript
            cds_ids = ','.join(f['cds_id'] for f in features)
            
            # BED12 format: chrom, chromStart, chromEnd, name, score, strand, 
            #               thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
            # Plus columns 13-14: CDS IDs, parent ID (gene name)
            bed_line = (f'{chrom}\t{tx_start}\t{tx_end}\t{gene_name}\t{score}\t{strand}\t'
                       f'{tx_start}\t{tx_end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\t'
                       f'{cds_ids}\tgene-{parent_id}\n')
            outfile.write(bed_line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python gff_to_bed.py <input.gff> <output.bed>")
        sys.exit(1)
    
    input_gff = sys.argv[1]
    output_bed = sys.argv[2]
    
    print(f"Converting {input_gff} to BED format...")
    gff_to_bed(input_gff, output_bed)
    print(f"Output written to {output_bed}")
