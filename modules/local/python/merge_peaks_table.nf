process MERGE_PEAKS_TABLE {
    tag "merge_peaks_table"
    label 'process_medium'

    conda "conda-forge::python=3.8.3 conda-forge::pandas=1.2.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path peak_beds

    output:
    path "merged_peaks_table.txt", emit: table
    path "merged_peaks.bed"      , emit: bed
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bed_files = peak_beds instanceof List ? peak_beds.join(' ') : peak_beds
    """
    #!/usr/bin/env python3
    
    import os
    import sys
    from collections import defaultdict
    
    # Read all peak files
    peaks_by_sample = {}
    all_peaks = []
    
    bed_files = "${bed_files}".split()
    
    for bed_file in bed_files:
        sample = os.path.basename(bed_file).replace('.macs2.peaks.cut.bed', '').replace('.bed', '')
        peaks_by_sample[sample] = []
        
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\\t')
                if len(fields) < 3:
                    continue
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                peak = (chrom, start, end)
                peaks_by_sample[sample].append(peak)
                all_peaks.append(peak)
    
    # Merge overlapping peaks across all samples
    def merge_peaks(peaks):
        if not peaks:
            return []
        
        # Sort peaks
        sorted_peaks = sorted(peaks, key=lambda x: (x[0], x[1], x[2]))
        merged = [sorted_peaks[0]]
        
        for current in sorted_peaks[1:]:
            last = merged[-1]
            # If same chromosome and overlapping
            if current[0] == last[0] and current[1] <= last[2]:
                # Merge by extending the end position
                merged[-1] = (last[0], last[1], max(last[2], current[2]))
            else:
                merged.append(current)
        
        return merged
    
    merged_peaks = merge_peaks(all_peaks)
    
    # Check presence of each merged peak in each sample
    def has_overlap(peak, sample_peaks, min_overlap=0.5):
        chrom, start, end = peak
        peak_len = end - start
        
        for sp in sample_peaks:
            if sp[0] != chrom:
                continue
            # Calculate overlap
            overlap_start = max(start, sp[1])
            overlap_end = min(end, sp[2])
            if overlap_start < overlap_end:
                overlap_len = overlap_end - overlap_start
                # Require at least min_overlap fraction
                if overlap_len >= peak_len * min_overlap or overlap_len >= (sp[2] - sp[1]) * min_overlap:
                    return True
        return False
    
    # Write merged peaks BED
    with open('merged_peaks.bed', 'w') as out:
        for i, (chrom, start, end) in enumerate(merged_peaks, 1):
            out.write(f"{chrom}\\t{start}\\t{end}\\tpeak_{i}\\t0\\t.\\n")
    
    # Create presence/absence table
    samples = sorted(peaks_by_sample.keys())
    
    with open('merged_peaks_table.txt', 'w') as out:
        # Header
        header = ['PeakID', 'Chr', 'Start', 'End', 'Length'] + samples + ['Total']
        out.write('\\t'.join(header) + '\\n')
        
        # Data rows
        for i, (chrom, start, end) in enumerate(merged_peaks, 1):
            peak_id = f"peak_{i}"
            length = end - start
            
            presence = []
            total = 0
            for sample in samples:
                if has_overlap((chrom, start, end), peaks_by_sample[sample]):
                    presence.append('1')
                    total += 1
                else:
                    presence.append('0')
            
            row = [peak_id, chrom, str(start), str(end), str(length)] + presence + [str(total)]
            out.write('\\t'.join(row) + '\\n')
    
    # Write versions
    with open('versions.yml', 'w') as v:
        v.write('"${task.process}":\\n')
        v.write(f'    python: "{sys.version.split()[0]}"\\n')
    """
}
