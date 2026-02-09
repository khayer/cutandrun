// Module: SUMMARIZE_HOMER_MOTIFS
// Author: Katharina Hayer
// Co-created with: GitHub Copilot (Claude Sonnet 4.5)
// Purpose: Generate comprehensive motif analysis reports and visualizations

process SUMMARIZE_HOMER_MOTIFS {
    tag "summarize_motifs"
    label 'process_single'

    conda "conda-forge::python=3.8.3 conda-forge::pandas=1.2.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path motif_dirs

    output:
    path "Known_Motifs_Summary.txt", emit: known_summary
    path "DeNovo_Motifs_Summary.txt", emit: denovo_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat > summarize_motifs.py << 'EOF'
    #!/usr/bin/env python3
    
    import os
    import re
    from pathlib import Path
    from collections import defaultdict
    
    def parse_known_results(file_path):
        # Parse knownResults.txt and extract top 10 motifs
        motifs = []
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                for line in lines[1:11]:  # Top 10 motifs
                    parts = line.strip().split('\\t')
                    if len(parts) >= 7:
                        motif_name = parts[0].split('(')[0].strip()
                        consensus = parts[1]
                        p_value = parts[2]
                        pct_target = parts[6].rstrip('%')
                        motifs.append({
                            'Motif': motif_name,
                            'Consensus': consensus,
                            'P-value': p_value,
                            '% Target': pct_target
                        })
        except Exception as e:
            print(f"Error parsing {file_path}: {e}")
        return motifs
    
    def parse_denovo_motifs(motif_dir):
        # Parse de novo motifs from homerResults directory
        motifs = []
        homer_results = motif_dir / 'homerResults'
        
        if not homer_results.exists() or not homer_results.is_dir():
            return motifs
        
        motif_files = sorted([f for f in homer_results.glob('motif*.motif') 
                             if 'RV' not in f.name and 'similar' not in f.name],
                            key=lambda x: int(re.search(r'motif(\\d+)', x.name).group(1)))
        
        for motif_file in motif_files[:10]:  # Top 10
            try:
                with open(motif_file, 'r') as f:
                    first_line = f.readline().strip()
                    parts = first_line.lstrip('>').split('\\t')
                    if len(parts) < 6:
                        continue
                    
                    consensus = parts[0]
                    motif_num = motif_file.stem.replace('motif', '')
                    motif_name = f"Motif{motif_num}"
                    
                    stats = parts[5]
                    target_match = re.search(r'T:[0-9.]+\\(([0-9.]+)%\\)', stats)
                    pvalue_match = re.search(r'P:([0-9e\\-]+)', stats)
                    
                    percent_target = target_match.group(1) if target_match else "0.00"
                    pvalue = pvalue_match.group(1) if pvalue_match else parts[2]
                    
                    motifs.append({
                        'Motif': motif_name,
                        'Consensus': consensus,
                        'P-value': pvalue,
                        '% Target': percent_target
                    })
            except Exception:
                continue
        
        return motifs
    
    def create_report(results_dict, report_type):
        # Create summary report
        lines = []
        lines.append("=" * 100)
        if report_type == 'known':
            lines.append("HOMER MOTIF ANALYSIS - KNOWN MOTIFS SUMMARY")
        else:
            lines.append("HOMER MOTIF ANALYSIS - DE NOVO MOTIFS SUMMARY")
        lines.append("=" * 100)
        lines.append("")
        lines.append("Top 10 motifs for each peak set")
        lines.append("")
        
        for group_name in sorted(results_dict.keys()):
            motifs = results_dict[group_name]
            lines.append("")
            lines.append(f"{'='*80}")
            lines.append(f"Peak Set: {group_name}")
            lines.append(f"{'='*80}")
            lines.append("")
            
            if motifs:
                lines.append(f"  {'Rank':<6} {'Motif':<40} {'Consensus':<20} {'P-value':<15} {'% Target':<10}")
                lines.append(f"  {'-'*6} {'-'*40} {'-'*20} {'-'*15} {'-'*10}")
                for i, motif in enumerate(motifs, 1):
                    lines.append(f"  {i:<6} {motif['Motif']:<40} {motif['Consensus']:<20} {motif['P-value']:<15} {motif['% Target']:<10}")
            else:
                lines.append("  No motifs found")
            lines.append("")
        
        return lines
    
    def compare_motifs(results_dict, report_type):
        # Compare motifs across groups
        lines = []
        lines.append("")
        lines.append("=" * 100)
        lines.append("MOTIF COMPARISON ACROSS GROUPS")
        lines.append("=" * 100)
        lines.append("")
        
        # Separate merged peaks from consensus peaks
        merged_motifs = {}
        consensus_motifs = {}
        
        for group_name, motifs in results_dict.items():
            motif_set = set()
            for motif in motifs:
                if report_type == 'denovo':
                    motif_id = f"{motif['Motif']} ({motif['Consensus']})"
                else:
                    motif_id = motif['Motif']
                motif_set.add(motif_id)
            
            if 'merged' in group_name.lower():
                merged_motifs[group_name] = motif_set
            else:
                consensus_motifs[group_name] = motif_set
        
        # Compare consensus peaks groups
        if len(consensus_motifs) > 1:
            lines.append("Consensus Peaks Comparison:")
            lines.append("-" * 80)
            
            groups = sorted(consensus_motifs.keys())
            for i in range(len(groups)):
                for j in range(i+1, len(groups)):
                    group1 = groups[i]
                    group2 = groups[j]
                    
                    common = consensus_motifs[group1] & consensus_motifs[group2]
                    unique1 = consensus_motifs[group1] - consensus_motifs[group2]
                    unique2 = consensus_motifs[group2] - consensus_motifs[group1]
                    
                    lines.append(f"\\n{group1} vs {group2}:")
                    lines.append(f"  Shared motifs: {len(common)}")
                    if common:
                        for motif in sorted(list(common)[:5]):
                            lines.append(f"    - {motif}")
                    
                    lines.append(f"  Unique to {group1}: {len(unique1)}")
                    if unique1:
                        for motif in sorted(list(unique1)[:3]):
                            lines.append(f"    - {motif}")
                    
                    lines.append(f"  Unique to {group2}: {len(unique2)}")
                    if unique2:
                        for motif in sorted(list(unique2)[:3]):
                            lines.append(f"    - {motif}")
        
        # Compare merged vs consensus groups
        if merged_motifs and consensus_motifs:
            lines.append("\\n")
            lines.append("=" * 80)
            lines.append("Merged Peaks vs Consensus Peaks:")
            lines.append("-" * 80)
            
            merged_all = set()
            for motifs in merged_motifs.values():
                merged_all.update(motifs)
            
            for group_name, consensus_set in sorted(consensus_motifs.items()):
                common = merged_all & consensus_set
                only_merged = merged_all - consensus_set
                only_consensus = consensus_set - merged_all
                
                lines.append(f"\\nMerged peaks vs {group_name}:")
                lines.append(f"  Shared motifs: {len(common)}")
                lines.append(f"  Only in merged: {len(only_merged)}")
                lines.append(f"  Only in {group_name}: {len(only_consensus)}")
        
        lines.append("")
        return lines
    
    # Main execution
    print("Collecting Homer motif results...")
    
    known_results = {}
    denovo_results = {}
    
    # Search for motif directories
    for item in Path('.').glob('**/*_motifs'):
        if not item.is_dir():
            continue
        
        group_name = item.name.replace('_motifs', '')
        
        # Parse known motifs
        known_file = item / 'knownResults.txt'
        if known_file.exists():
            motifs = parse_known_results(known_file)
            if motifs:
                known_results[group_name] = motifs
        
        # Parse de novo motifs
        denovo_motifs = parse_denovo_motifs(item)
        if denovo_motifs:
            denovo_results[group_name] = denovo_motifs
    
    # Generate reports
    if known_results:
        known_report = create_report(known_results, 'known')
        known_comparison = compare_motifs(known_results, 'known')
        
        with open('Known_Motifs_Summary.txt', 'w') as f:
            f.write('\\n'.join(known_report + known_comparison))
        
        print(f"Created Known_Motifs_Summary.txt with {len(known_results)} groups")
    else:
        # Create empty file
        with open('Known_Motifs_Summary.txt', 'w') as f:
            f.write("No known motif results found\\n")
    
    if denovo_results:
        denovo_report = create_report(denovo_results, 'denovo')
        denovo_comparison = compare_motifs(denovo_results, 'denovo')
        
        with open('DeNovo_Motifs_Summary.txt', 'w') as f:
            f.write('\\n'.join(denovo_report + denovo_comparison))
        
        print(f"Created DeNovo_Motifs_Summary.txt with {len(denovo_results)} groups")
    else:
        # Create empty file
        with open('DeNovo_Motifs_Summary.txt', 'w') as f:
            f.write("No de novo motif results found\\n")
    
    print("Motif summarization complete")
    EOF
    
    python summarize_motifs.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
