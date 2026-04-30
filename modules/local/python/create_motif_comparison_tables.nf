/*
 * Create motif comparison tables across conditions
 * 
 * Author: Katharina Hayer
 * Co-created with: GitHub Copilot (Claude Sonnet 4.5)
 */

process CREATE_MOTIF_COMPARISON_TABLES {
    tag "motif_comparison"
    label 'process_low'

    // Use conda specification to support local/conda runs and conda-to-container resolution.
    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::matplotlib=3.8.4 conda-forge::cairosvg=2.7.1"

    input:
    // Accept as individual paths in a collection, each will be available as a variable
    val(motif_dirs)  // This will be a list of directory paths

    output:
    path("Known_Motifs_Comparison_Table.tsv") , emit: known_table
    path("DeNovo_Motifs_Comparison_Table.tsv"), emit: denovo_table
    path("*_Motifs_Table.pdf")               , emit: motif_pdfs, optional: true
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Create directory structure expected by the script
    task_workdir="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"

    mkdir -p homer_motifs/consensus_peaks
    mkdir -p homer_motifs/merged_peaks

    # Process motif directories
    echo "Processing motif directories..."
    found_count=0
    
    # Read directories from input (motif_dirs is a list of paths)
    cat > /tmp/motif_dirs.txt <<'EOF'
${motif_dirs.join('\n')}
EOF
    
    while IFS= read -r dir; do
        # Skip empty lines
        [[ -z "\$dir" ]] && continue
        
        # Handle case where path might not exist
        if [[ ! -d "\$dir" ]]; then
            echo "WARNING: Motif directory does not exist or is not a directory: \$dir"
            continue
        fi
        
        echo "Processing motif directory: \$dir"
        found_count=\$((found_count + 1))
        
        # Resolve the real absolute path
        real_path="\$(cd "\$dir" 2>/dev/null && pwd)" || {
            echo "WARNING: Failed to resolve path for \$dir, skipping..."
            continue
        }
        
        if [[ "\$dir" == *"merged_peaks_motifs"* ]]; then
            echo "Linking merged peaks: \$real_path"
            ln -sfn "\$real_path" homer_motifs/merged_peaks/merged_peaks_motifs
        else
            # Extract condition name from directory (e.g., DRB_RI_26_motifs)
            condition=\$(basename "\$dir")
            echo "Linking consensus peaks (\$condition): \$real_path"
            ln -sfn "\$real_path" homer_motifs/consensus_peaks/\${condition}
        fi
    done < /tmp/motif_dirs.txt
    
    if [[ \$found_count -eq 0 ]]; then
        echo "ERROR: No valid motif directories found!"
        echo "This may indicate that HOMER processes did not produce output."
        exit 1
    fi
    
    echo "Successfully processed \$found_count motif directories"
    
    # Run the comparison table script
    create_motif_comparison_tables.py homer_motifs/ 5
    
    # Copy output files to the true task workdir where Nextflow validates outputs
    cp homer_motifs/Known_Motifs_Comparison_Table.tsv "\${task_workdir}/"
    cp homer_motifs/DeNovo_Motifs_Comparison_Table.tsv "\${task_workdir}/"
    cp homer_motifs/*_Motifs_Table.pdf "\${task_workdir}/" 2>/dev/null || true
    
    cat <<-END_VERSIONS > "\${task_workdir}/versions.yml"
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        cairosvg: \$(python -c "import importlib.util as u; m=u.find_spec('cairosvg'); print(__import__('cairosvg').__version__ if m else 'not_available')" 2>/dev/null || echo 'not_available')
    END_VERSIONS
    """
}
