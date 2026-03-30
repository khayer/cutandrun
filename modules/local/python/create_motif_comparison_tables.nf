/*
 * Create motif comparison tables across conditions
 * 
 * Author: Katharina Hayer
 * Co-created with: GitHub Copilot (Claude Sonnet 4.5)
 */

process CREATE_MOTIF_COMPARISON_TABLES {
    tag "motif_comparison"
    label 'process_low'

    // Run this step with conda so SVG logos can be rendered via Python cairosvg.
    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::matplotlib=3.8.4 conda-forge::cairosvg=2.7.1"

    input:
    path(motif_dirs, stageAs: 'motif_dirs/*')

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
    
    # Link all motif directories to appropriate locations
    for dir in motif_dirs/*; do
        if [[ \$dir == *"merged_peaks_motifs"* ]]; then
            ln -sfn "\$(readlink -f \$dir)" homer_motifs/merged_peaks/merged_peaks_motifs
        else
            # Extract condition name from directory (e.g., DRB_RI_26_motifs)
            condition=\$(basename \$dir)
            ln -sfn "\$(readlink -f \$dir)" homer_motifs/consensus_peaks/\${condition}
        fi
    done
    
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
    END_VERSIONS
    """
}
