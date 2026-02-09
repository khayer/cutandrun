/*
 * Create motif comparison tables across conditions
 * 
 * Author: Katharina Hayer
 * Co-created with: GitHub Copilot (Claude Sonnet 4.5)
 */

process CREATE_MOTIF_COMPARISON_TABLES {
    tag "motif_comparison"
    label 'process_low'

    conda "conda-forge::python=3.8.3 conda-forge::pandas=1.2.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path(motif_dirs)

    output:
    path("Known_Motifs_Comparison_Table.tsv") , emit: known_table
    path("DeNovo_Motifs_Comparison_Table.tsv"), emit: denovo_table
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Create directory structure expected by the script
    mkdir -p homer_motifs/consensus_peaks
    mkdir -p homer_motifs/merged_peaks
    
    # Link all motif directories to appropriate locations
    for dir in ${motif_dirs}; do
        if [[ \$dir == *"merged_peaks_motifs"* ]]; then
            ln -s "\$(readlink -f \$dir)" homer_motifs/merged_peaks/merged_peaks_motifs
        else
            # Extract condition name from directory (e.g., DRB_RI_26_motifs)
            condition=\$(basename \$dir)
            ln -s "\$(readlink -f \$dir)" homer_motifs/consensus_peaks/\${condition}
        fi
    done
    
    # Run the comparison table script
    create_motif_comparison_tables.py homer_motifs/
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
