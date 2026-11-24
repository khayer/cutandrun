// ---------------------------------------------------------------------------
// Module: SPECIES_COMPOSITION
// Authors: Katharina E. Hayer, co-created with Claude
// Version: 1.0.0
// Date:    2026-11-23
// Purpose: Aggregate samtools idxstats across samples into per-species counts,
//          a MultiQC bargraph and a two-panel figure.
//
// The combined reference carries decoy contigs (hRSV, mycoplasma) whose whole
// purpose is to be empty. A 100% stacked bar cannot show the difference between
// 0.4% and 0.0% of a decoy, so the figure also plots non-host counts on a log
// axis where a trace is visible.
// ---------------------------------------------------------------------------

process SPECIES_COMPOSITION {
    label 'process_single'

    // deeptools ships numpy and matplotlib, and this image is already used by
    // the nf-core deeptools modules in this repo, so it is cached and known to
    // pull. Do not swap it for a mulled hash without checking the hash exists.
    conda "bioconda::deeptools=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'biocontainers/deeptools:3.5.1--py_0' }"

    input:
    path idxstats
    val  groups
    val  host_label
    val  viral_labels
    val  title

    output:
    path "species_composition.tsv"    , emit: tsv
    path "species_composition_mqc.tsv", emit: mqc
    path "species_composition.png"    , emit: png
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/species_composition.py \\
        ${idxstats} \\
        --groups '${groups}' \\
        --host '${host_label}' \\
        --viral '${viral_labels}' \\
        --title '${title}' \\
        --prefix species_composition

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
