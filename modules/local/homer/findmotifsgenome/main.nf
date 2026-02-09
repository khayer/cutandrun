// Module: HOMER_FINDMOTIFSGENOME
// Author: Katharina Hayer
// Co-created with: GitHub Copilot (Claude Sonnet 4.5)
// Purpose: Homer motif discovery and enrichment analysis

process HOMER_FINDMOTIFSGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::homer=4.11"
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(bed)
    path fasta
    val size

    output:
    tuple val(meta), path("${prefix}/"), emit: motifs
    tuple val(meta), path("${prefix}/homerResults.html"), emit: html
    tuple val(meta), path("${prefix}/knownResults.txt"), emit: known_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    findMotifsGenome.pl \\
        $bed \\
        $fasta \\
        ${prefix}/ \\
        -size $size \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(homer2 -h 2>&1) | grep -o 'v[0-9.]*' | sed 's/v//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/
    touch ${prefix}/homerResults.html
    touch ${prefix}/knownResults.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(homer2 -h 2>&1) | grep -o 'v[0-9.]*' | sed 's/v//')
    END_VERSIONS
    """
}
