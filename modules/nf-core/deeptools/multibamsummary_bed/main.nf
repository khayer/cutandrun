process DEEPTOOLS_MULTIBAMSUMMARY_BED {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::deeptools=3.5.1"
    container "biocontainers/deeptools:3.5.1--py_0"

    input:
    tuple val(meta), path(bams), path(bais), val(labels)
    path bed

    output:
    tuple val(meta), path("*.npz")     , emit: matrix
    tuple val(meta), path("*.tab")     , emit: table
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def label = labels ? "--labels ${labels.join(' ')}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    multiBamSummary BED-file \\
        $args \\
        $label \\
        --BED $bed \\
        --bamfiles ${bams.join(' ')} \\
        --numberOfProcessors $task.cpus \\
        --outFileName ${prefix}.npz \\
        --outRawCounts ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBamSummary --version | sed -e "s/multiBamSummary //g")
    END_VERSIONS
    """
}
