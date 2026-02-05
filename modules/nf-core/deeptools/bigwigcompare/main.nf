process DEEPTOOLS_BIGWIGCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bigwig1), path(bigwig2)

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bigwigCompare \\
        --bigwig1 $bigwig1 \\
        --bigwig2 $bigwig2 \\
        --outFileName ${prefix}.log2ratio.bigWig \\
        --operation log2 \\
        --numberOfProcessors $task.cpus \\
        $args

    bigwigCompare \\
        --bigwig1 $bigwig1 \\
        --bigwig2 $bigwig2 \\
        --outFileName ${prefix}.subtract.bigWig \\
        --operation subtract \\
        --numberOfProcessors $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bigwigCompare --version | sed -e "s/bigwigCompare //g")
    END_VERSIONS
    """
}
