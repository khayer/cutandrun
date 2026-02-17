process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_single'
    label 'error_ignore'

    conda "bioconda::preseq=3.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.1.2--h445547b_2':
        'biocontainers/preseq:3.1.2--h445547b_2' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.lc_extrap.txt"), emit: lc_extrap
    tuple val(meta), path("*.log")          , emit: log
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Respect any user-supplied args, but by default attempt a high-accuracy run
    // and automatically fall back to defect mode (-D) if the estimator fails.
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    def args_str = args.toString()
    def has_defect = args_str.contains('-D')
    def has_user_highacc = args_str.contains('-x') || args_str.contains('-n') || args_str.contains('-e')
    def highacc_defaults = '-v -B -P -x 200 -n 1000 -e 100000000'

    """
    if [ "${has_defect}" = "true" ]; then
        # User requested defect mode explicitly — run once with provided args
        preseq lc_extrap \
            $args \
            $paired_end \
            -output ${prefix}.lc_extrap.txt \
            $bam > ${prefix}.preseq.log 2>&1 || true

        cp .command.err ${prefix}.command.log || true

    else

        if [ "${has_user_highacc}" = "true" ]; then
            # User supplied their own high-accuracy flags: attempt once, then fallback to -D if needed
            preseq lc_extrap \
                $args \
                $paired_end \
                -output ${prefix}.lc_extrap.txt \
                $bam > ${prefix}.preseq_highacc.log 2>&1 || true

            if grep -q "too many defects" ${prefix}.preseq_highacc.log || [ ! -s ${prefix}.lc_extrap.txt ]; then
                echo "High-accuracy run failed or defect test failed — rerunning in defect mode (-D)" >> ${prefix}.preseq_highacc.log
                preseq lc_extrap \
                    $args -D \
                    $paired_end \
                    -output ${prefix}.lc_extrap.txt \
                    $bam > ${prefix}.preseq_highacc_defect.log 2>&1 || true
                cp ${prefix}.preseq_highacc_defect.log ${prefix}.command.log || true
            else
                cp ${prefix}.preseq_highacc.log ${prefix}.command.log || true
            fi

        else
            # Default behaviour: try a high-accuracy attempt (our defaults), then -D if it fails
            preseq lc_extrap \
                $args ${highacc_defaults} \
                $paired_end \
                -output ${prefix}.lc_extrap.highacc.txt \
                $bam > ${prefix}.preseq_highacc.log 2>&1 || true

            if [ -s ${prefix}.lc_extrap.highacc.txt ] && ! grep -q "too many defects" ${prefix}.preseq_highacc.log; then
                mv ${prefix}.lc_extrap.highacc.txt ${prefix}.lc_extrap.txt
                cp ${prefix}.preseq_highacc.log ${prefix}.command.log || true
            else
                echo "High-accuracy attempt failed or defect test failed — rerunning with -D (defect mode)" >> ${prefix}.preseq_highacc.log
                preseq lc_extrap \
                    $args -D ${highacc_defaults} \
                    $paired_end \
                    -output ${prefix}.lc_extrap.txt \
                    $bam > ${prefix}.preseq_highacc_defect.log 2>&1 || true
                cp ${prefix}.preseq_highacc_defect.log ${prefix}.command.log || true
            fi
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
}
