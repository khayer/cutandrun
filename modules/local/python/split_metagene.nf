// ---------------------------------------------------------------------------
// Module: SPLIT_METAGENE
// Authors: Katharina E. Hayer, co-created with Claude
// Version: 1.1.0
// Date:    2026-11-24
// Purpose: TSS- and TES-anchored computeMatrix over every replicate bigWig of
//          one group and two gene groups, then the split metagene figure.
//
// The standalone runner globs a directory for <cond>_R*.<type>.bigWig. Here the
// replicates arrive already staged as a channel, so the labels are derived from
// the staged file names instead and no globbing of a results folder happens.
// ---------------------------------------------------------------------------

process SPLIT_METAGENE {
    tag "${meta.group}"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bigwigs)
    path bed_intron_containing
    path bed_intronless
    val  signal
    val  ci
    val  bin_size
    val  tss_up
    val  tss_down
    val  tes_up
    val  tes_down

    output:
    tuple val(meta), path("${meta.group}.split_metagene.png"), emit: png
    tuple val(meta), path("${meta.group}.split_metagene.pdf"), emit: pdf,  optional: true
    tuple val(meta), path("${meta.group}.*.mat.gz")          , emit: matrices
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bw = bigwigs instanceof List ? bigwigs : [bigwigs]
    def bw_list = bw.join(' ')
    // A Groovy slashy string does not treat \ as an escape, so /\\./ is a
    // literal backslash followed by any character and strips nothing.
    def labels = bw.collect { it.getName().replaceAll(/\.(bigWig|bw)$/, '') }.join(' ')
    """
    set -euo pipefail

    # Counts read to EOF deliberately: a presence test through a pipe whose
    # reader exits early inverts under `set -o pipefail`.
    n_ic=\$(wc -l < ${bed_intron_containing})
    n_il=\$(wc -l < ${bed_intronless})
    if [ "\$n_ic" -eq 0 ] || [ "\$n_il" -eq 0 ]; then
        echo "ERROR: a gene group is empty (intron-containing \$n_ic, intronless \$n_il)" >&2
        exit 1
    fi

    computeMatrix reference-point --referencePoint TSS \\
        -b ${tss_up} -a ${tss_down} \\
        -S ${bw_list} -R ${bed_intron_containing} ${bed_intronless} \\
        --binSize ${bin_size} --missingDataAsZero \\
        --samplesLabel ${labels} -p ${task.cpus} \\
        -o ${meta.group}.TSS.mat.gz

    computeMatrix reference-point --referencePoint TES \\
        -b ${tes_up} -a ${tes_down} \\
        -S ${bw_list} -R ${bed_intron_containing} ${bed_intronless} \\
        --binSize ${bin_size} --missingDataAsZero \\
        --samplesLabel ${labels} -p ${task.cpus} \\
        -o ${meta.group}.TES.mat.gz

    python3 ${workflow.projectDir}/split_metagene/plot_split_metagene.py \\
        --tss ${meta.group}.TSS.mat.gz \\
        --tes ${meta.group}.TES.mat.gz \\
        --signal ${signal} \\
        --ci ${ci} \\
        --title "${meta.group}" \\
        --out ${meta.group}.split_metagene

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
    END_VERSIONS
    """
}
