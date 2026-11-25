/*
 * Module: VIRAL_WINDOW_COUNTS
 * Authors: Katharina E. Hayer, co-created with Claude
 * Version: 1.0.0  |  Date: 2026-11-25
 *
 * Tiles one contig with fixed-width windows and counts reads per sample.
 *
 * bedtools multicov EXCLUDES duplicate-flagged reads unless -D is given
 * (verified: 10 unique + 30 duplicate-flagged pairs count as 20 reads without
 * -D and 80 with it). That matters here: with dedup_target_reads = false the
 * target BAMs keep their duplicates and the control BAMs have had theirs
 * removed, so counting all alignments would compare a duplicate-inflated
 * target against a deduplicated control. The default is what makes the two
 * sides comparable, so do not add -D.
 */
process VIRAL_WINDOW_COUNTS {
    tag "$contig"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    val  contig
    path fai
    path bams
    path bais
    val  sample_names
    val  sample_table
    val  window
    val  step
    val  min_mapq

    output:
    tuple val(contig), path("${contig}.counts.tsv"), path("${contig}.samples.tsv"), emit: counts
    path "versions.yml"                                                           , emit: versions

    script:
    def bam_list = bams.join(' ')
    def header   = (['chrom', 'start', 'end'] + sample_names).join('\t')
    """
    set -euo pipefail

    awk -v c="${contig}" 'BEGIN{FS=OFS="\\t"} \$1==c {print \$1, \$2}' ${fai} > ${contig}.genome
    if [ ! -s ${contig}.genome ]; then
        echo "ERROR: contig ${contig} is not in ${fai}" >&2
        exit 1
    fi

    bedtools makewindows -g ${contig}.genome -w ${window} -s ${step} > ${contig}.windows.bed

    # -p counts proper pairs only; duplicates are excluded by default (see header)
    printf '%s\\n' '${header}' > ${contig}.counts.tsv
    bedtools multicov -p -q ${min_mapq} -bams ${bam_list} -bed ${contig}.windows.bed \\
        >> ${contig}.counts.tsv

    cat > ${contig}.samples.tsv <<'SAMPLETABLE'
${sample_table}
SAMPLETABLE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
