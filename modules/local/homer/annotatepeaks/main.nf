// Module: HOMER_ANNOTATEPEAKS
// Author: Katharina Hayer
// Co-created with: GitHub Copilot (GPT-5.3-Codex)
// Purpose: Annotate replicate peak BED files with HOMER annotatePeaks.pl

process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::homer=4.11"
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(bed)
    path fasta
    path gtf
    val tss_dist

    output:
    tuple val(meta), path("${prefix}.annotatePeaks.txt"), emit: annot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotatePeaks.pl \
        $bed \
        $fasta \
        -gtf $gtf \
        -size given \
        -d $tss_dist \
        > ${prefix}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(homer2 -h 2>&1) | grep -o 'v[0-9.]*' | sed 's/v//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'EOF' > ${prefix}.annotatePeaks.txt
PeakID\tChr\tStart\tEnd\tStrand\tPeak Score\tFocus Ratio/Region Size\tAnnotation\tDetailed Annotation\tDistance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\tNearest Refseq\tNearest Ensembl\tGene Name\tGene Alias\tGene Description\tGene Type\tGC%\tCpG%\npeak_1\tchr1\t100\t200\t+\t10\t0\tIntergenic\tIntergenic\t10000\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t42.0\t2.0\nEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(homer2 -h 2>&1) | grep -o 'v[0-9.]*' | sed 's/v//')
    END_VERSIONS
    """
}
