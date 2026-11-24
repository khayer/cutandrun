// ---------------------------------------------------------------------------
// Module: MAKE_GENE_GROUPS
// Authors: Katharina E. Hayer, co-created with Claude
// Version: 1.0.0
// Date:    2026-11-23
// Purpose: Split the run's own GTF into intron-containing and intronless gene
//          BED6 files for the split metagene.
//
// "Intronless" counts `exon` features, not CDS: a gene is intronless only when
// EVERY transcript is a single exon. Counting CDS exons instead would call a
// gene intronless whenever its coding region is contiguous but its 5'UTR is
// spliced, which is common and which is still a spliced gene.
// ---------------------------------------------------------------------------

process MAKE_GENE_GROUPS {
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path gtf
    val  min_length
    val  biotype
    val  exclude_contigs

    output:
    path "genes.intron_containing.bed", emit: intron_containing
    path "genes.intronless.bed"       , emit: intronless
    path "gene_groups.counts.txt"     , emit: counts
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/split_metagene/make_gene_groups.py \\
        ${gtf} \\
        --out-prefix genes \\
        --min-length ${min_length} \\
        --biotype ${biotype} \\
        --exclude-contigs '${exclude_contigs}' \\
        2> gene_groups.counts.txt

    cat gene_groups.counts.txt >&2

    # The intronless set is small by construction (min-length plus the
    # all-transcripts-single-exon rule). Warn loudly rather than let a
    # three-gene average be read as a class difference.
    n_il=\$(wc -l < genes.intronless.bed)
    n_ic=\$(wc -l < genes.intron_containing.bed)
    echo "intronless genes: \$n_il    intron-containing genes: \$n_ic" >&2
    if [ "\$n_il" -lt 100 ]; then
        echo "WARNING: only \$n_il intronless genes passed the filters." >&2
        echo "WARNING: the intronless profile is an average over \$n_il genes;" >&2
        echo "WARNING: differences between the two curves may be gene sampling." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
