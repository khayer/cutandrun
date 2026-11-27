/*
 * Module: MAKE_HOMER_BACKGROUND
 * Authors: Katharina E. Hayer, co-created with Claude
 * Version: 1.1.0  |  Date: 2026-11-27
 *
 * Builds a region-matched background for findMotifsGenome.pl -bg.
 *
 * HOMER's automatic background is GC-matched but drawn from the GENOME, so a
 * motif can be reporting "peak-like region versus genome" rather than anything
 * about the factor. On the eCLIP side, matching the background to the same
 * region class cut the target-vs-background GC gap from 0.106 to 0.035 and
 * changed which motifs came out; see homer_motif_background.md, panels B and C.
 *
 * The background is "assayable but not called": sequence covered in at least
 * one CONTROL library, with every peak removed. Sequence no library covers
 * could never have been a peak, so including it reintroduces the genome-wide
 * comparison. Control libraries are used because they carry no target signal,
 * and because DEDUPLICATE_PICARD has already removed their duplicates.
 *
 * Coverage is unioned across controls rather than pooled: a region counts as
 * assayable if any one control reached the threshold there.
 */
process MAKE_HOMER_BACKGROUND {
    label 'process_medium'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    path control_bams
    path peak_beds
    val  min_cov
    val  width
    val  exclude_regex

    output:
    path "homer_background.bed", emit: bed
    path "homer_background.log", emit: log
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    for b in ${control_bams}; do
        bedtools genomecov -ibam \$b -bg \\
            | awk -v m=${min_cov} 'BEGIN{FS=OFS="\\t"} \$4 >= m {print \$1,\$2,\$3}'
    done | sort -k1,1 -k2,2n | bedtools merge -i - > covered.bed

    cat ${peak_beds} \\
        | awk 'BEGIN{FS=OFS="\\t"} NF>=3 && \$1 !~ /^#/ {print \$1,\$2,\$3}' \\
        | sort -k1,1 -k2,2n | bedtools merge -i - > allpeaks.bed

    bedtools subtract -a covered.bed -b allpeaks.bed \\
        | awk -v ex='${exclude_regex}' 'BEGIN{FS=OFS="\\t"} \$1 !~ ex && (\$3-\$2) >= 50' \\
        > pool.bed

    if [ ! -s pool.bed ]; then
        echo "ERROR: background pool is empty after removing peaks and excluded contigs" >&2
        exit 1
    fi

    bedtools makewindows -b pool.bed -w ${width} \\
        | awk -v w=${width} 'BEGIN{FS=OFS="\\t"} (\$3-\$2) >= w*0.8 {print \$1,\$2,\$3,"bg_"NR,0,"."}' \\
        > homer_background.bed

    # HOMER bins target and background together for GC normalisation; a
    # background much smaller than the target set leaves empty bins and
    # findKnownMotifs.pl dies with "Illegal division by zero". Observed
    # 2026-11-27: 334 background regions against 30,383 peaks. Fail here,
    # where the message is readable, not 18 minutes into motif finding.
    N_BG=\$(wc -l < homer_background.bed)
    N_PK=\$(wc -l < allpeaks.bed)
    MIN_BG=\$(( N_PK * 5 ))
    if [ "\$N_BG" -lt "\$MIN_BG" ]; then
        echo "ERROR: background too small for HOMER: \$N_BG regions against \$N_PK peaks (need >= \$MIN_BG)." >&2
        echo "       Lower --homer_background_min_cov (currently ${min_cov}; shallow IgG" >&2
        echo "       controls rarely reach 5x), or use --homer_background_mode genomic." >&2
        exit 1
    fi

    {
        echo "min coverage      : ${min_cov}"
        echo "window width      : ${width}"
        echo "excluded contigs  : ${exclude_regex}"
        echo "covered intervals : \$(wc -l < covered.bed)"
        echo "peaks removed     : \$(wc -l < allpeaks.bed)"
        echo "pool intervals    : \$(wc -l < pool.bed)"
        echo "background regions: \$(wc -l < homer_background.bed)"
    } > homer_background.log
    cat homer_background.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
