// Module: PYGENOMETRACKS_TOP10
// Purpose: Render top consensus peak snapshots per group using pyGenomeTracks.

process PYGENOMETRACKS_TOP10 {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::pygenometracks=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pygenometracks:3.9--pyhdfd78af_0' :
        'biocontainers/pygenometracks:3.9--pyhdfd78af_0' }"

    publishDir "${params.outdir}/04_consensus_peaks/pygenometracks_top10", mode: 'copy', saveAs: { filename -> "${task.tag}/${filename}" }

    input:
    tuple val(meta), path(consensus_bed), path(bigwigs), val(own_count), val(has_control), path(chipseeker_annot), path(homer_annot)
    path gtf
    val top_n
    val flank
    val window_bp
    val rank_mode
    val feature_types
    val feature_anchor_window
    val output_format

    output:
    tuple val(meta), path("${meta.id}"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bigwig_list = bigwigs instanceof List ? bigwigs.join(' ') : bigwigs
    """
    set -euo pipefail

    mkdir -p ${meta.id}

    python ${workflow.projectDir}/bin/plot_top10_peaks_pygenometracks.py \
        --consensus-bed ${consensus_bed} \
        --gtf ${gtf} \
        --chipseeker-annot ${chipseeker_annot} \
        --homer-annot ${homer_annot} \
        --group-id ${meta.id} \
        --top-n ${top_n} \
        --flank ${flank} \
        --window-bp ${window_bp} \
        --rank-mode ${rank_mode} \
        --feature-types "${feature_types}" \
        --anchor-window ${feature_anchor_window} \
        --out-regions-tsv ${meta.id}/top_regions.tsv \
        --out-top-bed ${meta.id}/top_peaks.bed

    cat > ${meta.id}/tracks.ini <<'EOF'
[x-axis]
where = top

[genes]
file = ${gtf}
file_type = gtf
title = Genes
height = 2
renderGenes = yes
arrowLinkColor = black
arrowLinkStyle = solid

[consensus_peaks]
file = ${meta.id}/top_peaks.bed
file_type = bed
title = Consensus peaks
height = 1.2
color = black
EOF

    i=1
    own_count=${own_count}
    has_control=${has_control}
    total_count=\$(echo ${bigwig_list} | wc -w)
    for bw in ${bigwig_list}; do
        name=\$(basename "\$bw")
        name=\${name%.bigWig}
        name=\${name%.bw}
        color="#1f78b4"
        if (( has_control == 1 && i == total_count )); then
            color="#808080"
            name="\${name} (control rep)"
        elif (( i > own_count )); then
            color="#33a02c"
            name="\${name} (other rep)"
        fi
        cat >> ${meta.id}/tracks.ini <<EOF

[signal_\${i}]
file = \$bw
file_type = bigwig
title = \$name
height = 1.5
color = \$color
EOF
        i=\$((i+1))
    done

    while IFS="\$(printf '\t')" read -r rank chrom peak_start peak_end window_start window_end nearest_gene distance_bp title file_stub; do
        region="\${chrom}:\${window_start}-\${window_end}"
        if [[ "${output_format}" == "png" || "${output_format}" == "both" ]]; then
            pyGenomeTracks --tracks ${meta.id}/tracks.ini --region "\$region" --title "\$title" --dpi 150 --outFileName ${meta.id}/\${file_stub}.png
        fi
        if [[ "${output_format}" == "pdf" || "${output_format}" == "both" ]]; then
            pyGenomeTracks --tracks ${meta.id}/tracks.ini --region "\$region" --title "\$title" --dpi 150 --outFileName ${meta.id}/\${file_stub}.pdf
        fi
    done < <(tail -n +2 ${meta.id}/top_regions.tsv)

    pygt_version=\$(pyGenomeTracks --version 2>&1 | awk '{print \$NF}')
    printf '"%s":\n' "${task.process}" > versions.yml
    printf '    pyGenomeTracks: "%s"\n' "\$pygt_version" >> versions.yml
    """
}
