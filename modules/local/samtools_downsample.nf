// Module: DOWNSAMPLE_BAM
// Author: Katharina Hayer
// Co-created with: GitHub Copilot (Claude Sonnet 4.5)
// Purpose: Downsample BAMs to a target coverage depth

process DOWNSAMPLE_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path chrom_sizes
    val target_coverage
    val seed

    output:
    tuple val(meta), path("${meta.id}.downsampled.bam"), path("${meta.id}.downsampled.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    genome_size=\$(awk '{sum+=\$2} END {print sum+0}' ${chrom_sizes})
    mapped_reads=\$(samtools view -c -F 260 ${bam})
    avg_read_len=\$(samtools view ${bam} | head -n 10000 | awk '{sum+=length(\$10)} END {print (NR>0)?sum/NR:0}')

    fraction=\$(awk -v cov=${target_coverage} -v genome=\$genome_size -v mapped=\$mapped_reads -v readlen=\$avg_read_len 'BEGIN {
        if (mapped == 0 || readlen == 0) { print 1; exit }
        frac = (cov * genome) / (mapped * readlen)
        if (frac > 1) frac = 1
        if (frac < 0) frac = 0
        print frac
    }')

    should_downsample=\$(awk -v f=\$fraction 'BEGIN { print (f < 0.999999) ? 1 : 0 }')

    if [[ \$should_downsample -eq 0 ]]; then
        ln -s ${bam} ${meta.id}.downsampled.bam
        ln -s ${bai} ${meta.id}.downsampled.bam.bai
    else
        samtools view -@ ${task.cpus} -b -s \${seed}.\${fraction} ${bam} > ${meta.id}.downsampled.bam
        samtools index ${meta.id}.downsampled.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed -e 's/samtools //g')
    END_VERSIONS
    """
}
