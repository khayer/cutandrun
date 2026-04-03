// Module: SUMMARIZE_PEAK_ANNOTATIONS
// Author: Katharina Hayer
// Co-created with: GitHub Copilot (GPT-5.3-Codex)
// Purpose: Summarize HOMER annotatePeaks outputs and generate replicate stacked bar plots

process SUMMARIZE_PEAK_ANNOTATIONS {
    tag "summarize_peak_annotations"
    label 'process_single'

    conda "bioconda::multiqc=1.33"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.33--pyhdfd78af_0' :
        'biocontainers/multiqc:1.33--pyhdfd78af_0' }"

    input:
    path annotate_tables

    output:
    path "homer_peak_annotation.raw_annotation_summary.tsv", emit: raw_summary
    path "homer_peak_annotation.feature_summary.tsv", emit: feature_summary
    path "homer_peak_annotation.sample_stats.tsv", emit: sample_stats
    path "homer_peak_annotation.feature_percent_table.tsv", emit: percent_table
    path "homer_peak_annotation.stacked_bar.png", emit: plot_png
    path "homer_peak_annotation.stacked_bar.pdf", emit: plot_pdf
    path "homer_peak_annotation.plot_descriptions.txt", emit: plot_descriptions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    def defaultParams = [
        normalisation_mode      : 'Spikein',
        normalisation_mode_dual : false,
        normalisation_c         : 10000,
        peakcaller              : 'macs2',
        bigwigcompare_binsize   : 50,
        skip_preseq             : false,
        publish_frip            : false,
        run_pygenometracks_top10: false
    ]

    def observedParams = [
        normalisation_mode      : params.normalisation_mode,
        normalisation_mode_dual : params.normalisation_mode_dual,
        normalisation_c         : params.normalisation_c,
        peakcaller              : params.peakcaller,
        bigwigcompare_binsize   : params.bigwigcompare_binsize,
        skip_preseq             : params.skip_preseq,
        publish_frip            : params.publish_frip,
        run_pygenometracks_top10: params.run_pygenometracks_top10
    ]

    def nonDefaultLines = observedParams
        .collect { key, value ->
            defaultParams.containsKey(key) && value != defaultParams[key] ? "- ${key}: ${value} (default: ${defaultParams[key]})" : null
        }
        .findAll { it != null }

    def nonDefaultText = nonDefaultLines ? nonDefaultLines.join('\n') : '- none'

    def preSteps = []
    preSteps << (!params.skip_fastqc ? 'FastQC was run before alignment.' : 'FastQC was skipped (--skip_fastqc).')
    preSteps << (!params.skip_trimming ? 'Read trimming was performed before alignment.' : 'Read trimming was skipped (--skip_trimming).')
    preSteps << (!params.skip_removeduplicates ? 'Duplicate removal / marking was performed.' : 'Duplicate removal was skipped (--skip_removeduplicates).')
    preSteps << (params.remove_linear_duplicates ? 'Linear duplicate removal was enabled (--remove_linear_duplicates).' : 'Linear duplicate removal was not enabled.')
    preSteps << "Peak calling used: ${params.peakcaller}."

    def preStepsText = preSteps.collect { "- ${it}" }.join('\n')

    script:
    """
    python -m venv .plotenv
    source .plotenv/bin/activate
    python -m pip install --no-cache-dir pandas matplotlib

    python ${workflow.projectDir}/bin/summarize_peak_annotations.py \
        --inputs *.annotatePeaks.txt \
        --out-prefix homer_peak_annotation

        cat <<-END_PLOT_DESCRIPTIONS > homer_peak_annotation.plot_descriptions.txt
        HOMER Peak Annotation Plot Guide
        =================================

        Short Summary
        -------------
        These plots summarize where called peaks fall across genomic features and how GC content is distributed for this dataset at the peak-annotation stage. Values reflect the selected normalization context and all upstream processing that has already been applied in this run. Use this file as a quick interpretation key and run-context record for reproducibility.

        Detailed Plot Descriptions
        --------------------------
        - homer_peak_annotation.stacked_bar.{png,pdf}
            Stacked percentage bars per replicate showing annotation feature composition (Promoter, UTR, Exon, Intron, Intergenic, etc.). Labels include peak counts and mean GC.

        - homer_peak_annotation.stacked_bar.slim_version.{png,pdf}
            Collapsed feature classes (for easier cross-sample comparison), useful when many categories are sparse.

        - homer_peak_annotation.gc_by_sample.{png,pdf}
            Boxplot of GC percent per replicate from per-peak GC values.

        - homer_peak_annotation.gc_by_sample_violin.{png,pdf}
            Violin plot of the same per-replicate GC distributions to highlight shape/density.

        - homer_peak_annotation.gc_by_condition.{png,pdf}
            Boxplot of GC percent aggregated by biological condition (replicates collapsed by removing trailing _R#).

        - homer_peak_annotation.gc_by_condition_violin.{png,pdf}
            Violin plot of GC distributions by biological condition.

        - homer_peak_annotation.*.tsv tables
            Machine-readable summaries backing the plots (raw annotation counts, feature percentages, sample stats, per-peak GC, sample-level GC, condition-level GC).

        Run Context (Normalization and Upstream Steps)
        ----------------------------------------------
        Normalization settings used in this run:
        - normalisation_mode: ${params.normalisation_mode}
        - normalisation_mode_dual: ${params.normalisation_mode_dual}
        - normalisation_c: ${params.normalisation_c}

        Upstream processing status before this plot stage:
        ${preStepsText}

        Non-default parameters detected in this run (tracked subset):
        ${nonDefaultText}

        Notes
        -----
        - This description file is generated during SUMMARIZE_PEAK_ANNOTATIONS.
        - "Condition" is inferred from sample IDs by removing trailing replicate suffixes matching _R<integer>.
        END_PLOT_DESCRIPTIONS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """
}
