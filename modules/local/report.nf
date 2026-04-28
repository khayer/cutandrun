process REPORT {
    label 'process_medium'
    container 'python:3.10-slim'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    val staging_dir_path

    output:
    path "CUT_and_RUN_report.html", emit: html
    path "versions.yml", emit: versions

    script:
    """
    pip install -q jinja2 pandas pyyaml > /dev/null 2>&1
    
    python3 ${workflow.projectDir}/bin/report_generator/generator.py \\
        --staging-dir "${staging_dir_path}" \\
        --out CUT_and_RUN_report.html
    
    cat > versions.yml <<-END_VERSIONS
        "report_generator":
            python: \$(python3 --version | awk '{print \$2}')
            jinja2: \$(python3 -c "import jinja2; print(jinja2.__version__)")
    END_VERSIONS
    """
}
