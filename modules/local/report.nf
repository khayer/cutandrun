process REPORT {
    label 'process_medium'
    container 'docker://python:3.10'

    input:
    val staging_dir_path

    output:
    path "CUT_and_RUN_report.html", emit: html
    path "versions.yml", emit: versions

    script:
    """
    set -x  # Print commands for debugging
    pip install -q jinja2 pandas pyyaml > /dev/null 2>&1
    
    echo "========== REPORT PROCESS DEBUG =========="
    echo "Staging dir path: ${staging_dir_path}"
    echo "Current directory: \$(pwd)"
    echo "Python version: \$(python3 --version)"
    echo "========== Directory Contents =========="
    ls -la "${staging_dir_path}/" 2>&1 | head -30 || echo "ERROR: Cannot access staging directory!"
    echo "========== Starting Report Generation =========="
    
    python3 ${workflow.projectDir}/bin/report_generator/generator.py \\
        --staging-dir "${staging_dir_path}" \\
        --out CUT_and_RUN_report.html 2>&1
    
    echo "========== Report Generation Complete =========="
    
    cat > versions.yml <<-END_VERSIONS
        "report_generator":
            python: \$(python3 --version | awk '{print \$2}')
            jinja2: \$(python3 -c "import jinja2; print(jinja2.__version__)")
    END_VERSIONS
    """
}
