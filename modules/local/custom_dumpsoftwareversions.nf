process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "bioconda::multiqc=1.33"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.33--pyhdfd78af_0' :
        'biocontainers/multiqc:1.33--pyhdfd78af_0' }"

    input:
    path versions

    output:
    path "software_versions.yml"           , emit: yml
    path "software_versions_mqc.yml"       , emit: mqc_yml
    path "software_versions_unique_mqc.yml", emit: mqc_unique_yml
    path "local_versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    import sys
    import yaml
    import platform
    from textwrap import dedent

    def _make_versions_html(versions):
        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions">
                    <thead>
                        <tr>
                            <th> Process Name </th>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]
        for process, tmp_versions in sorted(versions.items()):
            html.append("<tbody>")
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                html.append(
                    dedent(
                        f'''\\
                        <tr>
                            <td><samp>{process if (i == 0) else ''}</samp></td>
                            <td><samp>{tool}</samp></td>
                            <td><samp>{version}</samp></td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)

    def _make_versions_unique_html(versions):
        unique_versions = []

        for process, tmp_versions in sorted(versions.items()):
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                tool_version = tool + "=" + version
                if tool_version not in unique_versions:
                    unique_versions.append(tool_version)

        unique_versions.sort()

        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions-unique tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions-unique">
                    <thead>
                        <tr>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]

        for tool_version in unique_versions:
            tool_version_split = tool_version.split('=')
            html.append("<tbody>")
            html.append(
                dedent(
                    f'''\\
                    <tr>
                        <td><samp>{tool_version_split[0]}</samp></td>
                        <td><samp>{tool_version_split[1]}</samp></td>
                    </tr>
                    '''
                )
            )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)

    def _load_workflow_versions(path):
        with open(path) as f:
            raw = f.read()

        try:
            parsed = yaml.load(raw, Loader=yaml.BaseLoader)
            return parsed if isinstance(parsed, dict) else {}
        except yaml.YAMLError as e:
            # Some module version snippets can occasionally include malformed YAML.
            # Fall back to a best-effort parser so the pipeline can still finish.
            print(f"WARN: Failed to parse {path} as YAML, using best-effort parser: {e}", file=sys.stderr)

        recovered = {}
        current_process = None
        for line in raw.splitlines():
            if not line.strip() or line.lstrip().startswith('#'):
                continue

            if not line.startswith(' ') and line.rstrip().endswith(':'):
                current_process = line.strip()[:-1]
                recovered.setdefault(current_process, {})
                continue

            if current_process and line.startswith('  ') and ':' in line:
                key, value = line.strip().split(':', 1)
                recovered[current_process][key.strip()] = value.strip().strip('"').strip("'")

        return recovered

    module_versions = {}
    module_versions["${task.process}"] = {
        'python': platform.python_version(),
        'yaml': yaml.__version__
    }

    workflow_versions = _load_workflow_versions("$versions")
    workflow_versions.update(module_versions)

    workflow_versions["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version"
    }

    versions_mqc = {
        'parent_id': 'software_versions',
        'parent_name': 'Software Versions',
        'parent_description': 'Details software versions used in the pipeline run',
        'id': 'software-versions-by-process',
        'section_name': '${workflow.manifest.name} software versions by process',
        'section_href': 'https://github.com/${workflow.manifest.name}',
        'plot_type': 'html',
        'description': 'are collected at run time from the software output.',
        'data': _make_versions_html(workflow_versions)
    }

    versions_mqc_unique = {
        'parent_id': 'software_versions',
        'parent_name': 'Software Versions',
        'parent_description': 'Details software versions used in the pipeline run',
        'id': 'software-versions-unique',
        'section_name': '${workflow.manifest.name} Software Versions',
        'section_href': 'https://github.com/${workflow.manifest.name}',
        'plot_type': 'html',
        'description': 'are collected at run time from the software output.',
        'data': _make_versions_unique_html(workflow_versions)
    }

    with open("software_versions.yml", 'w') as f:
        yaml.dump(workflow_versions, f, default_flow_style=False)

    with open("software_versions_mqc.yml", 'w') as f:
        yaml.dump(versions_mqc, f, default_flow_style=False)

    with open("software_versions_unique_mqc.yml", 'w') as f:
        yaml.dump(versions_mqc_unique, f, default_flow_style=False)

    with open('local_versions.yml', 'w') as f:
        yaml.dump(module_versions, f, default_flow_style=False)
    """
}
