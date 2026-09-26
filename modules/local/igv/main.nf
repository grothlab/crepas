process IGV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c6135cbef4e9b7044424e79d7edfc6a44bb0178174e24c3bd6e5fc1f742d1ecf/data':
        'community.wave.seqera.io/library/coreutils_python:4f8a678b09a371d6' }"

    input:
    tuple path(files), val(outpaths), val(colors)
    tuple path(fasta), val(fasta_outpath)

    output:
    path "*files.txt", emit: txt
    path "*.xml", emit: xml
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "igv_session"
    def outpath_list = outpaths.join(";")
    def color_list = colors.join(";")
    """
    # save outpaths and colors to a tab-separated file with one line per file
    paste \
        <(echo "${outpath_list}" | tr ';' '\\n') \
        <(echo "${color_list}" | tr ';' '\\n') \
        > ${prefix}.files.txt

    igv_files_to_session.py \
        ${args} \
        --file_list ${prefix}.files.txt \
        --genome_fasta ${fasta_outpath} \
        --xml_output ${prefix}.xml
    """

    stub:
    def prefix = task.ext.prefix ?: "igv_session"
    """
    touch ${prefix}.files.txt
    touch ${prefix}.xml
    """
}
