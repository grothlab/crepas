process UMITRANSFER {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi-transfer:1.5.0--h715e4b3_0' :
        'quay.io/biocontainers/umi-transfer:1.5.0--h715e4b3_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    tuple val("${task.process}"), val('umitransfer'), eval("umi-transfer external --version 2>&1 | tail -1 | sed 's/^umi-transfer-external //'"), emit: versions_umitransfer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        umi-transfer \\
            external \\
            --in ${reads[0]} \\
            --in2 ${reads[0]} \\
            --umi ${reads[1]} \\
            --threads $task.cpus \\
            --out ${prefix}.umitransfer.fastq.gz \\
            --out2 ${prefix}.umitransfer_discard.fastq.gz \\
            --gzip \\
            $args \\
            > ${prefix}.umi_extract.log

        rm -f ${prefix}.umitransfer_discard.fastq.gz
        """
    }  else {
        """
        umi-transfer \\
            external \\
            --in ${reads[0]} \\
            --in2 ${reads[1]} \\
            --umi ${reads[2]} \\
            --threads $task.cpus \\
            --out ${prefix}.umitransfer_1.fastq.gz \\
            --out2 ${prefix}.umitransfer_2.fastq.gz \\
            --gzip \\
            $args \\
            > ${prefix}.umitransfer.log
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}.umitransfer.fastq.gz"
    } else {
        output_command = "echo '' | gzip > ${prefix}.umitransfer_1.fastq.gz ;"
        output_command += "echo '' | gzip > ${prefix}.umitransfer_2.fastq.gz"
    }
    """
    touch ${prefix}.umitransfer.log
    ${output_command}
    """
}
