process BEDGRAPH_SIGNAL_OVER_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed4799c386bc676b388fb01df7e97792ed9dd5e84c8f20ef568feca5dc2cccdb/data' :
        'community.wave.seqera.io/library/gawk:5.4.0--0877e7a42fa88325' }"

    input:
    tuple val(meta), path(ip_bedgraph), path(control_bedgraph)

    output:
    tuple val(meta), path("*.bedgraph") , emit: bedgraph
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_signal = args2.contains('--min_signal ') ? args2.split('--min_signal ')[1].split(' ')[0].toBigDecimal() : 0
    def min_input = args2.contains('--min_input ') ? args2.split('--min_input ')[1].split(' ')[0].toBigDecimal() : 0.001
    min_input = min_input ?: 0.001 // to prevent division by zero
    """
    awk \\
        $args \\
        'NR==FNR {key=\$1 FS \$2 FS \$3; value[key]=\$4; next}
            {key=\$1 FS \$2 FS \$3; signal=\$4; input=(key in value ? value[key] : 0);
            if (signal >= $min_signal && input >= $min_input) {
                ratio = signal / input;
            } else {
                ratio = "NaN";
            }
            print \$1, \$2, \$3, ratio
            }' OFS="\\t" \\
        $control_bedgraph \\
        $ip_bedgraph \\
    > ${prefix}.bedgraph
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bedgraph
    """
}
