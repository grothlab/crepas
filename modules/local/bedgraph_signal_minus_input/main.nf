process BEDGRAPH_SIGNAL_MINUS_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed4799c386bc676b388fb01df7e97792ed9dd5e84c8f20ef568feca5dc2cccdb/data' :
        'community.wave.seqera.io/library/gawk:5.4.0--0877e7a42fa88325' }"

    input:
    tuple val(meta), path(signal_bdg), path(input_bdg)

    output:
    tuple val(meta), path("*.bedgraph") , emit: bedgraph
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def norm_factor = args2.contains('--norm_factor ') ? args2.split('--norm_factor ')[1].split(' ')[0] : '1'
    def min_count = args2.contains('--min_count ') ? args2.split('--min_count ')[1].split(' ')[0] : 0
    min_count = min_count.toInteger()
    """
    awk \\
        ${args} \\
        'BEGIN { FS=OFS="\\t" } NR==FNR { key=\$1 FS \$2; input_val[key]=\$4; next }
        { key=\$1 FS \$2; if (key in input_val) {
            signal=\$4; input=input_val[key];
            if (signal >= ${min_count} && input >= ${min_count}) {
                diff=signal-input;
                \$4 = (diff > ${norm_factor}) ? diff : ${norm_factor};
                print
            }
        }}' \\
    ${input_bdg} \\
    ${signal_bdg} \\
    > ${prefix}.bedgraph
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bedgraph
    """
}
