process BEDGRAPH_SIGNAL_OVER_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(ip_bedgraph), path(control_bedgraph)

    output:
    tuple val(meta), path("*.bedgraph") , emit: bedgraph
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_count = args2.contains('--min_count ') ? args2.split('--min_count ')[1].split(' ')[0] : 1
    min_count = min_count.toInteger() 
    """
    awk \\
        $args \\
        'NR==FNR {key=\$1 FS \$2 FS \$3; value[key]=\$4; next} 
            {key=\$1 FS \$2 FS \$3; chip=\$4; input=(key in value ? value[key] : 0); 
            if (chip >= $min_count && input >= $min_count) { 
                ratio = chip / input; 
            } else { 
                ratio = "NaN"; 
            }
            print \$1, \$2, \$3, ratio 
            }' OFS="\\t" \\
        $control_bedgraph \\
        $ip_bedgraph \\
    > ${prefix}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
