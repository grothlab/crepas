process SAMTOOLS_STATS_TRANSPOSE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'nf-core/ubuntu:22.04'}"

    input:
    tuple val(meta), file(stats)

    output:
    tuple val(meta), path("*.tsv"), emit: t_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${stats.baseName}.transposed"
    """
    basename=\$(basename ${stats} .stats)

    grep ${args} ^SN ${stats} | cut -f 2-3 | sed 's/:\\t/\\t/' > ${prefix}.tmp
    
    # transpose
    awk ${args2} '
    BEGIN { FS=OFS="\\t" }
    {
        for (rowNr=1;rowNr<=NF;rowNr++) {
            cell[rowNr,NR] = \$rowNr
        }
        maxRows = (NF > maxRows ? NF : maxRows)
        maxCols = NR
    }
    END {
        for (rowNr=1;rowNr<=maxRows;rowNr++) {
            for (colNr=1;colNr<=maxCols;colNr++) {
                printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
            }
        }
    }' ${prefix}.tmp > ${prefix}.tmp2

    # prepend column with ID as first row and the basename as second row
    awk ${args3} -v basename="\$basename" 'NR==1{print "ID\\t"\$0} NR==2{print basename"\\t"\$0} NR>2{print}' ${prefix}.tmp2 > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
