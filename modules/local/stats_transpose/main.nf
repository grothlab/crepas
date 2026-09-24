process STATS_TRANSPOSE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data'
        : 'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9'}"

    input:
    tuple val(meta), file(stats)

    output:
    tuple val(meta), path("*.tsv"), emit: t_stats
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('sed'), eval("sed --version 2>&1 | head -1 | sed 's/^sed (GNU sed) //'"), emit: versions_sed, topic: versions
    tuple val("${task.process}"), val('coreutils'), eval("sort --version 2>&1 | head -1 | sed 's/^sort (GNU coreutils) //'"), emit: versions_coreutils, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${stats.baseName}.transposed"
    """
    touch ${prefix}.tsv
    """
}
