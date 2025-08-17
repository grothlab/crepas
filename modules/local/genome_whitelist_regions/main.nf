/*
 * Prepare genome intervals for filtering by removing regions in blacklist file
 */
process GENOME_WHITELIST_REGIONS {
    tag "$sizes"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(sizes)
    path blacklist

    output:
    tuple val(meta), path("*.bed")     , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${sizes.simpleName}.whitelist"
    if (blacklist) {
        """
        sortBed ${args} -i $blacklist -g $sizes | complementBed ${args2} -i stdin -g $sizes > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        awk ${args3} '{print \$1, '0' , \$2}' OFS='\t' $sizes > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${sizes.simpleName}.whitelist"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
