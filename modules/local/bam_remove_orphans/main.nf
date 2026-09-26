
/*
 * Remove orphan reads from paired-end BAM file
 * Downloaded from: https://github.com/nf-core/chipseq/blob/76e2382b6d443db4dc2396e6831d1243256d80b0/modules/local/bam_remove_orphans.nf
 * Adapted by Samuel Ruiz-Pérez <samper@cancer.dk>.
 */
process BAM_REMOVE_ORPHANS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/0654633a56b94c91d27e119ff13755b1c605c96ccb8e392fbb6abeb0c859343d/data' :
        'community.wave.seqera.io/library/pysam_samtools:80fdc084a2b4ffc3' }"

    input:
    tuple val(meta), path(bam)
    val skip_name_sort

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('pysam'), eval("python3 -c 'import pysam; print(pysam.__version__)'"), emit: versions_pysam, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    if (!meta.single_end) {

        if (!skip_name_sort) {
            """
            if [ "\$(samtools view -c ${bam})" -eq 0 ]; then
                ln -s ${bam} ${prefix}.bam
            else
                samtools sort -n -@ ${task.cpus} -o ${prefix}.nsorted.bam -T ${prefix}.nsorted ${bam}
                bampe_rm_orphan.py ${prefix}.nsorted.bam ${prefix}.bam ${args}
            fi
            """
        } else {
            """
            if [ "\$(samtools view -c ${bam})" -eq 0 ]; then
                ln -s ${bam} ${prefix}.bam
            else
                bampe_rm_orphan.py ${bam} ${prefix}.bam ${args}
            fi
            """
        }

    } else {
        """
        ln -s ${bam} ${prefix}.bam
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
