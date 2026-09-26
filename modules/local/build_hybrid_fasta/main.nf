process BUILD_HYBRID_FASTA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed4799c386bc676b388fb01df7e97792ed9dd5e84c8f20ef568feca5dc2cccdb/data' :
        'community.wave.seqera.io/library/gawk:5.4.0--0877e7a42fa88325' }"

    input:
    tuple val(meta), path(fasta, stageAs: "genome/*"), val(genome)
    tuple val(meta2), path(spikein_fasta, stageAs: "spikein_genome/*"), val(spikein_genome)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    {
        cat "${fasta}"
        awk -v suffix="_${spikein_genome}" '/^>/{print \$0 suffix; next} 1' "${spikein_fasta}"
    } > "${genome}_${spikein_genome}.fa"
    """

    stub:
    """
    touch ${genome}_${spikein_genome}.fa
    """
}
