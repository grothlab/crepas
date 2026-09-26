process DANPOS2_DREGION {
    tag "${meta.id}"
    label 'process_low_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/dd9856b8473bf4b2e67285cc76dda3f9184e9bfe8432c8b667cde9ce1f971ef3/data'
        : 'community.wave.seqera.io/library/danpos:2.2.2--34b9ae5a26a3f9b4'}"

    input:
    tuple val(meta),
          path(treatment,       stageAs: "treatment/*"),
          path(treatment_input, stageAs: "treatment_input/*"),
          val(treatment_count),
          path(control,         stageAs: 'control/*'),
          path(control_input,   stageAs: 'control_input/*'),
          val(control_count)

    output:
    tuple val(meta), path("*.regions.integrative.xls"), emit: integrative_regions, optional: true
    tuple val(meta), path("result/diff/*local_gain.refregions.xls"), emit: local_gain, optional: true
    tuple val(meta), path("result/diff/*local_loss.refregions.xls"), emit: local_loss, optional: true
    tuple val(meta), path("result/diff/*.wig"), emit: diff_wig, optional: true
    tuple val(meta), path("result/pooled/*.refregions.xls"), emit: pooled_refregions, optional: true
    tuple val(meta), path("result/pooled/*.regions.xls"), emit: pooled_regions, optional: true
    tuple val(meta), path("result/pooled/*.wig"), emit: pooled_wig, optional: true
    tuple val("${task.process}"), val('danpos'), eval("danpos.py dregion -h 2>&1 | grep -oP 'danpos\\s+\\K[\\d.]+' | head -1"), emit: versions_danpos, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_arg = control ? "treatment/:control/" : "treatment/"
    def treatment_input_arg = treatment_input ? "treatment/:treatment_input/" : ""
    def control_input_arg = control_input ? "control/:control_input/" : ""
    def input_arg = treatment_input || control_input ? "--bg ${treatment_input_arg + control_input_arg}" : ""
    def treatment_count_arg = treatment_count ? "treatment/${treatment_count}" : ""
    def control_count_arg = control_count ? "control/${control_count}" : ""
    def count_arg = treatment_count || control_count ? "--count ${treatment_count_arg + control_count_arg}" : ""
    def se = meta.single_end ? "" : "--paired 1"
    """
    danpos.py dregion \\
        ${file_arg} \\
        ${args} \\
        ${input_arg} \\
        ${count_arg} \\
        ${se} \\
        --out result/

    # replace any occurrence of 'treatment' in filenames with the prefix
    find ./ -type f -name '*treatment*' -exec bash -c 'f="{}"; mv "\$f" "\${f//treatment/${prefix}}"' \\;
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p result/pooled result/diff
    touch ${prefix}.regions.integrative.xls
    touch result/diff/${prefix}.local_gain.refregions.xls
    touch result/diff/${prefix}.local_loss.refregions.xls
    touch result/diff/${prefix}.wig
    touch result/pooled/${prefix}.refregions.xls
    touch result/pooled/${prefix}.regions.xls
    touch result/pooled/${prefix}.wig
    """
}
