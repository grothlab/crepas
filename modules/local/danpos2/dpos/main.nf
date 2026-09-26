process DANPOS2_DPOS {
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
    tuple val(meta), path("result/pooled/*.smooth.wig"), emit: pooled_treat_wig
    tuple val(meta), path("result/pooled/*.smooth.positions.xls"), emit: pooled_smooth_positions
    tuple val(meta), path("result/pooled/*_input.wig"), emit: pooled_input_wig, optional: true
    tuple val(meta), path("*.positions.integrative.xls"), emit: integrative_pos, optional: true
    tuple val(meta), path("*.reference_positions.xls"), emit: reference_pos, optional: true
    tuple val(meta), path("*.positions.ref_adjust.xls"), emit: pos_ref_adjust, optional: true
    tuple val(meta), path("result/diff/*.wig"), emit: diff_wig, optional: true
    tuple val("${task.process}"), val('danpos'), eval("danpos.py dpos -h 2>&1 | grep -oP 'danpos\\s+\\K[\\d.]+' | head -1"), emit: versions_danpos, topic: versions

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
    danpos.py dpos \\
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
    touch result/pooled/${prefix}.smooth.wig
    touch result/pooled/${prefix}.smooth.positions.xls
    touch result/pooled/${prefix}_input.wig
    touch ${prefix}.positions.integrative.xls
    touch ${prefix}.reference_positions.xls
    touch ${prefix}.positions.ref_adjust.xls
    touch result/diff/${prefix}.wig
    """
}
