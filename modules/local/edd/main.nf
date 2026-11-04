process EDD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ff12cfef6551c56963bb0c98777c8eac9021d0c6150317c5d7991cd5a7007480/data'
        : 'community.wave.seqera.io/library/edd:1.1.19--56ecc27069bc9001'}"

    input:
    tuple val(meta), path(ip_bam), path(input_bam)
    tuple val(meta2), path(chrom_sizes)
    tuple val(meta3), path(unalignable_regions)

    output:
    tuple val(meta), path("*_peaks.bed"), emit: peaks
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*_bin_score.bedgraph"), optional: true, emit: bin_score
    tuple val(meta), path("*_log_ratio*.bedgraph"), optional: true, emit: log_ratios
    tuple val("${task.process}"), val('edd'), eval("edd --version 2>&1 | sed -e 's/edd //g'"), topic: versions, emit: versions_edd
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed -e 's/Python //g'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Config-file settings: move them from ext.args (e.g. `--max_CI_value 0.3`) into a copy of
    // EDD's default config file, so settings that are not given keep EDD's defaults.
    def config_keys = ['max_CI_value', 'required_fraction_of_informative_bins', 'p_hat_CI_method', 'log_ratio_bin_size']
    def config = [:]
    config_keys.each { key ->
        def match = java.util.regex.Pattern.compile("(?i)--${key}[=\\s]+(\\S+)").matcher(args)
        if (match.find()) {
            config[key] = match.group(1)
            args = match.replaceAll('').replaceAll('\\s+', ' ').trim()
        }
    }
    def config_file = "${prefix}.edd.conf"
    def write_config = config ? (
        ["cat \"\$(python -c 'import eddlib, os; print(os.path.join(os.path.dirname(eddlib.__file__), \"default_parameters.conf\"))')\" > ${config_file}"] +
        config.collect { key, value -> "sed -i 's|^${key} *=.*|${key} = ${value}|' ${config_file}" }
    ).join('\n    ') : ''
    def config_arg = config ? "--config-file ${config_file}" : ''
    // EDD names its outputs after the output directory; move them out of it
    """
    ${write_config}

    edd \\
        ${args} \\
        ${config_arg} \\
        --nprocs ${task.cpus} \\
        ${chrom_sizes} \\
        ${unalignable_regions} \\
        ${ip_bam} \\
        ${input_bam} \\
        ./${prefix}/

    mv ${prefix}/log.txt ${prefix}.log
    mv ${prefix}/* .
    rmdir ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_peaks.bed
    touch ${prefix}.log
    """
}
