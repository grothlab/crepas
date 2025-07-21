ch_bam_bai = channel.of(

    [[ id: 'exp1_IP',    control: 'exp1_INPUT', genome: 'mm10', is_control: false, exp_type: 'chipseq', antibody: 'H3K27me3', total_mapped_reads: 900000 ], 'bam1.bam', 'bam1.bai' ],
    [[ id: 'exp1_IP',    control: 'exp1_INPUT', genome: 'dm6',  is_control: false, exp_type: 'chipseq', antibody: 'H3K27me3', total_mapped_reads: 500000 ], 'bam2.bam', 'bam2.bai' ],
    [[ id: 'exp1_INPUT', control: null,         genome: 'mm10', is_control: true,  exp_type: 'chipseq', antibody: '',         total_mapped_reads: 1100000 ], 'bam3.bam', 'bam3.bai' ],

    [[ id: 'exp1_INPUT', control: null,         genome: 'dm6',  is_control: true,  exp_type: 'chipseq', antibody: '',         total_mapped_reads: 600000  ], 'bam4.bam', 'bam4.bai' ],
    [[ id: 'exp2_IP',    control: 'exp2_INPUT', genome: 'mm10', is_control: false, exp_type: 'chipseq', antibody: 'H3K27me3', total_mapped_reads: 700000 ], 'bam5.bam', 'bam5.bai' ],
    [[ id: 'exp3_IP',    control: 'exp3_INPUT', genome: 'mm10', is_control: false, exp_type: 'chipseq', antibody: 'H3K27ac',  total_mapped_reads: 600000 ], 'bam6.bam', 'bam6.bai' ],
    [[ id: 'exp3_INPUT', control: null,         genome: 'mm10', is_control: true,  exp_type: 'chipseq', antibody: '',         total_mapped_reads: 850000 ], 'bam7.bam', 'bam7.bai' ],
    [[ id: 'exp4_IP',    control: 'exp4_INPUT', genome: 'dm6',  is_control: false, exp_type: 'chipseq', antibody: 'H3K27me3', total_mapped_reads: 800000 ], 'bam8.bam', 'bam8.bai' ],
    [[ id: 'exp4_INPUT', control: null,         genome: 'dm6',  is_control: true,  exp_type: 'chipseq', antibody: '',         total_mapped_reads: 900000 ], 'bam9.bam', 'bam9.bai' ],
    [[ id: 'exp4_IP',    control: 'exp4_INPUT', genome: 'mm10', is_control: false, exp_type: 'chipseq', antibody: 'H3K27me3',  total_mapped_reads: 300000 ], 'bam10.bam', 'bam10.bai' ],
    [[ id: 'exp4_INPUT', control: null,         genome: 'mm10', is_control: true,  exp_type: 'chipseq', antibody: '',         total_mapped_reads: 400000 ], 'bam11.bam', 'bam11.bai' ]

)


        // Get the minimum number of exogenous reads among the ChIPs and inputs per experiment and antibody: [ exp_type, antibody, min_exo ]
        // First, modify the controls' metas to add their corresponding ChIP's antibody
        ch_bam_bai_genome_type.exo_control
            .combine(ch_bam_bai_genome_type.exo_ip, by: 0)
            .map { control_id, exo_control_meta, exo_control_bam, exo_control_bai, exo_ip_meta, exo_ip_bam, exo_ip_bai ->
                def meta_clone = exo_control_meta.clone()
                meta_clone.antibody = exo_ip_meta.antibody
                [ control_id, meta_clone, exo_control_bam, exo_control_bai ]
            }
            // remove duplicates based on control_meta and filename (basically the control_meta.antibody we added above)
            // Because we don't need the same input downsampled in the same way for multiple times
            .unique()
            // Now we mix the control and ip channels to evaluate them together below
            .mix( ch_bam_bai_genome_type.exo_ip )
            .map { control_id, meta, bam, bai ->
                [ meta.exp_type, meta.antibody, meta.is_control, meta.total_mapped_reads, meta, bam, bai ]
            }
            // Thus, this groups ChIPs and inputs together by experiment type, antibody and whether control or not
            .groupTuple(by: [0,1,2])
            .map { exp_type, antibody, is_control, totals, metas, bams, bais ->
                def min_exo = totals.min()
                def min_exo_id = metas[totals.indexOf(min_exo)].id
                [ exp_type, antibody, is_control, min_exo, min_exo_id, metas, bams, bais ]
            }
            // transpose back
            .transpose()
            .map { exp_type, antibody, is_control, min_exo, min_exo_id, meta, bam, bai ->
                def meta_clone = meta.clone()
                meta_clone.downsampling_prob = min_exo / meta_clone.total_mapped_reads
                meta_clone.downsampling_ref_count = min_exo
                meta_clone.downsampling_ref_id = min_exo_id
                [ meta_clone.id, meta_clone, bam, bai ]
            }
            .set { ch_bam_bai_exo }

        // Now we copy the downsampling probability by meta.id to the exogenous bams and bais:
        // We want the ratio of endo and exo reads in each sample to stay the same after downsampling
        ch_bam_bai_exo
            .join(ch_bam_bai_genome_type.endo_ip.mix(ch_bam_bai_genome_type.endo_control), by: 0)
            .map { id, exo_meta, exo_bam, exo_bai, endo_meta, endo_bam, endo_bai ->
                def meta_clone = endo_meta.clone()
                meta_clone.downsampling_prob = exo_meta.downsampling_prob
                meta_clone.downsampling_ref_count = exo_meta.downsampling_ref_count
                meta_clone.downsampling_ref_id = exo_meta.downsampling_ref_id
                [ meta_clone.id, meta_clone, endo_bam, endo_bai ]
            }
            .mix( ch_bam_bai_exo )
            // Remove the antibody from the control metas
            .map { id, meta, bam, bai ->
                def meta_clone = meta.clone()
                meta_clone.antibody = meta.is_control ? null : meta.antibody
                [ meta_clone, bam, bai ]
            }
            .set { ch_bam_bai_to_ds }

            
            ch_bam_bai_to_ds.view()
