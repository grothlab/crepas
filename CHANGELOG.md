# grothlab/crepas: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.2.0dev](https://github.com/grothlab/crepas/releases/tag/1.2.0dev)]

Development version of grothlab/crepas.

### `Added`

- LanceOtron peak caller (`--peak_callers lanceotron`, [#38](https://github.com/grothlab/crepas/issues/38)). It runs on a 1 bp RPKM coverage track of each BAM file, made with deepTools bamCoverage, with `lanceotron callPeaksInput` for IPs with an input control and `lanceotron callPeaks` for IPs without one. The peaks with an overall peak score of at least `--lanceotron_score_cutoff` (default 0.5) are kept, in a separate step, and saved as narrowPeak and as LanceOtron's table. They get FRiP scores, peak counts, HOMER annotation, peak QC plots (overall, shape and enrichment scores and p-value), MultiQC sections and an IGV track, like the other peak callers. New parameters: `--lanceotron_score_cutoff`, `--lanceotron_threshold`, `--lanceotron_window` and `--save_lanceotron_bigwig`.
- QC of the TEcount and TElocal count tables (`TE_COUNTING_QC`): the per-sample tables are merged, normalised with DESeq2 size factors and variance-stabilised, and plotted as a PCA and a sample distance heatmap (genes and TEs), also added to the MultiQC report. The merged raw count matrix, an unfiltered DESeqDataSet and the variance-stabilised object are saved for downstream analysis. Each PCA compares the IPs of one antibody with their input controls, separately per experiment type, filtering stage, downsampling method and strand. Disable with `skip_te_counting_qc`.
- epic2 peak QC plots (fold enrichment, p-value, q-value and peak length per sample), as for MACS3 and Genrich. epic2 peaks are now also saved in ENCODE broadPeak format (`*.epic2.broadPeak`), next to the HOMER-compatible `*.epic2.diffusePeak`; both are written by the new `EPIC2_TO_PEAK` module.
- Peak QC plots for SEACR (total and maximum signal, per bin size and normalisation) and MACE (p-value), so every peak caller with FRiP scores and annotation now has them. All callers use the new `PLOT_PEAK_QC` module (formerly `PLOT_MACS3_QC`), which plots the peak count, the peak length and the score columns of each caller set in `ext.args` (`--columns`); custom configs that select `PLOT_MACS3_QC` by name need updating.
- Output documentation for every peak caller (MACS3, epic2, Genrich, SEACR, MACE, EDD, DANPOS2, Consenrich and ROCCO), as subsections of "Peak calling" in `docs/output.md`.
- EDD configuration file settings as parameters: `--edd_max_ci_value`, `--edd_required_fraction_of_informative_bins`, `--edd_p_hat_ci_method` and `--edd_log_ratio_bin_size` ([#32](https://github.com/grothlab/crepas/issues/32)). The `EDD` module takes them from `ext.args` (e.g. `--max_CI_value 0.3`) and sets them in a copy of EDD's default configuration file.

### `Changed`

- The run's RO-Crate (`ro-crate-metadata_<timestamp>.json`) and the files nf-prov copies next to it (`README.md`, `main.nf`, `nextflow.config`, `nextflow_schema.json` and the input samplesheet) are saved in `pipeline_info/provenance/` instead of `pipeline_info/`.
- FRiP scores, peak counts, HOMER annotation and peak QC plots are done once for all peak callers, in the new `BAM_PEAKS_QC_ANNOTATE_HOMER` subworkflow, instead of separately in each peak caller's subworkflow. The peak caller subworkflows are renamed `BAM_PEAKS_CALL_<CALLER>` (e.g. `BAM_PEAKS_CALL_MACS3`), so custom configs that select their processes by name (`withName: '.*:BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER:...'`) need updating. Output file names and directories are unchanged.
- SEACR also runs on 1 bp coverage of the BAM files (`bedtools genomecov`), scaled by the factor of each normalisation, besides the coverage tracks at `--coverage_bin_size`. SEACR finds peaks as blocks of contiguous signal and is designed for 1 bp coverage; on 250 bp bins nearly the whole genome formed a few thousand long blocks, which made its threshold jump between very different values (rpm peak counts dropped by about 90% compared with raw for some replicates and not at all for others). The new peak files are named `*.1.<NORM>.stringent.bed`; the SEACR FRiP score, peak count and annotation files now also have the bin size in their names (e.g. `*.seacr.250.rpm.peaks.annotatePeaks.txt`), and each bin size and normalisation has its own FRiP score and peak count sections in the MultiQC report (e.g. "SEACR (1 bp, rpm) FRiP score").

### `Fixed`

- Long sample names in the consensus peak DESeq2 QC plots: the PCA labels went past the plot and the sample distance heatmap was unreadable, because `--sample_suffix` (left over from nf-core/chipseq) never matched the crepas BAM names. Samples are now named by their ID (e.g. `ChIP-seq_SPT5_T0_bRep_1`), also in the PCA and distance tables and the `size_factors/` file names.
- Intersection sizes in the consensus peak UpSet plot are now vertical.
- MACS3 peaks and consensus peaks in the IGV session pointing at missing files: the peak files were given MACS3's `*_peaks.<type>Peak` name instead of the published `*.peaks.<type>Peak`, and the consensus files lost the `/` before their antibody directory.
- SEACR peaks missing from the IGV session: every SEACR peak file (each bin size and normalisation) is now a track.
- IGV session pointing at missing files with `--multimap_allocation_method chromap`: the paths of the coverage tracks and of the MACS3, EDD, Genrich, MACE, epic2, Consenrich and ROCCO files had a `cm_allo/` directory that the files are not published in.
- Partition plots (`PARTITION_PLOT_SAMPLE`, `PARTITION_PLOT_GROUP`) missing for SCAR-seq and eSPAN since 1.1.0: the filter meant to keep OK-seq out of the plot channel (`!exp_type == 'OK-seq'`) always evaluated to false and dropped every sample; it now keeps every sample that is not OK-seq.
- `PICARD_ADDORREPLACEREADGROUPS` failing with `CIGAR_MAPS_OFF_REFERENCE` after chromap, which can report primary and secondary alignments that run past the end of a contig; it now runs with `--VALIDATION_STRINGENCY LENIENT`, like the other Picard steps.
- Institutional profiles from `--custom_config_base` (e.g. `grothlab/configs`) not being found with Nextflow 25.10, because lineage tracking was enabled; lineage is now disabled by default (enable it with `lineage.enabled = true` in a custom config).
- `Unrecognized config option` warnings for `prov.formats.wrroc.agent` and `prov.formats.wrroc.organization` with Nextflow 26.04; the pipeline no longer sets a fixed agent and organization in the run RO-Crate.
- MultiQC failing when `epic2` is in `--peak_callers`, because the gawk version reported by `EPIC2` included gawk's whole licence text and made `grothlab_crepas_software_mqc_versions.yml` invalid YAML. The `CONSENRICH` version (reported as `Version: 0.1.7b1`, which also broke the YAML) and the `EDD` version (empty, because `edd --version` prints to stderr) are fixed too.
- epic2 FRiP scores, peak counts and peak annotation missing from the MultiQC report.
- High-resolution Repli-seq chromosome heatmaps (`*.hr_repliseq.plots.pdf` and the MultiQC heatmap) closing the gaps left by bins removed when building the array, which shifted the heatmap by up to several Mb against the early/late track and the feature tracks below it. The heatmap is now drawn at the bins' own coordinates with gaps left grey, and the early/late line breaks at gaps instead of bridging them ([#163](https://github.com/grothlab/crepas/issues/163)).
- epic2 and Consenrich failing with `input file name collision` when replicates of a sample share an input control (e.g. replicates 1 and 3 both using control replicate 1). epic2 pools all controls, so each control is now passed once; Consenrich matches controls to replicates by position, so each control is staged in its own directory and repeated where it is shared.
- Genrich and Consenrich silently pairing replicates with the wrong input control when only some replicates of a sample have one; the pipeline now stops with an error naming the sample.
- epic2, Genrich and MACE FRiP scores and peak counts not being calculated, missing from the MultiQC report.
- Sample labels in the Genrich (and epic2) peak QC plots including the file extension.
- Genrich peak QC plots labelling the area under the curve (the signalValue column of Genrich peaks) as fold enrichment.
- Overlapping sample labels in the peak QC plots: the labels are file names without spaces, so they never wrapped; they now wrap after the dots.
- epic2, Genrich and MACE peak QC and peak annotation plots being written outside `downsampled/` with `--bam_downsampling_method`, and missing parts of the file name prefix (e.g. `.dedup`, `.flTbl`). They now follow the MACS3 plots. Samples in the Genrich and MACE annotation plots are now in a fixed order.
- MACE being run on input controls, as if they were samples. MACE has no input control option, so its peaks for input controls meant nothing, and on input controls with few reads it stopped the pipeline (`Only N border pairs found! MACE will be terminated`). Input controls are now left out.
- SEACR FRiP scores and peak counts missing from the MultiQC report for all but one normalisation of each sample, because MultiQC trims sample names back to the sample ID and the rows of the different normalisations got the same name. Each normalisation now has its own FRiP score and peak count sections. The SEACR sections were also missing from the MultiQC section order.
- SEACR ignoring the first bin of every bedGraph, which it skips as a header line. Its bedGraphs now start with a `track` line.
- TE counting (TEcount, TElocal and their name sorting) running once for each antibody on an input control shared by several antibodies, with identical outputs written over each other, and `TE_COUNTING_QC` then failing with `input file name collision` ([#152](https://github.com/grothlab/crepas/issues/152)). Without downsampling each input control is now counted once and its table is added to the QC of each antibody; with downsampling each antibody's copy is still counted, since it is downsampled separately.
- TEcount and TElocal tables of the final BAM files missing `.flT3` in their names for antibodies in `--use_flT2_as_total_ref` (`H3K9me3` and `H3K27me3` by default). With `--skip_flTbl`, these tables had the same names as the `pre_flT3` ones and were written over them.
- Stub runs with TE counting (`--skip_te_counting false`) failing with `input file name collision`, because the stubs of the TEcount gene and TE indices had the same name, and stub runs with `--bam_downsampling_method` failing on the empty read total of the downsampled BAM files.
- EDD stopping the pipeline with a missing output error whenever it ran successfully, because EDD writes its results into a subdirectory that the module's outputs did not look in.
- Partition/RFD summary plots drawing their grey quadrant backgrounds once per plotted data row, which left hundreds of stacked grey squares in the PDFs ([#149](https://github.com/grothlab/crepas/issues/149)).
- nf-test shards on GitHub Actions failing with `no space left on device` when pulling container images; the runner's unused preinstalled software is now removed before the tests.
- Smoothed OK-seq RFD and SCAR-seq/eSPAN partition values (in the `6_collect` tables, bedGraphs and bigWigs), and the derivatives, zero crossings and initiation zones called from them, being calculated over windows in the wrong order since 1.0.0. The windows reach `PARTITION_OR_RFD_SMOOTH` sorted by name (`chr1_1`, `chr1_10`, `chr1_100`, ...), and the smoothing averaged each window with its neighbours in that order, which mixed in windows megabases away (on mm39 OK-seq, the smoothed RFD of 12% of the windows had the wrong sign). `partition_or_rfd_smooth.pl` now orders each chromosome's windows by their number itself, and stops with an error if the forward and reverse tables list different windows. `COLLECT_PARTITIONS` now also stops with an error if its input tables, which it joins line by line, list different windows on a line; before, it only checked that they had the same number of lines. This also fixes two smaller errors in the script, given here for the default `--smooth_radius 30`: the raw value of 29 windows near the start of every chromosome being taken from another window, and wrong smoothed values on contigs of 61 to 121 windows.

### `Removed`

- `--seacr_norm`. SEACR now normalises the input control to the IP (`norm`) for raw tracks, where IP and control are not on the same scale, and compares the normalised rpm, srpm and cisrpm tracks as they are (`non`). Before, `--seacr_norm` defaulted to `non` for all tracks, so raw IP and control signal were compared without correcting for their different depths.

## [[1.1.0](https://github.com/grothlab/crepas/releases/tag/1.1.0)] - Venusian Pyrite - 2026-09-26

### `Added`

- Full `meta.yml` documentation, `nf-test` tests, and topic-based version reporting for the local modules that previously lacked them.
- Comprehensive stubs for all local modules, enabling end-to-end dry-runs of the pipeline (`-stub`) ([#61](https://github.com/grothlab/crepas/issues/61)).
- [`minibwa`](https://github.com/lh3/minibwa) as an additional read-alignment option.
- `hybrid_fasta` parameter and improved spike-in genome handling.
- Per-run provenance / RO-Crate ([WRROC](https://www.researchobject.org/workflow-run-crate/)) generation via the [`nf-prov`](https://github.com/nextflow-io/nf-prov) plugin.
- GEFION HPC configuration profile.
- Support for the `eSPAN` experiment type ([#28](https://github.com/grothlab/crepas/issues/28)). Includes a `test_espan` profile and test samplesheet.
- Grouped partition plots, combining all samples sharing an antibody into a single plot, with a `skip_partition_group_plot` parameter to disable them.
- A reduced `*.totals.tsv` summary table emitted alongside the full samtools stats summary, holding only the per-processing-step read counts.
- Experiment-type-specific alignment settings for CUT&RUN, CUT&Tag, TIP-seq and eSPAN, applied to every aligner that exposes the corresponding options (bowtie2, bowtie, HISAT2, Chromap and minimap2).
- `peak_callers`, a comma-separated list of peak callers that are all run in the same pipeline run (e.g. `--peak_callers macs3,genrich`), replacing the single-valued `peak_caller`. The requested caller is carried in each sample's `meta`, so the results of different callers no longer overwrite each other.
- `skip_peak_calling`, to switch peak calling off entirely without changing the `peak_callers` selection.
- `partition_iz_rm_overlap_range`, which sets the window used to discard overlapping initiation zones independently of `partition_plot_range`, so that the plotted range can be widened without discarding more initiation zones.

- Raw, unsmoothed partition/RFD plots drawn from the `RFD_raw_mean`/`RFD_raw_sd` columns, in both plain and SD-ribbon variants, alongside the smoothed plots.
- GPU support in the GEFION profile.
- MultiQC report sections for: per-stage SAMTools stats (library, merged, multimapper-allocated, deduplicated, filtered, blacklist-filtered and downsampled), Picard MarkDuplicates and CollectMultipleMetrics, phantompeakqualtools, MACS3 peak calling, UMI-tools extract/dedup, UMICollapse, and the deepTools plotPCA, plotCorrelation, fingerprint and gene-body/consensus-peak profile outputs.
- A parameter check that the `plotProfile` region lengths are multiples of `--coverage_bin_size`, which deepTools' `computeMatrix` requires; the mismatch previously surfaced only as a `computeMatrix` failure part-way through a run.

- Support for the `Repli-seq` experiment type, in two designs selected by a new `rt_fraction` samplesheet column (`early`, `mid`, `late`, `G1`, `S1`..`S16`); one sample may carry both. Early/late Repli-seq gives a log2(early/late) replication-timing track per sample, CPM- and optionally quantile-normalized, replicates paired or pooled, loess- or rolling-mean-smoothed, written as bedGraph and bigWig with a per-sample `*.qc.txt`, plus a replication-timing index track when three or more fractions are given. The track is segmented into domains of constant timing with DNAcopy, and genes are classified by the fraction with the highest read density over their gene body. Adds a `test_repliseq` profile and the `repliseq_*`/`save_repliseq_intermeds` parameters.

- Support for high-resolution (16-fraction) Repli-seq, following [Zhao, Sasaki & Gilbert (2020)](https://doi.org/10.1186/s13059-020-01983-8). Reads are counted per 50-kb bin and assembled into a Gaussian-smoothed, column-scaled array; a `G1` control, when present, drops the bins it does not cover (`hr_repliseq_normalization`). Initiation zones, timing transition regions, breakages, termination sites and late constant-timing regions are called from the array, each as its own BED. The published description admits more than one reading of how breakages relate to transition regions, and the readings move most of the genome between the two classes, so all three are emitted (`overlap`, `disjoint` and `flanked`), each with a partition BED assigning every analysed bin to one class. Fork speed within each transition region is estimated from its length and the stretch of S phase it spans (`hr_repliseq_s_phase_hours`).
- The smoothed, normalized coverage track of each early/late Repli-seq replicate (`*.coverage.smooth.bedGraph`, plus a bigWig), written when smoothing is applied per replicate (`--repliseq_smooth_stage replicate`, the default).
- Repli-seq outputs in the IGV session: the E/L ratio and replication-timing index bigWigs, the per-replicate smoothed coverage bigWigs, the RT domains, the per-class gene classification BEDs, and the high-resolution Repli-seq initiation zones, transition regions, breakages, termination sites, constant-timing regions and partition.
- FRiP scores, HOMER peak annotation, annotation QC plots and MultiQC sections (peak count, FRiP score and peak annotation) for SEACR peaks, as for MACS3. SEACR peaks are called once per normalization of the input bedGraph, and each set is scored and annotated separately.
- A diagram in `docs/usage.md` showing how `strandedness` relates Read 1 and Read 2 to the nascent strands.
- A pipeline website ([crepas.grothlab.workers.dev](https://crepas.grothlab.workers.dev/)). Links to the usage and output docs now point to the website.

### `Changed`

- Replaced generic `ubuntu:22.04` module containers with purpose-built Wave / BioContainers images.
- Migrated version reporting from per-module `versions.yml` files to the `versions` topic channel.
- Replaced several local modules with their nf-core equivalents (`gtf2bed` now uses ea-utils, plus `SAMTOOLS_REHEADER` and `BEDTOOLS_GENOMECOV`).
- Updated the Chromap nf-core module and reworked `fasta`/`fai` input-channel handling.
- Removed the pipeline-wide `sort_bam` parameter; alignment sorting, indexing and stats are now always handled by `BAM_SORT_STATS_SAMTOOLS`, with the aligners emitting an unsorted BAM.
- Updated numerous nf-core modules (aligners, deepTools, MACS3, HOMER, bedtools, MultiQC, FastQC, …) to their latest versions.
- Added Apptainer container definitions to the local modules.
- Updated the pipeline to the nf-core tools 4.1.0 template, bumped `nf-schema` to 2.7.3 and `nf-prov` to 1.7.0, and regenerated the parameter documentation from the schema.
- `BAM_SPLIT_BY_STRAND` now takes `exp_type` and `strandedness` as explicit inputs and selects reads with `samtools view --expr` filter expressions, excluding unmapped records from every output.
- Renamed the samtools stats summary outputs to `*.all.tsv` (every column) and `*.totals.tsv` (read counts only), so that the two tables can be matched separately.
- Renamed the multimapping-reads column of the samtools stats summary to `<column>_minus_<column>`, since the filtering step it is derived from removes more than just multimapping reads.
- The GEFION profile now prepends to the container `PATH` instead of replacing it, so images that install elsewhere (e.g. Wave/pixi images) keep working.
- Expanded `CITATIONS.md.
- SEACR is no longer restricted to CUT&RUN, CUT&Tag and TIP-seq samples; it now runs on every sample for which it is requested through `peak_callers`.
- Updated the pipeline metro map.
- The samplesheet examples in `docs/usage.md` now match the example samplesheets shipped in `assets/test-datasets/`.
- Removed the DAN System usage guide (`docs/ku_sund_danhead_crepas_usage.md`).
- Renamed `min_reps_consensus` to `consensus_min_replicates_per_sample`, and the corresponding `macs3_merged_expand.py` argument from `--min_replicates` to `--min_replicates_per_sample`.
- Removed some channel dumps written to `<outdir>/.debug/BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2/`.

- Renamed the partition plot outputs: the previous `*_plot_raw` files are now `*_plot_smoothed`, and the previous `*_plot_smoothed` files are now `*_plot_smoothed.gam`.
- `strandedness` is now optional for eSPAN samples and defaults to `forward`, with a warning when it is left unset. The column descriptions in `assets/schema_input.json` and `docs/usage.md` were unified.
- TE counting now also runs on pre-flT3 BAMs, in addition to the pre- and post-blacklist-filtering ones.
- `featureCounts` over the consensus peaks now runs separately for single- and paired-end libraries, because `-p` applies to a whole run and changes the unit counted (reads vs fragments). The quantification and DESeq2 QC outputs gain a `.se`/`.pe` suffix accordingly.
- `--gtf` is no longer a required parameter, so an annotation can be supplied through `--gff` alone, which the pipeline already converted with `gffread`.
- `--partition_iz_rm_overlap_range 0` now skips the removal of overlapping initiation zones entirely, keeping every initiation zone, instead of being treated as a zero-width window.
- The partition and RFD plots now assign every sample its own color when there are more samples than the 12-colour palette they previously used. Plots with more than 15 samples are also drawn 3 inches wider, so that a long legend does not squeeze the panel.
- The bowtie2 test profiles (`test_cutandrun`, `test_espan`, `test_hr_repliseq` and `test_repliseq`) use the prebuilt iGenomes Bowtie2 indices instead of building one on every run.
- The `test_repliseq` samplesheet has extra samples.
- The grothlab institutional profiles (`dcai_gefion`, `ku_cpromegate`, `ku_sund_danhead_mod`) moved out of the pipeline to [grothlab/configs](https://github.com/grothlab/configs); select them with `--custom_config_base https://raw.githubusercontent.com/grothlab/configs/master`.
- Local module containers are pulled over `https://` instead of `oras://`.

### `Fixed`

- flT3 orphan removal is no longer skipped, so TE counting can run on pre-blacklist-filtered BAMs.
- Chromap segmentation-fault and memory issues.
- TEtranscripts and TElocal container definitions.
- `test_cutandrun` was mapped against the wrong genome.
- Missing `versions` workflow output parameter.
- Incorrect `fasta`/chromsizes selection when both `hybrid_fasta` and `fasta` were provided.
- `strobealign` ignoring the spike-in genome index.
- Division-by-zero when a sample had no spike-in reads (now guarded and logged).
- Endogenous BAM normalization incorrectly skipped when the exogenous/spike-in BAM had no reads.
- `params.seq_platform` not being applied to the BAM read group.
- Parameter validation guards that tested `params.containsKey()` on parameters declared with a `null` default, so they were unreachable (blacklist, initiation zones / OK-seq RFD file, GTF/GFF and TE index checks). As a result, enabling blacklist filtering without a blacklist silently emptied every downstream channel while the run still reported success.
- CISRPM normalization crashing with `Cannot invoke method div() on null object` when the input control had no spike-in reads: the filter guarded the input control's endogenous total but the normalization factor used its exogenous one.
- Duplicated `--maxins` in the bowtie2 arguments, which silently overrode the CUT&RUN, CUT&Tag and TIP-seq insert-size limit with the default value.
- Bugs in `process_stats_summary.R`.
- An empty `final_samtools_stats_summary` table, and a stale `meta` reference in `STATS_CAT`.
- Add missing `skip_flTbl` param in profiles to avoid.
- `PARTITION_PLOT` failing for samples without an input control.
- Partition and RFD plot labels for sample names containing dots (e.g. `H3.3`), which were truncated at the first dot.
- The samtools stats summary for sample names containing dots (e.g. `H3.3`), which were truncated at the first dot so that several samples collapsed into one, leaving the summary table with list-columns and failing with `non-numeric argument to binary operator`.
- File name collisions in `GENRICH` when the same input control is shared between several IP samples; treatment and control BAMs are now staged into separate directories.
- Consensus peak `plotProfile` outputs of different experiment types sharing an antibody being written to the same directory; the output path and file prefix now include the experiment type.
- `featureCounts` quantification and DESeq2 QC of the consensus peaks not running.

- A bin-size mismatch between the initiation zones and the partition table is now a warning that skips only the scatter plot, instead of aborting the whole plotting step.
- MultiQC `path_filters` that matched directories.
- `CALL_PEAKS` never emitting its MultiQC channel, so FRiP scores, peak counts, HOMER annotation, consensus `featureCounts` and the DESeq2 QC plots were collected and then discarded.
- MultiQC search patterns that did not match the pipeline's own file names, so Picard MarkDuplicates metrics and the phantompeakqualtools cross-correlation scores were never parsed, and the five CollectMultipleMetrics programs were skipped.
- `plotPCA`, `plotCorrelation` and the plotFingerprint quality metrics never reaching MultiQC, and the gene-body and consensus-peak `plotProfile` outputs sharing a single section.
- Stub runs (`-stub`) stopping early or failing, because the pipeline parsed the empty stub output of some modules as real content: the khmer genome size in `PREPARE_GENOME`, the TrimGalore read count (which silently filtered out every sample), and the `BAM_FLAGSTAT_MAPPED` total. `STATS_TRANSPOSE`'s stub also named its output differently from the real run, so outputs of different stages collided in `STATS_CAT`.
- Nextflow warnings when optional reference inputs (blacklist, sparse BED, active regions, OK-seq RFD file) are not provided; their placeholder channels are now value channels.
- The nf-test suite failing on GitHub Actions: `tests/nextflow.config` now caps resources to the runner limits, and the tests exposed by that were fixed (stale module input shapes, include paths and snapshots; BAM outputs are now compared by read checksum instead of file md5).
- The MACE Docker container tag, which did not exist, and the bigtools Docker container and conda environment, which did not match the Singularity image version.
- `TRIMGALORE_HARDTRIM` output files colliding with its own input; its outputs are now prefixed `<sample>_hardtrim`.
- `denopa` is marked as unsupported for now: its only container installs DeNOPA under `/root`, which is not readable when the container runs as a regular user (Singularity, Apptainer, and Docker with the pipeline's run options). The pipeline stops early if `denopa` is requested, and its module tests are skipped.
- HISAT2 alignment failing under Singularity and Apptainer, because the `hisat2/align` module's `oras://` image could no longer be pulled; the module is updated to its current nf-core version, which uses an `https://` image.
- The software versions file for MultiQC being named `grothlabcrepas_software_mqc_versions.yml`; it is now `grothlab_crepas_software_mqc_versions.yml`.
- The MACS3 and HOMER peak QC summary tables listing samples in a random order between runs; their input files are now sorted by name.

## [[1.0.0](https://github.com/grothlab/crepas/releases/tag/1.0.0)] - Mercurian Cinnabar - 2026-06-21

Initial release of grothlab/crepas.
