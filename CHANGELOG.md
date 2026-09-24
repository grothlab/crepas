# grothlab/crepas: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.1.0dev](https://github.com/grothlab/crepas/releases/tag/1.1.0dev)]

Development version of grothlab/crepas.

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
- `denopa` failing under Singularity and Apptainer; the pipeline now stops early and asks to remove `denopa` from `--peak_callers` or to run with Docker.
- HISAT2 alignment failing under Singularity and Apptainer, because the `hisat2/align` module's `oras://` image could no longer be pulled; the module is updated to its current nf-core version, which uses an `https://` image.
- `denopa` failing under Docker, because its container was given as `docker://hepingshiming2007/denopa:latest`, which Docker rejects.
- The software versions file for MultiQC being named `grothlabcrepas_software_mqc_versions.yml`; it is now `grothlab_crepas_software_mqc_versions.yml`.

## [[1.0.0](https://github.com/grothlab/crepas/releases/tag/1.0.0)] - Mercurian Cinnabar - 2026-06-21

Initial release of grothlab/crepas.
