# grothlab/glseq: Output

> [!IMPORTANT]
> Please read this documentation on the grothlab/glseq repository: [https://github.com/grothlab/glseq/blob/dev/docs/output.md](https://github.com/grothlab/glseq/blob/dev/docs/output.md)

## Table of Contents

1. [Introduction](#introduction)
2. [Pipeline overview](#pipeline-overview)
3. [Library-level analysis](#library-level-analysis)
    - [Raw read QC](#raw-read-qc)
    - [UMI extraction/transfer](#umi-extractiontransfer)
    - [Adapter trimming](#adapter-trimming)
    - [Alignment](#alignment)
        - [Unmapped reads](#unmapped-reads)
        - [STAR logs](#star-logs)
4. [Merged library-level analysis](#merged-library-level-analysis)
    - [Alignment merging](#alignment-merging)
    - [Preseq](#preseq)
    - [UMI-based alignment deduplication](#umi-based-alignment-deduplication)
    - [Duplicate marking](#duplicate-marking)
    - [Filtering](#filtering)
    - [Splitting alignments by genome (spike-in normalization)](#splitting-alignments-by-genome-spike-in-normalization)
    - [Allocation of multimapping reads](#allocation-of-multimapping-reads)
    - [Final filtering of BAM files](#final-filtering-of-bam-files)
    - [Collection of multiple metrics](#collection-of-multiple-metrics)
    - [Read shifting (ATAC-seq)](#read-shifting-atac-seq)
    - [phantompeakqualtools](#phantompeakqualtools)
    - [Normalized bigWig files](#normalized-bigwig-files)
    - [deepTools plots](#deeptools-plots)
    - [Peak calling](#peak-calling)
    - [Create and quantify consensus set of peaks](#create-and-quantify-consensus-set-of-peaks)
  ---

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report generated from the [full-sized test dataset](https://github.com/nf-core/test-datasets/tree/chipseq#full-test-dataset-origin) for the pipeline using a command similar to the one below:

```console
nextflow run nf-core/chipseq -profile test_full,<docker/singularity/institute>
```

The directories listed below will be created in the output directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/). See [`main README.md`](../README.md) for a condensed overview of the steps in the pipeline, and the bioinformatics tools used at each step.

See [Illumina website](https://emea.illumina.com/techniques/sequencing/dna-sequencing/chip-seq.html) for more information regarding the ChIP-seq protocol, and for an extensive list of publications.

## Library-level analysis

The initial QC and alignments are performed at the library-level e.g. if the same library has been sequenced more than once to increase sequencing depth. This has the advantage of being able to assess each library individually, and the ability to process multiple libraries from the same sample in parallel.

### Raw read QC

<details markdown="1" open>
    <summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics for read 1 (_and read2 if paired-end_) **before** adapter trimming.
- `fastqc/zips/`
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### UMI extraction/transfer

<details markdown="1" open>
    <summary>Output files</summary>

- `umitransfer/`: If UMIs are provided in a separate FastQ file, the UMI sequence will be transferred to the read name using [umi-transfer](https://github.com/SciLifeLab/umi-transfer).
  - `*umitransfer.fastq.gz`: FastQ files (single-end) containing the UMI sequence transferred to the read name.
  - `*umitransfer_1.fastq.gz`, `*umitransfer_2.fastq.gz`: FastQ files (paired-end) containing the UMI sequence transferred to the read name.

- `umitools_extract/`: If UMIs are provided in the read sequence, the UMI sequence will be extracted using [UMI-tools](https://github.com/CGATOxford/UMI-tools).
  - `*umi_extract.fastq.gz`: FastQ files (single-end) containing the UMI sequence extracted from the read.
  - `*umi_extract_1.fastq.gz`, `*umi_extract_2.fastq.gz`: FastQ files (paired-end) containing the UMI sequence extracted from the read.

</details>

Unique molecular identifiers (UMIs) are short sequences that are added to the 5' or 3' end of the read sequence. They are used to identify and remove PCR duplicates from the sequencing data. UMIs can be provided in a separate FastQ file or can be part of the read sequence in one of the main FastQ files. The pipeline supports both scenarios and will transfer ([umi-transfer](https://github.com/SciLifeLab/umi-transfer)) or extract ([UMI-tools](https://github.com/CGATOxford/UMI-tools)) the UMI sequence to the read name for further processing.

### Adapter trimming

<details markdown="1" open>
    <summary>Output files</summary>

- `trimgalore/`
  - `*fastq.gz`: If `--save_trimmed` is specified, FastQ files **after** adapter trimming will be placed in this directory.
- `trimgalore/logs/`
  - `*.log`: Log file generated by Trim Galore!.
- `trimgalore/fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics for read 1 (_and read2 if paired-end_) **after** adapter trimming.
- `trimgalore/fastqc/zips/`
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) is a wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files. By default, Trim Galore! will automatically detect and trim the appropriate adapter sequence. See [`usage.md`](usage.md) for more details about the trimming options.

![MultiQC - Cutadapt trimmed sequence plot](images/mqc_cutadapt_plot.png)

### Alignment

The pipeline has been written in a way where all the files generated downstream of the alignment are placed in the same directory as specified by `--aligner` e.g. if `--aligner bwa` is specified then all the downstream results will be placed in the `bwa/` directory. This helps with organising the directory structure and more importantly, allows the end-user to get the results from multiple aligners by simply re-running the pipeline with a different `--aligner` option along the `-resume` parameter. It also means that results won't be overwritten when resuming the pipeline and can be used for benchmarking between alignment algorithms if required. Thus, `<aligner>` in the directory structure below corresponds to the aligner set when running the pipeline.

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/library/`
  - `*.bam`: The files resulting from the alignment of individual libraries are not saved by default so this directory will not be present in your results. You can override this behaviour with the use of the `--save_align_intermeds` flag in which case it will contain the coordinate sorted alignment files in [`*.bam`](https://samtools.github.io/hts-specs/SAMv1.pdf) format.
  - `<sample>.sorted.bam.flagstat`, `<sample>.sorted.bam.idxstats` and `<sample>.sorted.bam.stats` are the SAMtools stat files generated from the alignment files.

> [!NOTE]
> File names in the resulting directory (i.e. `<aligner>/library/`) will have the '`.Lb.`' suffix.

</details>

Adapter-trimmed reads are mapped to the reference assembly using the aligner set by the `--aligner` parameter. Available aligners are [BWA](http://bio-bwa.sourceforge.net/bwa.shtml), [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Chromap](https://github.com/haowenz/chromap) (default) and [STAR](https://github.com/alexdobin/STAR). A genome index is required to run any of this aligners so if this is not provided explicitly using the corresponding parameter (e.g. `--bwa_index`), then it will be created automatically from the genome fasta input. The index creation process can take a while for larger genomes so it is possible to use the `--save_reference` parameter to save the indices for future pipeline runs, reducing processing times.

![MultiQC - SAMtools stats plot](images/mqc_samtools_stats_plot.png)

#### Unmapped reads

The `--save_unaligned` parameter enables to obtain FastQ files containing unmapped reads (only available for STAR and Bowtie2).

<details markdown="1" open>
    <summary>Output files</summary>
    
- `<aligner>/library/unmapped/`
  - `*.fastq.gz`: If `--save_unaligned` is specified, FastQ files containing unmapped reads will be placed in this directory.

</details>

#### STAR logs

<details markdown="1" open>
    <summary>Output files</summary>

- `star/library/log/`
  - `*.SJ.out.tab`: File containing filtered splice junctions detected after mapping the reads.
  - `*.Log.final.out`: STAR alignment report containing the mapping results summary.
  - `*.Log.out` and `*.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes.

</details>

## Merged library-level analysis

The library-level alignments associated with the same sample are merged and subsequently used for the downstream analyses.

### Alignment merging

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/`

  - `*.bam`: Merged library-level, coordinate sorted `*.bam` files. The file suffix for the final merged files will be `*.mLb.sorted.bam`. If you specify the `--save_align_intermeds` parameter then the unsorted merged files will be present in the directory with the suffix `*.mLb.bam`.

  - SAMtools `*.flagstat`, `*.idxstats` and `*.stats` files generated from the merged alignment files.

</details>

[Picard MergeSamFiles](https://broadinstitute.github.io/picard/command-line-overview.html) is used to merge the alignments. If you only have one library for any given replicate then the merging step is not carried out because the library-level and merged library-level BAM files will be exactly the same.



![MultiQC - Picard insert size plot](images/mqc_picard_insert_size_plot.png)



### Preseq

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/preseq/`

  - `*.lc_extrap.txt`: Preseq expected future yield file.

</details>

The [Preseq](http://smithlabresearch.org/software/preseq/) package is aimed at predicting and estimating the complexity of a genomic sequencing library, equivalent to predicting and estimating the number of redundant reads from a given sequencing depth and how many will be expected from additional sequencing using an initial sequencing experiment. The estimates can then be used to examine the utility of further sequencing, optimize the sequencing depth, or to screen multiple libraries to avoid low complexity samples. The dashed line shows a perfectly complex library where total reads = unique reads. Note that these are predictive numbers only, not absolute. The MultiQC plot can sometimes give extreme sequencing depth on the X axis - click and drag from the left side of the plot to zoom in on more realistic numbers.

![MultiQC - Preseq library complexity plot](images/mqc_preseq_plot.png)

### UMI-based alignment deduplication

a



### Duplicate marking

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/picard_markduplicates/`

  - `*.bam`: Merged library-level, coordinate sorted `*.bam` files after the marking of duplicates.

  - `*.metrics.txt`: Metrics file from MarkDuplicates.

</details>

For samples for which UMIs have not been provided, UMI-based deduplication is not possible. Thus, the pipeline will automatically use the [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html) tool to *mark* their duplicate alignments. These samples will then be specifically *filtered* for duplicates in the downstream [filtering step](#filtering) (in addition to the standard filtering criteria). The pipeline will also generate a MultiQC plot showing the percentage of duplicates in each sample.

![MultiQC - Picard deduplication stats plot](images/mqc_picard_deduplication_plot.png)


### Filtering

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/bam_filter/`
  - `*.bam`: Merged library-level, coordinate sorted `*.bam` files after filtering. The file suffix for the final filtered files will be `*.flT1.sorted.bam`. If you specify the `--save_align_intermeds` parameter then the unsorted filtered files will be present in the directory with the suffix `*.flT1.bam`.

  - `*.{bai,csi,crai}`: Index files for the filtered BAM files.

  - SAMtools `*.flagstat`, `*.idxstats` and `*.stats` files generated from the alignment files.

</details>

Alignments are then filtered using [SAMBAMBA](https://github.com/biod/sambamba) to remove:

  - Duplicates (if not already removed with UMI-based deduplication)

  - Improper pairs (in the case of paired-end samples)

  - Unmapped reads
  
### Splitting alignments by genome (spike-in normalization)

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/spikein_split/`

  - `*.bam`: Merged library-level, coordinate sorted BAM files split by genome and refiltered. The file suffix for the final filtered files will be `*.<genome>.flT2.sorted.bam` and `*.<spikein_genome>.flT2.sorted.bam`, e.g. `*.mm10.flT2.sorted.bam` and `*.dm6.flT2.sorted.bam`. If you specify the `--save_spikein_intermeds` parameter then the unsorted files will be present in the directory with the suffix `*.flT2.bam`.
  
  - `*.{bai,csi,crai}`: Index files for the split and refiltered BAM files.
  
  - SAMtools `*.flagstat`, `*.idxstats` and `*.stats` files generated from the split and refiltered files.

</details>

### Allocation of multimapping reads

Multimapping reads are reads that map to multiple locations in the genome. The `--allocation_method` parameter allows you to choose the method to use for allocating these reads. 

As with the choice of aligner, the pipeline has been written in a way where all the files generated downstream of the allocation are placed in the same directory as specified by `--allocation_method` e.g. if `--allocation_method 'allo'` is specified then all the downstream results will be placed in the `<aligner>/mergedLibrary/allo/` directory. This helps with organising the directory structure and more importantly, allows the end-user to get the results from multiple allocation methods by simply re-running the pipeline with a different `--allocation_method` option along the `-resume` parameter. It also means that results won't be overwritten when resuming the pipeline and can be used for benchmarking between allocation algorithms if required. 

Thus, `<allocation_method>` in the directory structure below corresponds to the allocation method set when running the pipeline. If multimapper allocation is disabled (by leaving the parameter `--allocate_n_multimappers 0` as it is by default) then the `--allocation_method` parameter will be ignored and the downstream directories will be placed in the `<aligner>/mergedLibrary/` directory.

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/<allocation_method>/`

  - `*.bam`: Merged library-level, coordinate sorted BAM files after the allocation of multimapping reads. The file suffix for the final filtered files will be `*.<allocation_method>.sorted.bam`. If you specify the `--save_align_intermeds` parameter then the unsorted files will be present in the directory with the suffix `*.<allocation_method>.sorted.bam`.

  - `*.{bai,csi,crai}`: Index files for the allocated BAM files.

  - SAMtools `*.flagstat`, `*.idxstats` and `*.stats` files generated from the allocated alignment files.

> [!NOTE]
> These files will not be present if the `--allocation_method` is set to `'chromap'`, since Chromap's allocation algorithm is part of the alignment step. In such case, the BAM files with allocated multimappers will be found in the `<aligner>/mergedLibrary/` directory with the suffix `*.cm_allo.sorted.bam`, which means these will go through the same merging, deduplication, filtering and spike-in splitting steps described above.

</details>

### Final filtering of BAM files

a

### Collection of multiple metrics

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/*/picard_metrics/`

  - `*_metrics`: Alignment QC files from picard CollectMultipleMetrics.

- `<aligner>/mergedLibrary/*/picard_metrics/pdf/`

  - `*.pdf`: Alignment QC plot files from picard CollectMultipleMetrics.

</details>

### Read shifting (ATAC-seq)

a

### phantompeakqualtools

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/phantompeakqualtools/`

  - `*.spp.out`, `*.spp.pdf`: phantompeakqualtools output files.

  - `*_mqc.tsv`: MultiQC custom content files.

</details>

[phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools) plots the strand cross-correlation of aligned reads for each sample. In a strand cross-correlation plot, reads are shifted in the direction of the strand they map to by an increasing number of base pairs and the Pearson correlation between the per-position read count for each strand is calculated. Two cross-correlation peaks are usually observed in a ChIP experiment, one corresponding to the read length ("phantom" peak) and one to the average fragment length of the library. The absolute and relative height of the two peaks are useful determinants of the success of a ChIP-seq experiment. A high-quality IP is characterized by a ChIP peak that is much higher than the "phantom" peak, while often very small or no such peak is seen in failed experiments.

![MultiQC - spp strand-correlation plot](images/mqc_spp_strand_correlation_plot.png)

Normalized strand coefficient (NSC) is the normalized ratio between the fragment-length cross-correlation peak and the background cross-correlation. NSC values range from a minimum of 1 to larger positive numbers. 1.1 is the critical threshold. Datasets with NSC values much less than 1.1 (< 1.05) tend to have low signal to noise or few peaks (this could be biological e.g. a factor that truly binds only a few sites in a particular tissue type OR it could be due to poor quality). ENCODE cut-off: **NSC > 1.05**.

![MultiQC - spp NSC plot](images/mqc_spp_nsc_plot.png)

Relative strand correlation (RSC) is the ratio between the fragment-length peak and the read-length peak. RSC values range from 0 to larger positive values. 1 is the critical threshold. RSC values significantly lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mis-mappings, shallow sequencing depth (significantly below saturation) or a combination of these. Like the NSC, datasets with few binding sites (< 200), which is biologically justifiable, also show low RSC scores. ENCODE cut-off: **RSC > 0.8**.

![MultiQC - spp RSC plot](images/mqc_spp_rsc_plot.png)

### Normalized bigWig files

Coverage tracks are generated for the final filtered BAM files with [deepTools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html). The coverage is calculated as the number of reads per bin ($\alpha$), where bins are short consecutive counting windows of a defined size (`--coverage_bin_size <Int>` parameter). The lengths of the reads are extended to better reflect the actual fragment length (`--coverage_extend_reads <Int>` parameter). In the case of paired-end data, reads with mates are always extended to match the fragment size defined by the two read mates, so the user-provided fragment length is only used as a fallback for singletons or mate reads that map too far apart (with a distance greater than four times the fragment length or are located on different chromosomes).

Additionally, the following normalization methods (e.g., to account for input control or spike-in reads) are available in the pipeline; note that $\alpha$ corresponds to the read count in each bin:

| Method   | Description | Formula | Output | References |
| -------- | ----------- | -------- | ------ | ---------- |
| Raw     | No normalization | $$\alpha \times 1$$ | <ul><li>$$\text{endogenous ChIP } \alpha$$</li><li>$$\text{exogenous ChIP } \alpha$$ </li><li>$$\text{endogenous input } \alpha$$</li><li>$$\text{exogenous input } \alpha$$</li></ul> | - |
| RPM     | **R**eads **P**er **M**illion mapped reads | $$\alpha \times \frac{10^6}{\text{total mapped reads}}$$ | <ul><li>$$\text{endogenous ChIP } \alpha_{\text{RPM}}$$</li><li>$$\text{exogenous ChIP } \alpha_{\text{RPM}}$$ </li><li>$$\text{endogenous input } \alpha_{\text{RPM}}$$</li><li>$$\text{exogenous input } \alpha_{\text{RPM}}$$</li></ul> | - |
| SRPM    | **S**pike-in-normalized **R**eads **P**er **M**illion mapped reads | For the endogenous ChIP: <br><br> $$\alpha \times \frac{10^6}{\text{total mapped exogenous ChIP reads}}$$ <br><br> For the endogenous input: <br><br> $$\alpha \times \frac{10^6}{\text{total mapped exogenous input reads}}$$  | <ul><li>$$\text{endogenous ChIP }\alpha_{\text{SRPM}}$$</li><li>$$\text{endogenous input }\alpha_{\text{SRPM}}$$</li></ul> | [Petryk et al. (2021)](https://doi.org/10.1038/s41596-021-00585-3) |
| CISRPM | **C**hIP-and-**I**nput-**S**pike-in-normalized **R**eads **P**er **M**illion mapped reads | For the endogenous ChIP: <br><br> $$\alpha \times \frac{10^6}{\text{total mapped exogenous ChIP reads}} \times \frac{\text{total mapped exogenous input reads}}{\text{total mapped endogenous input reads}}$$ <br><br> For the endogenous input: <br><br> $$\alpha \times \frac{10^6}{\text{total mapped exogenous input reads}} \times \frac{\text{total mapped exogenous input reads}}{\text{total mapped endogenous input reads}}$$ | <ul><li>$$\text{endogenous ChIP }\alpha_{\text{CISRPM}}$$</li><li>$$\text{endogenous input }\alpha_{\text{CISRPM}}$$</li></ul> | [Flury et al. (2023)](https://doi.org/10.1016/j.cell.2023.01.007) |
| CISRPM-SOI | **CISRPM** **S**ignal (ChIP) **O**ver **I**nput | Per bin: <br><br> <picture> <source media="(prefers-color-scheme: dark)" srcset="./images/cisrpm_soi_eq_w.png"> <source media="(prefers-color-scheme: light)" srcset="./images/cisrpm_soi_eq.png"> <img alt="Shows a black logo in light color mode and a white one in dark color mode." src="./images/cisrpm_soi_eq.png"> </picture> | $$\alpha_{\text{CISRPM-SOI}}$$ | - |

<!-- $$\alpha_{\text{CISRPM-SOI}} = \begin{cases} 0, & \text{if } \alpha_{\text{CISRPM input}} = 0 \\ \frac{\alpha_{\text{CISRPM ChIP}}}{\alpha_{\text{CISRPM input}}}, & \text{otherwise} \end{cases}$$ | $$\alpha_{\text{CISRPM-SOI}}$$  -->

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/*/<exp_type>/coverage`

  - `*.raw.<bin_size>.bigWig`: Raw bigWig files.
  - `*.rpm.<bin_size>.bigWig`: RPM bigWig files.
  - `*.srpm.<bin_size>.bigWig`: SRPM bigWig files.
  - `*.cisrpm.<bin_size>.bigWig`: CISRPM bigWig files.
  - `*.cisrpm.<bin_size>.soi.bigWig`: CISRPM-SOI bigWig files.

</details>

<!-- 
<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/*/genomecov/`
  - `*.bigWig`: Normalised bigWig files scaled to 1 million mapped reads.

</details>

The [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format is in an indexed binary format useful for displaying dense, continuous data in Genome Browsers such as the [UCSC](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](http://software.broadinstitute.org/software/igv/). This mitigates the need to load the much larger BAM files for data visualisation purposes which will be slower and result in memory issues. The coverage values represented in the bigWig file can also be normalised in order to be able to compare the coverage across multiple samples - this is not possible with BAM files. The bigWig format is also supported by various bioinformatics software for downstream processing such as meta-profile plotting. -->

### deepTools plots

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/deepTools/plotFingerprint/`

  - `*.plotFingerprint.pdf`, `*.plotFingerprint.qcmetrics.txt`, `*.plotFingerprint.raw.txt`: plotFingerprint output files.

- `<aligner>/mergedLibrary/deepTools/plotProfile/`

  - `*.computeMatrix.mat.gz`, `*.computeMatrix.vals.mat.tab`, `*.plotProfile.pdf`, `*.plotProfile.tab`, `*.plotHeatmap.pdf`, `*.plotHeatmap.mat.tab`: plotProfile output files.

</details>

[deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) plotFingerprint is a useful QC for ChIP-seq data in order to see the relative enrichment of the IP samples with respect to the controls on a genome-wide basis. The results, however, are expected to look different for example when comparing narrow marks such as transcription factors and broader marks such as histone modifications (see [plotFingerprint docs](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html)).

![MultiQC - deepTools plotFingerprint plot](images/mqc_deeptools_plotFingerprint_plot.png)

The results from deepTools plotProfile gives you a quick visualisation for the genome-wide enrichment of your samples at the TSS, and across the gene body. During the downstream analysis, you may want to refine the features/genes used to generate these plots in order to see a more specific condition-related effect.

![MultiQC - deepTools plotProfile plot](images/mqc_deeptools_plotProfile_plot.png)

### Peak calling

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/macs3/<PEAK_TYPE>/`
  - `*.xls`, `*.broadPeak` or `*.narrowPeak`, `*.gappedPeak`, `*summits.bed`: MACS3 output files - the files generated will depend on whether MACS3 has been run in _narrowPeak_ or _broadPeak_ mode.
  - `*.annotatePeaks.txt`: HOMER peak-to-gene annotation file.
- `<aligner>/mergedLibrary/macs3/<PEAK_TYPE>/qc/`
  - `macs3_peak.plots.pdf`: QC plots for MACS3 peaks.
  - `macs3_annotatePeaks.plots.pdf`: QC plots for peak-to-gene feature annotation.
  - `*.FRiP_mqc.tsv`, `*.peak_count_mqc.tsv`, `annotatepeaks.summary_mqc.tsv`: MultiQC custom-content files for FRiP score, peak count and peak-to-gene ratios.

> **NB:** `<PEAK_TYPE>` in the directory structure above corresponds to the type of peak that you have specified to call with MACS3 i.e. `broadPeak` or `narrowPeak`. If you so wish, you can call both narrow and broad peaks without redoing the preceding steps in the pipeline such as the alignment and filtering. For example, if you already have broad peaks then just add `--narrow_peak -resume` to the command you used to run the pipeline, and these will be called too! However, resuming the pipeline will only be possible if you have not deleted the `work/` directory generated by the pipeline.

</details>

[MACS3](https://github.com/macs3-project/MACS) is one of the most popular peak-calling algorithms for ChIP-seq data. By default, the peaks are called with the MACS3 `--broad` parameter. If, however, you would like to call narrow peaks then please provide the `--narrow_peak` parameter when running the pipeline. See [MACS3 outputs](https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md#output-files) for a description of the output files generated by MACS3.

![MultiQC - MACS3 total peak count plot](images/mqc_macs3_peak_count_plot.png)

[HOMER annotatePeaks.pl](http://homer.ucsd.edu/homer/ngs/annotation.html) is used to annotate the peaks relative to known genomic features. HOMER is able to use the `--gtf` annotation file which is provided to the pipeline. Please note that some of the output columns will be blank because the annotation is not provided using HOMER's in-built database format. However, the more important fields required for downstream analysis will be populated i.e. _Annotation_, _Distance to TSS_ and _Nearest Promoter ID_.

![MultiQC - HOMER annotatePeaks peak-to-gene feature ratio plot](images/mqc_annotatePeaks_feature_percentage_plot.png)

Various QC plots per sample including number of peaks, fold-change distribution, [FRiP score](https://genome.cshlp.org/content/22/9/1813.full.pdf+html) and peak-to-gene feature annotation are also generated by the pipeline. Where possible these have been integrated into the MultiQC report.

![MultiQC - MACS3 peaks FRiP score plot](images/mqc_frip_score_plot.png)

### Create and quantify consensus set of peaks

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/macs3/<PEAK_TYPE>/consensus/<ANTIBODY>/`
  - `*.bed`: Consensus peak-set across all samples in BED format.
  - `*.saf`: Consensus peak-set across all samples in SAF format. Required by featureCounts for read quantification.
  - `*.featureCounts.txt`: Read counts across all samples relative to consensus peak-set.
  - `*.annotatePeaks.txt`: HOMER peak-to-gene annotation file for consensus peaks.
  - `*.boolean.annotatePeaks.txt`: Spreadsheet representation of consensus peak-set across samples **with** gene annotation columns. The columns from individual peak files are included in this file along with the ability to filter peaks based on their presence or absence in multiple replicates/conditions.
  - `*.boolean.txt`: Spreadsheet representation of consensus peak-set across samples **without** gene annotation columns. Same as file above but without annotation columns.
  - `*.boolean.intersect.plot.pdf`, `*.boolean.intersect.txt`: [UpSetR](https://cran.r-project.org/web/packages/UpSetR/README.html) files to illustrate peak intersection.

</details>

In order to perform the differential binding analysis we need to be able to carry out the read quantification for the same intervals across **all** of the samples in the experiment. To this end, the individual peak-sets called per sample have to be merged together in order to create a consensus set of peaks.

Using the consensus peaks it is possible to assess the degree of overlap between the peaks from a set of samples e.g. _Which consensus peaks contain peaks that are common/unique to a given set of samples?_. This may be useful for downstream filtering of peaks based on whether they are called in multiple replicates/conditions. Please note that it is possible for a consensus peak to contain multiple peaks from the same sample. Unfortunately, this is sample-dependent but the files generated by the pipeline do have columns that report such instances and allow you to factor them into any further analysis.

![R - UpSetR peak intersection plot](images/r_upsetr_intersect_plot.png)

By default, the peak-sets are not filtered, therefore, the consensus peaks will be generated across the union set of peaks from all samples. However, you can increment the `--min_reps_consensus` parameter appropriately if you are confident you have good reproducibility amongst your replicates to create a "reproducible" set of consensus of peaks. In future iterations of the pipeline more formal analyses such as [IDR](https://projecteuclid.org/euclid.aoas/1318514284) may be implemented to obtain reproducible and high confidence peak-sets with which to perform this sort of analysis.

The [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) tool is used to count the number of reads relative to the consensus peak-set across all of the samples. This essentially generates a file containing a matrix where the rows represent the consensus intervals, the columns represent all of the samples in the experiment, and the values represent the raw read counts.

![MultiQC - featureCounts consensus peak read assignment plot](images/mqc_featureCounts_assignment_plot.png)

### Read counting and differential binding analysis

<details markdown="1" open>
    <summary>Output files</summary>

- `<aligner>/mergedLibrary/macs3/<PEAK_TYPE>/consensus/<ANTIBODY>/deseq2/`
  - `*.sample.dists.txt`: Spreadsheet containing sample-to-sample distance across each consensus peak.
  - `*.plots.pdf`: File containing PCA and hierarchical clustering plots.
  - `*.dds.RData`: File containing R `DESeqDataSet` object generated by DESeq2, with either
    an rlog or vst `assay` storing the variance-stabilised data.
  - `*.rds`: Alternative version of the RData file suitable for
    `readRDS` to give user control of the eventual object name.
  - `*pca.vals.txt`: Matrix of values for the first 2 principal components.
  - `R_sessionInfo.log`: File containing information about R, the OS and attached or loaded packages.
  - `<aligner>/mergedLibrary/macs3/<PEAK_TYPE>/consensus/<ANTIBODY>/sizeFactors/`
  - `*.txt`, `*.RData`: Files containing DESeq2 sizeFactors per sample.

</details>

[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) is more commonly used to perform differential expression analysis for RNA-seq datasets. However, it can also be used for ChIP-seq differential binding analysis, in which case you can imagine that instead of counts per gene for RNA-seq data we now have counts per bound region.

**This pipeline uses a standardised DESeq2 analysis script to get an idea of the reproducibility within the experiment, and to assess the overall differential binding. Please note that this will not suit every experimental design, and if there are other problems with the experiment then it may not work as well as expected.**

For larger experiments, it is recommended to use the `vst` transformation instead of the `rlog` option. This is the default behaviour and can be controlled with the `--deseq2_vst` parameter. See [DESeq2 docs](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization) for a more detailed explanation.

![MultiQC - DESeq2 PCA plot](images/mqc_deseq2_pca_plot.png)

![MultiQC - DESeq2 sample similarity plot](images/mqc_deseq2_sample_similarity_plot.png)

## Aggregate analysis

### Present QC for the raw read, alignment, peak and differential binding results

<details markdown="1" open>
    <summary>Output files</summary>

- `multiqc/<PEAK_TYPE>/`
  - `multiqc_report.html`: A standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: Directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: Directory containing static images from the report in various formats.

</details>

[MultiQC](https://multiqc.info/docs/) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available within the report data directory.

Results generated by MultiQC collate pipeline QC from FastQC, TrimGalore, samtools flagstat, samtools idxstats, samtools stats, picard CollectMultipleMetrics, picard MarkDuplicates, Preseq, deepTools plotProfile, deepTools plotFingerprint, phantompeakqualtools and featureCounts. The default [`multiqc config file`](../assets/multiqc_config.yaml) also contains the provision for loading custom-content to report peak counts, FRiP scores, peak-to-gene annnotation proportions, spp NSC coefficient, spp RSC coefficient, PCA plots and sample-similarity heatmaps.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Create IGV session file

<details markdown="1" open>
    <summary>Output files</summary>

- `igv/<PEAK_TYPE>/`
  - `igv_session.xml`: Session file that can be directly loaded into IGV.
  - `igv_files.txt`: File containing a listing of the files used to create the IGV session.

</details>

An [IGV](https://software.broadinstitute.org/software/igv/UserGuide) session file will be created at the end of the pipeline containing the normalised bigWig tracks, per-sample peaks, consensus peaks and differential sites. This avoids having to load all of the data individually into IGV for visualisation.

The genome fasta file required for the IGV session will be the same as the one that was provided to the pipeline. This will be copied into `genome/` to overcome any loading issues. If you prefer to use another path or an in-built genome provided by IGV just change the `genome` entry in the second-line of the session file.

The file paths in the IGV session file will only work if the results are kept in the same place on your storage. If the results are moved or for example, if you prefer to load the data over the web then just replace the file paths with others that are more appropriate.

Once installed, open IGV, go to `File > Open Session` and select the `igv_session.xml` file for loading.

![IGV screenshot](images/igv_screenshot.png)

## Other results

### Reference genome files

<details markdown="1" open>
    <summary>Output files</summary>

- `genome/`

  A number of genome-specific files are generated by the pipeline in order to aid in the filtering of the data, and because they are required by standard tools such as BEDTools. These can be found in this directory along with the genome fasta file which is required by IGV.

- `genome/index/`

  - `bwa/`: Directory containing BWA indices.

  - `bowtie2/`: Directory containing Bowtie2 indices.

  - `chromap/`: Directory containing Chromap indices.

  - `star/`: Directory containing STAR indices.

  - If the `--save_reference` parameter is provided then the alignment indices generated by the pipeline will be saved in this directory. This can be quite a time-consuming process so it permits their reuse for future runs of the pipeline or for other purposes.

</details>

Reference genome-specific files can be useful to keep for the downstream processing of the results.

### Pipeline information

<details markdown="1" open>
    <summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to trouble-shoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
