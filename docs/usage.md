# grothlab/glseq: Usage

> [!IMPORTANT]
> Please read this documentation on the grothlab/glseq repository: [https://github.com/grothlab/glseq/blob/dev/docs/usage.md](https://github.com/grothlab/glseq/blob/dev/docs/usage.md)

## Table of Contents

1. [Samplesheet input](#samplesheet-input)
2. [Reference genome files](#reference-genome-files)
3. [Running the pipeline](#running-the-pipeline)
    1. [Parameters](#parameters)
    2. [Updating the pipeline](#updating-the-pipeline)
    3. [Reproducibility](#reproducibility)
    4. [Core Nextflow arguments](#core-nextflow-arguments)
    5. [Custom configuration](#custom-configuration)
    6. [Running in the background](#running-in-the-background)
    7. [Nextflow memory requirements](#nextflow-memory-requirements)
---

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location. It has to be a comma-separated (`.csv`) file with with 11 columns and a header row as explained below.

```bash
--input <path_to_samplesheet_csv>
```

| Column   | Description |
| -------- | ----------- |
| `sample` |  Custom sample name. This identifier should be identical when you have multiple replicates from the same experimental group; just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be `1`. Avoid including the experiment type in this identifier, since it will be parsed from the `exp_type` column and prepended to the sample name by default. |
| `fastq_1` | **Full path** to FASTQ file for reads 1. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`. |
| `fastq_2` | **Full path** to FASTQ file for reads 2. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`. Leave empty for single-end data. |
| `fastq_umi` | **Full path** to the corresponding UMI FASTQ file for deduplication. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`. Leave empty if a separate UMI file is not available. |
| `okseq_part_file` |**Full path** to the corresponding OK-seq partition file. Leave empty if OK-seq data is not available. Only for SCAR-seq data. |
| `replicate` | Integer representing replicate number. This will be identical for re-sequenced libraries (technical replicates). Must start from `1..<number of replicates>`. |
| `exp_type` | One of `chipseq`, `atacseq`, `scarseq`, `chorseq`, `ChIP-exo`, or `OK-seq` |
| `strandedness` | Either `forward` or `reverse`: <br><br><ul><li>If the library was prepared using NGS indexed PentAdapter™ adapters (PentaBase ApS, Denmark), as in the [SCAR-seq paper](https://doi.org/10.1038/s41596-021-00585-3), set this to `forward`.</li><li>If the library was prepared using [xGen™ UDI-UMI Adapters](https://eu.idtdna.com/page/products/next-generation-sequencing/ngs-adapters-indexing-primers) (Integrated DNA Technologies, Inc.), the [insert strandedness is flipped](https://eu.idtdna.com/pages/support/faqs/can-the-xgen-unique-dual-index-umi-adapters-be-used-for-rna-seq), so set this to `reverse`.</li></ul><br>This field is only relevant for SCAR-seq and OK-seq; leave it empty for unstranded data (ChIP-seq, ChIP-exo, ChOR-seq, or ATAC-seq). |
| `antibody` | This column is required to separate the downstream consensus peak merging for different antibodies. It is not advisable to generate a consensus peak set across different antibodies especially if their binding patterns are inherently different e.g. narrow transcription factors and broad histone marks. It should be empty in the case of input control rows. It is required when the `control` field is specified. |
| `control` | This column should contain the `sample` identifier of the corresponding input control for the IP. It should be empty in the case of input control rows. It is required when the `antibody` field is specified. |
| `control_replicate` | Integer representing the replicate number for the corresponding input control sample. It should be empty in the case of input control rows. |


### Example 1: Multiple biological replicates

This is an example of a samplesheet for a ChIP-seq experiment with one condition and two biological replicates for each antibody:

| sample | fastq_1 | fastq_2 | fastq_umi | okseq_part_file | replicate | exp_type | strandedness | antibody | control | control_replicate |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| condition_1_H3K9me3 | condition_1_bRep1_H3K9me3_R1.fastq.gz | condition_1_bRep1_H3K9me3_R3.fastq.gz | condition_1_bRep1_H3K9me3_R2.fastq.gz | | 1 | chipseq | | H3K9me3 | condition_1_INPUT | 1 |
| condition_1_H3K9me3 | condition_1_bRep2_H3K9me3_R1.fastq.gz | condition_1_bRep2_H3K9me3_R3.fastq.gz | condition_1_bRep2_H3K9me3_R2.fastq.gz | | 2 | chipseq | | H3K9me3 | condition_1_INPUT | 2 |
| condition_1_H3K27ac | condition_1_bRep1_H3K27ac_R1.fastq.gz | condition_1_bRep1_H3K27ac_R3.fastq.gz | condition_1_bRep1_H3K27ac_R2.fastq.gz | | 1 | chipseq | | H3K27ac | condition_1_INPUT | 1 |
| condition_1_H3K27ac | condition_1_bRep2_H3K27ac_R1.fastq.gz | condition_1_bRep2_H3K27ac_R3.fastq.gz | condition_1_bRep2_H3K27ac_R2.fastq.gz | | 2 | chipseq | | H3K27ac | condition_1_INPUT | 2 |
| condition_1_INPUT | condition_1_bRep1_INPUT_R1.fastq.gz | condition_1_bRep1_INPUT_R3.fastq.gz | condition_1_bRep1_INPUT_R2.fastq.gz | | 1 | chipseq | | | | |
condition_1_INPUT | condition_1_bRep2_INPUT_R1.fastq.gz | condition_1_bRep2_INPUT_R3.fastq.gz | condition_1_bRep2_INPUT_R2.fastq.gz | | 2 | chipseq | | | | |

> [!NOTE]
> You can download this example samplesheet [here](../assets/samplesheets/ex1_multiBioRep_samplesheet.csv) or copy and save the cell below:

```csv
sample,fastq_1,fastq_2,fastq_umi,okseq_part_file,replicate,exp_type,strandedness,antibody,control,control_replicate
condition_1_H3K9me3,condition_1_bRep1_H3K9me3_R1.fastq.gz,condition_1_bRep1_H3K9me3_R3.fastq.gz,condition_1_bRep1_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,condition_1_INPUT,1
condition_1_H3K9me3,condition_1_bRep2_H3K9me3_R1.fastq.gz,condition_1_bRep2_H3K9me3_R3.fastq.gz,condition_1_bRep2_H3K9me3_R2.fastq.gz,,2,chipseq,,H3K9me3,condition_1_INPUT,2
condition_1_H3K27ac,condition_1_bRep1_H3K27ac_R1.fastq.gz,condition_1_bRep1_H3K27ac_R3.fastq.gz,condition_1_bRep1_H3K27ac_R2.fastq.gz,,1,chipseq,,H3K27ac,condition_1_INPUT,1
condition_1_H3K27ac,condition_1_bRep2_H3K27ac_R1.fastq.gz,condition_1_bRep2_H3K27ac_R3.fastq.gz,condition_1_bRep2_H3K27ac_R2.fastq.gz,,2,chipseq,,H3K27ac,condition_1_INPUT,2
condition_1_INPUT,condition_1_bRep1_INPUT_R1.fastq.gz,condition_1_bRep1_INPUT_R3.fastq.gz,condition_1_bRep1_INPUT_R2.fastq.gz,,1,chipseq,,,,
condition_1_INPUT,condition_1_bRep2_INPUT_R1.fastq.gz,condition_1_bRep2_INPUT_R3.fastq.gz,condition_1_bRep2_INPUT_R2.fastq.gz,,2,chipseq,,,,
```

### Example 2: Multiple runs of the same library (technical replicates)

Both the `sample` and `replicate` identifiers have to be the same when you have sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will perform the alignments in parallel, and subsequently merge them before further analysis.

This is an example of a samplesheet for a ChIP-seq experiment with one condition, two biological replicates for each antibody, and two technical replicates for each biological replicate:

| sample | fastq_1 | fastq_2 | fastq_umi | okseq_part_file | replicate | exp_type | strandedness | antibody | control | control_replicate |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| condition_1_H3K9me3 | condition_1_bRep1_tRep1_H3K9me3_R1.fastq.gz | condition_1_bRep1_tRep1_H3K9me3_R3.fastq.gz | condition_1_bRep1_tRep1_H3K9me3_R2.fastq.gz | | 1 | chipseq | | H3K9me3 | condition_1_INPUT | 1 |
| condition_1_H3K9me3 | condition_1_bRep1_tRep2_H3K9me3_R1.fastq.gz | condition_1_bRep1_tRep2_H3K9me3_R3.fastq.gz | condition_1_bRep1_tRep2_H3K9me3_R2.fastq.gz | | 1 | chipseq | | H3K9me3 | condition_1_INPUT | 1 |
| condition_1_H3K9me3 | condition_1_bRep2_tRep1_H3K9me3_R1.fastq.gz | condition_1_bRep2_tRep1_H3K9me3_R3.fastq.gz | condition_1_bRep2_tRep1_H3K9me3_R2.fastq.gz | | 2 | chipseq | | H3K9me3 | condition_1_INPUT | 2 |
| condition_1_H3K9me3 | condition_1_bRep2_tRep2_H3K9me3_R1.fastq.gz | condition_1_bRep2_tRep2_H3K9me3_R3.fastq.gz | condition_1_bRep2_tRep2_H3K9me3_R2.fastq.gz | | 2 | chipseq | | H3K9me3 | condition_1_INPUT | 2 |
| condition_1_H3K27ac | condition_1_bRep1_tRep1_H3K27ac_R1.fastq.gz | condition_1_bRep1_tRep1_H3K27ac_R3.fastq.gz | condition_1_bRep1_tRep1_H3K27ac_R2.fastq.gz | | 1 | chipseq | | H3K27ac | condition_1_INPUT | 1 |
| condition_1_H3K27ac | condition_1_bRep1_tRep2_H3K27ac_R1.fastq.gz | condition_1_bRep1_tRep2_H3K27ac_R3.fastq.gz | condition_1_bRep1_tRep2_H3K27ac_R2.fastq.gz | | 1 | chipseq | | H3K27ac | condition_1_INPUT | 1 |
| condition_1_H3K27ac | condition_1_bRep2_tRep1_H3K27ac_R1.fastq.gz | condition_1_bRep2_tRep1_H3K27ac_R3.fastq.gz | condition_1_bRep2_tRep1_H3K27ac_R2.fastq.gz | | 2 | chipseq | | H3K27ac | condition_1_INPUT | 2 |
| condition_1_H3K27ac | condition_1_bRep2_tRep2_H3K27ac_R1.fastq.gz | condition_1_bRep2_tRep2_H3K27ac_R3.fastq.gz | condition_1_bRep2_tRep2_H3K27ac_R2.fastq.gz | | 2 | chipseq | | H3K27ac | condition_1_INPUT | 2 |
| condition_1_INPUT | condition_1_bRep1_tRep1_INPUT_R1.fastq.gz | condition_1_bRep1_tRep1_INPUT_R3.fastq.gz | condition_1_bRep1_tRep1_INPUT_R2.fastq.gz | | 1 | chipseq | | | | |
| condition_1_INPUT | condition_1_bRep1_tRep2_INPUT_R1.fastq.gz | condition_1_bRep1_tRep2_INPUT_R3.fastq.gz | condition_1_bRep1_tRep2_INPUT_R2.fastq.gz | | 1 | chipseq | | | | |
| condition_1_INPUT | condition_1_bRep2_tRep1_INPUT_R1.fastq.gz | condition_1_bRep2_tRep1_INPUT_R3.fastq.gz | condition_1_bRep2_tRep1_INPUT_R2.fastq.gz | | 2 | chipseq | | | | |
| condition_1_INPUT | condition_1_bRep2_tRep2_INPUT_R1.fastq.gz | condition_1_bRep2_tRep2_INPUT_R3.fastq.gz | condition_1_bRep2_tRep2_INPUT_R2.fastq.gz | | 2 | chipseq | | | | |

> [!NOTE]
> You can download this example samplesheet [here](../assets/samplesheets/ex2_multiTechRep_samplesheet.csv) or copy and save the cell below:

```csv
sample,fastq_1,fastq_2,fastq_umi,okseq_part_file,replicate,exp_type,strandedness,antibody,control,control_replicate
condition_1_H3K9me3,condition_1_bRep1_tRep1_H3K9me3_R1.fastq.gz,condition_1_bRep1_tRep1_H3K9me3_R3.fastq.gz,condition_1_bRep1_tRep1_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,condition_1_INPUT,1
condition_1_H3K9me3,condition_1_bRep1_tRep2_H3K9me3_R1.fastq.gz,condition_1_bRep1_tRep2_H3K9me3_R3.fastq.gz,condition_1_bRep1_tRep2_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,condition_1_INPUT,1
condition_1_H3K9me3,condition_1_bRep2_tRep1_H3K9me3_R1.fastq.gz,condition_1_bRep2_tRep1_H3K9me3_R3.fastq.gz,condition_1_bRep2_tRep1_H3K9me3_R2.fastq.gz,,2,chipseq,,H3K9me3,condition_1_INPUT,2
condition_1_H3K9me3,condition_1_bRep2_tRep2_H3K9me3_R1.fastq.gz,condition_1_bRep2_tRep2_H3K9me3_R3.fastq.gz,condition_1_bRep2_tRep2_H3K9me3_R2.fastq.gz,,2,chipseq,,H3K9me3,condition_1_INPUT,2
condition_1_H3K27ac,condition_1_bRep1_tRep1_H3K27ac_R1.fastq.gz,condition_1_bRep1_tRep1_H3K27ac_R3.fastq.gz,condition_1_bRep1_tRep1_H3K27ac_R2.fastq.gz,,1,chipseq,,H3K27ac,condition_1_INPUT,1
condition_1_H3K27ac,condition_1_bRep1_tRep2_H3K27ac_R1.fastq.gz,condition_1_bRep1_tRep2_H3K27ac_R3.fastq.gz,condition_1_bRep1_tRep2_H3K27ac_R2.fastq.gz,,1,chipseq,,H3K27ac,condition_1_INPUT,1
condition_1_H3K27ac,condition_1_bRep2_tRep1_H3K27ac_R1.fastq.gz,condition_1_bRep2_tRep1_H3K27ac_R3.fastq.gz,condition_1_bRep2_tRep1_H3K27ac_R2.fastq.gz,,2,chipseq,,H3K27ac,condition_1_INPUT,2
condition_1_H3K27ac,condition_1_bRep2_tRep2_H3K27ac_R1.fastq.gz,condition_1_bRep2_tRep2_H3K27ac_R3.fastq.gz,condition_1_bRep2_tRep2_H3K27ac_R2.fastq.gz,,2,chipseq,,H3K27ac,condition_1_INPUT,2
condition_1_INPUT,condition_1_bRep1_tRep1_INPUT_R1.fastq.gz,condition_1_bRep1_tRep1_INPUT_R3.fastq.gz,condition_1_bRep1_tRep1_INPUT_R2.fastq.gz,,1,chipseq,,,,
condition_1_INPUT,condition_1_bRep1_tRep2_INPUT_R1.fastq.gz,condition_1_bRep1_tRep2_INPUT_R3.fastq.gz,condition_1_bRep1_tRep2_INPUT_R2.fastq.gz,,1,chipseq,,,,
condition_1_INPUT,condition_1_bRep2_tRep1_INPUT_R1.fastq.gz,condition_1_bRep2_tRep1_INPUT_R3.fastq.gz,condition_1_bRep2_tRep1_INPUT_R2.fastq.gz,,2,chipseq,,,,
condition_1_INPUT,condition_1_bRep2_tRep2_INPUT_R1.fastq.gz,condition_1_bRep2_tRep2_INPUT_R3.fastq.gz,condition_1_bRep2_tRep2_INPUT_R2.fastq.gz,,2,chipseq,,,,
```
<!-- 
### Example 3: Full design

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet.

A final design file may look something like the one below. This is for two antibodies and associated controls, where the second replicate of the `WT_BCATENIN_IP` and `NAIVE_BCATENIN_IP` samples have been sequenced twice:

> [!NOTE]
> You can download this example samplesheet [here]() or copy and save the cell below: -->

## Reference genome files

### Genome index, FASTA, GTF/GFF and gene BED files

The minimum reference genome requirements are a FASTA and a GTF file, all other files required to run the pipeline can be generated from these files. However, it is more storage and compute friendly if you are able to re-use reference genome files as efficiently as possible. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices (e.g. those unavailable on [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space. You can then either provide the appropriate reference genome files on the command-line via the appropriate parameters (e.g. `--bwa_index '/path/to/bwa/index/'`) or via a [custom config file](https://nf-co.re/usage/configuration#custom-configuration-files).

<!-- - If `--genome` is provided then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`. -->
- If `--gene_bed` is not provided then it will be generated from the GTF file.

> **NB:** Compressed reference files are also supported by the pipeline i.e. standard files with the `.gz` extension and indices folders with the `tar.gz` extension.

### Genome blacklist regions

The blacklist bed files where obtained using the commands below:

```console
cd ..
mkdir -p v1.0
cd v1.0
wget -L https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz && gunzip ENCFF001TDO.bed.gz && mv ENCFF001TDO.bed hg19-blacklist.v1.bed

mkdir -p assets/blacklists/v2.0/
cd assets/blacklists/v2.0/
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz && gunzip ce10-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz && gunzip ce11-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz && gunzip dm3-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz && gunzip dm6-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz && gunzip hg19-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz && gunzip hg38-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz && gunzip mm10-blacklist.v2.bed.gz

cd ..
mkdir -p v3.0
cd v3.0
wget -L https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz && gunzip ENCFF356LFX.bed.gz && mv ENCFF356LFX.bed hg38-blacklist.v3.bed
```

> **NB:** A detailed description of the different versions of the files can be found [here](https://sites.google.com/site/anshulkundaje/projects/blacklists). Also, to to see which blacklist bed files are assigned by default to the respective reference genome check the [igenomes.config](conf/igenomes.config).


### Initiation zones
-

### OK-seq partitions
-

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run grothlab/glseq \
  -profile docker \
  --input samplesheet.csv \
  --genome GRCh37 \
  --outdir <OUTDIR>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
> The above pipeline run specified with a params file in yaml format:

```bash
nextflow run grothlab/glseq -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
input: 'data'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Parameters

#### Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row. See [usage docs](https://github.com/grothlab/glseq/blob/main/docs/usage.md).</small></details>| `string` |  | true |  |
| `fragment_size` | Estimated fragment size used to extend single-end reads. | `integer` | 150 |  |  |
| `seq_platform` | Platform/technology used to produce the reads. Corresponds to the `PL` tag in the SAM/BAM file header. <details><summary>Help</summary><small>See the [SAM format specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf). Valid values: CAPILLARY, DNBSEQ (MGI/BGI), ELEMENT, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT (Oxford Nanopore), PACBIO (Pacific Biosciences), SINGULAR, SOLID, and ULTIMA. This field should be omitted when the technology is not in this list (though the PM field may still be present in this case) or is unknown.</small></details>| `string` |  |  |  |
| `seq_center` | Sequencing center information to be added to read group of BAM files. | `string` |  |  |  |
| `read_length` | Read length used to calculate MACS3 genome size for peak calling if `--macs_gsize` isn't provided. | `integer` | 50 |  |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | true |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

#### Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genome` | Name of reference genome (e.g., `GRCh38`, `hg19`, `mm10`). | `string` |  | true |  |
| `fasta` | Path to FASTA genome file. <details><summary>Help</summary><small>If you don't have the appropriate alignment index available this will be generated for you automatically based on this file. Combine with `--save_reference` to save alignment index for future runs.</small></details>| `string` |  | true |  |
| `gtf` | Path to GTF annotation file. | `string` |  | true |  |
| `gff` | Path to GFF3 annotation file. | `string` |  |  |  |
| `bwa_index` | Path to directory or tar.gz archive for pre-built BWA index. | `string` |  |  |  |
| `bowtie2_index` | Path to directory or tar.gz archive for pre-built Bowtie2 index. | `string` |  |  |  |
| `chromap_index` | Path to directory or tar.gz archive for pre-built Chromap index. | `string` |  |  |  |
| `star_index` | Path to directory or tar.gz archive for pre-built STAR index. | `string` |  |  |  |
| `hisat2_index` | Path to directory or tar.gz archive for pre-built HISAT2 index. | `string` |  |  |  |
| `gene_bed` | Path to BED file containing gene intervals. This will be created from the GTF file if not specified. | `string` |  |  |  |
| `initiation_zones` | Path to BED file containing initiation zones for SCAR-seq analysis. | `string` |  |  |  |
| `macs_gsize` | Effective genome size parameter required by MACS3. <details><summary>Help</summary><small>[Effective genome size](https://github.com/taoliu/MACS#-g--gsize) parameter required by MACS3.</small></details>| `number` |  |  |  |
| `blacklist` | Path to blacklist regions in BED format, used for filtering alignments. <details><summary>Help</summary><small>If provided, alignments that overlap with the regions in this file will be filtered out (see [ENCODE blacklists](https://sites.google.com/site/anshulkundaje/projects/blacklists)). The file should be in BED format.</small></details>| `string` |  |  |  |
| `sparsebed` | Path to BED file with sparse segments for Consenrich. <details><summary>Help</summary><small>Defines the annotation file of 'sparse segments' (genomic regions) from the complement of an 'inflated' BED file of known active regions. See https://github.com/nolan-h-hamilton/Consenrich/blob/5e7f287d81ac4a9bc098ae29e825a5ded07af3b8/consenrich/consenrich.py#L58C54-L58C174</small></details>| `string` |  |  |  |
| `active_regions` | Path to BED file with active regions for Consenrich. <details><summary>Help</summary><small>Defines the annotation file of 'active regions' (genomic regions) for Consenrich. See https://github.com/nolan-h-hamilton/Consenrich/blob/5e7f287d81ac4a9bc098ae29e825a5ded07af3b8/consenrich/consenrich.py#L67</small></details>| `string` |  |  |  |
| `rocco_params` | Path to ROCCO parameters file. <details><summary>Help</summary><small>CSV file containing chromosome-specific optimization parameters. Each supported genome has a custom `--params` file packaged with ROCCO, but a custom file has to be manually passed to ROCCO within this pipeline. Consider using this argument to set chromosome-specific budget and gamma values custom to your data and preferences. See the [ROCCO documentation](https://nolan-h-hamilton.github.io/ROCCO/#rocco-constants) for more details.</small></details>| `string` |  |  |  |
| `splicesites` | Path to splice sites file for HISAT2. <details><summary>Help</summary><small>If using HISAT2 as the aligner, you can provide a list of known splice sites, which HISAT2 makes use of to align reads with small anchors. If this file is not provided, it will be generated from the GTF file.</small></details>| `string` |  |  |  |
| `save_reference` | If generated by the pipeline save the BWA index in the results directory. <details><summary>Help</summary><small>If the BWA index is generated by the pipeline use this parameter to save it to your results folder. These can then be used for future pipeline runs, reducing processing times.</small></details>| `boolean` | true |  |  |
| `igenomes_base` | Directory / URL base for iGenomes references. | `string` | s3://ngi-igenomes/igenomes/ |  | true |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` | true |  | true |

#### UMI extraction and deduplication options

Options to adjust UMI extraction and deduplication criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_umi` | Set this parameter to enable UMI extraction and deduplication. | `boolean` | false |  |  |
| `skip_umi_extract` | Skip UMI extraction step. <details><summary>Help</summary><small>Use this flag to skip the UMI extraction/transfer step. This is useful if you have already extracted UMIs from your data but want to run the deduplication</small></details>| `boolean` | true |  |  |
| `umitools_extract_method` | UMI pattern to use. Can be either `string` (default) or `regex`. <details><summary>Help</summary><small>More details can be found in the [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#extract-method).<br></small></details>| `string` | string |  |  |
| `umitools_bc_pattern` | The UMI barcode pattern to use e.g. 'NNNNNN' indicates that the first 6 nucleotides of the read are from the UMI. <details><summary>Help</summary><small>More details can be found in the [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#extract-method).</small></details>| `string` | NNNNNNNNN |  |  |
| `umitools_bc_pattern2` | The UMI barcode pattern to use if the UMI is located in read 2. | `string` | NNNNNNNNN |  |  |
| `umi_separator` | Separator used to separate the UMI from the read name. | `string` | _ |  |  |
| `umi_grouping_method` | Method to use to determine read groups by subsuming those with similar UMIs. All methods start by identifying the reads with the same mapping position, but treat similar yet nonidentical UMIs differently. | `string` | percentile |  |  |
| `umi_discard_read` | After UMI barcode extraction discard either R1 or R2 by setting this parameter to `1` or `2`, respectively. | `integer` |  |  |  |
| `get_dedup_stats` | Get deduplication statistics. | `boolean` | true |  |  |
| `umi_dedup_tool` | UMI deduplication tool to use. | `string` | umitools |  |  |
| `save_dedup_intermeds` | Save intermediate files from the UMI deduplication step. | `boolean` |  |  |  |

#### Adapter trimming options

Options to adjust adapter trimming criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `clip_r1` | Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads). | `integer` |  |  |  |
| `clip_r2` | Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only). | `integer` |  |  |  |
| `three_prime_clip_r1` | Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed. | `integer` |  |  |  |
| `three_prime_clip_r2` | Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed. | `integer` |  |  |  |
| `trim_full_length` | Instructs Trim Galore to search for the full length (33 bp) TruSeq (ChIP-seq, ChOR-seq, SCAR-seq) or Nextera (ATAC-seq) adapter sequences for both complement and reverse complements. <details><summary>Help</summary><small>The default (false) setting will mean that the auto-detect function of Trim Galore will be used. See [https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#adapter-auto-detection](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#adapter-auto-detection)</small></details>| `boolean` | true |  |  |
| `trim_nextseq` | Instructs Trim Galore to apply the `--nextseq=X` option, to trim based on quality after removing poly-G tails. <details><summary>Help</summary><small>This enables the Cutadapt `--nextseq-trim=3'CUTOFF` option via Trim Galore, which will set a quality cutoff (that is normally given with -q instead), but qualities of G bases are ignored. This trimming is in common for the NextSeq- and NovaSeq-platforms, where basecalls without any signal are called as high-quality G bases. Read the [Cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#nextseq-trim) for more details.</small></details>| `integer` | 20 |  |  |
| `trim_stringency` | Instructs Trim Galore to use this stringency (minimum overlap) when searching for adapter sequences. <details><summary>Help</summary><small>In other words, to reduce the number of falsely trimmed bases, the alignment algorithm in Trim Galore requires that at least this number of bases of the adapter are aligned to the read. This parameter corresponds to the `--stringency` option in Trim Galore and the `--overlap` option in [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html#random-matches).</small></details>| `integer` | 3 |  |  |
| `trim_error_rate` | Instructs Trim Galore to use this error rate when trimming. | `number` | 0.2 |  |  |
| `trim_q_cutoff` | Instructs Trim Galore to use this quality cutoff when trimming. | `integer` |  |  |  |
| `trim_minimum_length` | Instructs Trim Galore to remove reads shorter than this length after adapter/quality trimming. | `integer` | 5 |  |  |
| `skip_trimming` | Skip the adapter trimming step. <details><summary>Help</summary><small>Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no adapter contamination in your data.</small></details>| `boolean` |  |  |  |
| `save_trimmed` | Save the trimmed FastQ files in the results directory. <details><summary>Help</summary><small>By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.</small></details>| `boolean` |  |  |  |
| `min_trimmed_reads` | Minimum number of reads required after trimming. | `integer` | 1 |  |  |

#### Hard trimming options

Options to adjust hard trimming criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `hardtrim5_length` | Instructs Trim Galore to hard-trim (shorten) reads to this length from their 3' end. This is performed after adapter/quality trimming. | `integer` |  |  |  |
| `hardtrim3_length` | Instructs Trim Galore to hard-trim (shorten) reads to this length from their 5' end. This is performed after adapter/quality trimming. | `integer` |  |  |  |
| `save_hardtrimmed` | Whether to save the hard-trimmed reads (FASTQs) in the output directory. | `boolean` | true |  |  |

#### Alignment options

Options to adjust parameters and filtering criteria for read alignments.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `aligner` | Specifies the alignment algorithm to use. The available options are `bwa`, `bowtie2`, `chromap`, `star`, and `hisat2`. | `string` | chromap |  |  |
| `bwa_min_score` | Don't output BWA MEM alignments with score lower than this parameter. | `integer` |  |  |  |
| `sort_bam` | Should the BAM files be sorted by coordinate? Only relevant if Bowtie2 or BWA are used as aligners. | `boolean` | true |  |  |
| `save_align_intermeds` | Save the intermediate BAM files from the alignment step. <details><summary>Help</summary><small>By default, intermediate BAM files will not be saved. The final BAM files created after the appropriate filtering step are always saved to limit storage usage. Set this parameter to also save other intermediate BAM files.</small></details>| `boolean` |  |  |  |
| `save_unaligned` | Where possible, save unaligned reads from either STAR or HISAT2 to the results directory. <details><summary>Help</summary><small>This may either be in the form of FastQ or BAM files depending on the options available for that particular tool.</small></details>| `boolean` |  |  |  |
| `map_n_multimappers` | Number of multimappers to map. | `integer` |  |  |  |

#### Spike-in options

Options to adjust spike-in splitting criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `spikein_genome` | Name of the spike-in genome (if used in the experiment). The spike-in (exogenous) chromosome names in the `--fasta` and `<aligner>_index` files should have the format `<chromosome_name>_<spikein_genome>`, e.g. `chr1_dm6`. | `string` |  |  |  |
| `save_spikein_intermeds` | Save intermediate files during spike-in splitting and filtering of BAMs. These include the split and filtered BAM files before coordinate-sorting. | `boolean` | false |  |  |

#### Multimapper allocation options

Options to adjust multimapper allocation criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `multimap_allocation_method` | Method to allocate multimappers. Options are `allo`, `chromap`, and `mmr`. Not setting this parameter means that multimappers will not be allocated (the default). | `string` |  |  |  |
| `allo_mixed_cnn` | When using Allo for multimapper allocation, use the CNN trained on histone ChIP-seq datasets with mixed peaks, otherwise use the default (narrow). See [allo](https://github.com/seqcode/allo?tab=readme-ov-file#options) for more details. | `boolean` | true |  |  |
| `save_allocation_intermeds` | Save the intermediate BAM files from the multimapper allocation step. | `boolean` | true |  |  |
| `allocate_exogenous` | Whether to also allocate multimappers in the spike-in (exogenous) BAM files. | `boolean` | true |  |  |

## Alignment shifting options

Options to adjust alignment shifting criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `save_shifting_intermeds` | Whether to save the intermediate files from the alignment shifting (only done on ATAC-seq data). | `boolean` | false |  |  |

#### Filtering options

Options to adjust filtering criteria for BAM files.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_flT3` | Whether to skip the flT3 filtering step (filtering by mapping quality, which also removes multimapping reads in most cases). | `boolean` | false |  |  |

#### Coverage and normalization options

Options to adjust coverage and normalization criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_srpm` | Whether to skip the SRPM normalization step. See the [output documentation](https://github.com/grothlab/glseq/blob/dev/docs/output.md#normalized-bigwig-files) for more details. | `boolean` | false |  |  |
| `skip_cisrpm` | Whether to skip the CISRPM normalization step. See the [output documentation](https://github.com/grothlab/glseq/blob/dev/docs/output.md#normalized-bigwig-files) for more details. | `boolean` | false |  |  |
| `skip_cisrpmsoi` | Whether to skip the CISRPM SOI step. See the [output documentation](https://github.com/grothlab/glseq/blob/dev/docs/output.md#normalized-bigwig-files) for more details. | `boolean` | false |  |  |
| `soi_min_count` | Minimum CISRPM value required in a bin to be included in the signal over input (SOI) calculation. | `integer` | 0 |  |  |
| `coverage_extend_reads` | This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception. <details><summary>Help</summary><small>Single-end: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended. Paired-end: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads). Read the [deepTools documentation](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#read-processing-options) for more details.</small></details>| `integer` | 250 |  |  |
| `coverage_bin_size` | Size of the bins, in bases, for the output of the coverage calculation. | `integer` | 250 |  |  |
| `rpm_use_flT2_total` | Comma-separated list of histone marks for which RPM normalization factors should be calculated using the total number of reads before the quality filtering in flT3. The total mapped reads after flT1 will be used if flT2 was not performed (i.e., when there is no spike-in splitting). For marks/antibodies not in this list the total mapped reads after flT3 will be used for normalization factor calculation. | `string` | H3K9me3,H3K27me3 |  |  |
| `srpm_use_flT2_total` | Comma-separated list of histone marks for which SRPM normalization factors should be calculated using the total number of reads before the quality filtering in flT3. The total mapped reads after flT1 will be used if flT2 was not performed (i.e., when there is no spike-in splitting). For marks/antibodies not in this list the total mapped reads after flT3 will be used for normalization factor calculation. | `string` | H3K9me3,H3K27me3 |  |  |
| `cisrpm_use_flT2_total` | Comma-separated list of histone marks for which CISRPM normalization factors should be calculated using the total number of reads before the quality filtering in flT3. The total mapped reads after flT1 will be used if flT2 was not performed (i.e., when there is no spike-in splitting). For marks/antibodies not in this list the total mapped reads after flT3 will be used for normalization factor calculation. | `string` | H3K9me3,H3K27me3 |  |  |
| `save_coverage_intermeds` | Whether to save the intermediate files from the coverage calculation step. These include the coverage files in bedGraph format. | `boolean` | false |  |  |

#### Downsampling options

Options to adjust downsampling criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `bam_downsampling_method` | Method to downsample BAM files. Leave empty to disable downsampling. <details><summary>Help</summary><small>The `min_endo_across` method will downsample all IPs and inputs of the same antibody to the minimum number of endogenous reads across all IPs and inputs of the same antibody.<br>The `min_exo_across` method will downsample all IPs and inputs of the same antibody to the minimum number of exogenous reads across all IPs and inputs of the same antibody.<br>The `min_endo_by_type` method will downsample all IPs to the minimum number of endogenous reads across all IPs of the same antibody, and all the inputs to the minimum number of endogenous reads across all inputs of the same antibody.<br>The `min_exo_by_type` method will downsample all IPs to the minimum number of exogenous reads across all IPs of the same antibody, and all the inputs to the minimum number of exogenous reads across all inputs of the same antibody.</small></details>| `string` |  |  |  |
| `downsampling_endo_threshold` | Threshold of endogenous mapped reads for downsampling. For downsampling methods that search for the minimum endogenous mapped reads across a set of samples, this parameter will determine the lower limit of that minimum. In other words, the minimum mapped endogenous reads should be the minimum number above the `downsampling_endo_threshold`. If there are no totals above this threshold, then the (true) minimum is used. | `integer` | 10000000 |  |  |
| `downsampling_exo_threshold` | Threshold of exogenous mapped reads for downsampling. For downsampling methods that search for the minimum exogenous mapped reads across a set of samples, this parameter will determine the lower limit of that minimum. In other words, the minimum mapped exogenous reads should be the minimum number above the `downsampling_exogenous_threshold`. If there are no totals above this threshold, then the (true) minimum is used. | `integer` | 300000 |  |  |
| `dSp_use_flT2_total` | Comma-separated list of histone marks for which downsampling probabilities should be calculated using the total number of reads before the quality filtering in flT3. The total mapped reads after flT1 will be used if flT2 was not performed (i.e., when there is no spike-in splitting). For marks/antibodies not in this list, the total mapped reads after flT3 will be used for downsampling probability calculation. | `string` | H3K9me3,H3K27me3 |  |  |
| `bam_downsampling_strategy` | Strategy to downsample BAM files. See [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360036431292-DownsampleSam-Picard) for more details. | `string` | Chained |  |  |
| `bam_downsampling_accuracy` | Accuracy of downsampling BAM files. See [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360036431292-DownsampleSam-Picard) for more details. | `number` | 1e-06 |  |  |
| `save_downsampling_intermeds` | Save the intermediate BAM files from the downsampling step. | `boolean` | true |  |  |

#### Peak calling options

Options to adjust peak calling criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `narrow_peak` | Whether to call narrow peaks with MACS3 and Genrich. Default is to call broad peaks. | `boolean` | false |  |  |
| `peak_cutoff_method` | Whether to use q-value or p-value as the cutoff for peak calling. | `string` | qvalue |  |  |
| `narrow_cutoff` | Specifies narrow cutoff value (p-value or q-value) for MACS3 and Genrich peak calling. | `number` | 0.01 |  |  |
| `broad_cutoff` | Specifies broad cutoff value (p-value or q-value) for MACS3 and Genrich peak calling. | `number` | 0.1 |  |  |
| `min_reps_consensus` | Number of biological replicates required from a given condition for a peak to contribute to a consensus peak. <details><summary>Help</summary><small>If you are confident you have good reproducibility amongst your replicates then you can increase the value of this parameter to create a 'reproducible' set of consensus peaks. For example, a value of 2 will mean peaks that have been called in at least 2 replicates will contribute to the consensus set of peaks, and as such peaks that are unique to a given replicate will be discarded.</small></details>| `integer` | 1 |  |  |
| `macs_extsize` |  | `integer` | 200 |  |  |
| `macs_slocal` |  | `integer` | 1500 |  |  |
| `save_macs_pileup` | Instruct MACS3 to create bedGraph files normalised to signal per million reads. | `boolean` | true |  |  |
| `skip_peak_qc` | Skip MACS3 peak QC plot generation. | `boolean` |  |  |  |
| `skip_peak_annotation` | Skip annotation of MACS3 and consensus peaks with HOMER. | `boolean` |  |  |  |
| `skip_consensus_peaks` | Skip consensus peak generation, annotation and counting. | `boolean` |  |  |  |
| `skip_bdgcmp` | Skip bdgcmp calculation with MACS3. | `boolean` | false |  |  |
| `save_macs_bdgcmp_intermeds` |  | `boolean` |  |  |  |
| `skip_edd` |  | `boolean` | true |  |  |
| `edd_bin_size` |  | `string` |  |  |  |
| `edd_gap_penalty` |  | `string` |  |  |  |
| `edd_num_trials` |  | `integer` | 10000 |  |  |
| `edd_fdr` |  | `number` | 0.05 |  |  |
| `genrich_narrow_gap` | Genrich narrow gap parameter. | `integer` | 100 |  |  |
| `genrich_broad_gap` | Genrich broad gap parameter. | `integer` | 100 |  |  |
| `save_genrich_pileup` | Instruct Genrich to save pileup files. | `boolean` | true |  |  |
| `save_genrich_pvalues` | Instruct Genrich to save p-values. | `boolean` | true |  |  |
| `save_genrich_intervals` | Instruct Genrich to save intervals. | `boolean` | true |  |  |
| `skip_genrich` | Skip Genrich peak calling. | `boolean` |  |  |  |
| `epic2_binsize` | Size of the windows to scan the genome in epic2. Bin size is the smallest possible island. | `integer` | 150 |  |  |
| `epic2_fdr` | In epic2 peak calling, remove all islands with an FDR above this cutoff. | `number` | 0.01 |  |  |
| `skip_epic2` | Skip epic2 peak calling. | `boolean` | false |  |  |
| `skip_consenrich` | Skip Consenrich (Genome-wide 'consensus' epigenomic state estimates and uncertainty metrics). | `boolean` |  |  |  |
| `skip_mace` | Skip MACE peak calling. | `boolean` | false |  |  |
| `mace_kmer_size` | Set the k-mer size for MACE peak calling. Kmer size [6,12] to correct nucleotide composition bias. kmerSize < 0.5*read_length. Larger KmerSize might make program slower. Set kmerSize = 0 to turn off nucleotide compsition bias correction. default=6. See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `integer` |  |  |  |
| `mace_bin` | Chromosome chunk size. Each chomosome will be cut into small chunks of this size. Decrease chunk size will save more RAM. default=100000 (bp). See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `integer` |  |  |  |
| `mace_depth` | Reference reads count (default = 10 million). Sequencing depth will be normailzed to this number, so that wig files are comparable between replicates. See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `integer` |  |  |  |
| `mace_qcut` | phred scaled mapping quality threshhold to determine “uniqueness” of alignments. default=30. See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `number` |  |  |  |
| `mace_norm_method` | methods (“EM”, “AM”, “GM”, or “SNR”) used to consolidate replicates and reduce noise. “EM” = Entropy weighted mean, “AM”=Arithmetic mean, “GM”=Geometric mean, “SNR”=Signal-to-noise ratio. default=EM. See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `string` |  |  |  |
| `mace_pvalue` | Pvalue cutoff for border detection and subsequent border pairing. default=0.05. See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `number` |  |  |  |
| `mace_max_dist` | Maximum distance allowed for border pairing. default=100 (bp). See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `integer` |  |  |  |
| `mace_fuzzy_size` | Peaks (i.e. borders) located closely within this window will be merged. default=5 (bp). See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `integer` |  |  |  |
| `mace_window_size` | Background window size used to determine background signal level. default=100 (bp). See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `integer` |  |  |  |
| `mace_signal_fold` | Minimum coverage signal used to build model (i.e. estimate optimal peak pair size). default=2.0. See the [MACE documentation](https://chipexo.sourceforge.net/) for more details. | `number` |  |  |  |
| `save_mace_intermeds` | Whether to save the intermediate files from the MACE peak calling step. These include the MACE peak files in bedGraph format. | `boolean` | false |  |  |

#### Process skipping options

Options to skip various steps within the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_fastqc` | Skip FastQC. | `boolean` |  |  |  |
| `skip_picard_metrics` | Skip Picard CollectMultipleMetrics. | `boolean` |  |  |  |
| `skip_preseq` | Skip Preseq. | `boolean` | true |  |  |
| `deseq2_vst` | Use vst transformation instead of rlog with DESeq2. <details><summary>Help</summary><small>See [DESeq2 docs](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization).</small></details>| `boolean` | true |  |  |
| `skip_plot_profile` | Skip deepTools plotProfile. | `boolean` |  |  |  |
| `skip_plot_fingerprint` | Skip deepTools plotFingerprint. | `boolean` |  |  |  |
| `skip_spp` | Skip Phantompeakqualtools. | `boolean` |  |  |  |
| `skip_deseq2_qc` | Skip DESeq2 PCA and heatmap plotting. | `boolean` |  |  |  |
| `skip_igv` | Skip IGV. | `boolean` |  |  |  |
| `skip_multiqc` | Skip MultiQC. | `boolean` |  |  |  |
| `skip_qc` | Skip all QC steps except for MultiQC. | `boolean` |  |  |  |

#### SCAR-seq (partition) and OK-seq (RFD) analysis options

Options to adjust SCAR-seq and OK-seq analysis criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `partition_window_size` |  | `integer` | 1000 |  |  |
| `partition_step_size` |  | `integer` | 1000 |  |  |
| `smooth_radius` |  | `integer` | 30 |  |  |
| `derivative_radius` |  | `integer` | 30 |  |  |
| `zero_crossing_radius` |  | `integer` | 1 |  |  |
| `partition_cpm_cutoff` |  | `number` | 0.3 |  |  |
| `partition_plot_range` |  | `integer` | 100 |  |  |
| `partition_exclude_chromosomes` |  | `string` | chrX,chrY,chrM |  |  |
| `save_partition_intermeds` |  | `boolean` |  |  |  |

#### Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | true |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | true |
| `config_profile_name` | Institutional config name. | `string` |  |  | true |
| `config_profile_description` | Institutional config description. | `string` |  |  | true |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | true |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | true |

#### Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | true |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | true |
| `fingerprint_bins` | Number of genomic bins to use when calculating deepTools fingerprint plot. | `integer` | 500000 |  | true |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | true |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | true |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | true |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | true |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | true |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | true |  | true |

#### Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `multiqc_logo` |  | `string` |  |  |  |
| `multiqc_methods_description` |  | `string` |  |  |  |
| `hook_url` |  | `string` |  |  |  |
| `version` |  | `boolean` |  |  |  |
| `pipelines_testdata_base_path` |  | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  |  |
| `local_testdata_base_path` |  | `string` | /maps/projects/dan1/data/Groth_group/SRP/glseq_testdata/ |  |  |
| `validationFailUnrecognisedParams` |  | `boolean` |  |  |  |
| `validationLenientMode` |  | `boolean` |  |  |  |
| `validationShowHiddenParams` |  | `boolean` |  |  |  |


### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull grothlab/glseq
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

<!-- First, go to the [nf-core/chipseq releases page](https://github.com/nf-core/chipseq/releases) and find the latest pipeline version - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`. Of course, you can switch to another version by changing the number after the `-r` flag. -->

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

### Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Custom configuration

#### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

#### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

#### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

#### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

#### Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```