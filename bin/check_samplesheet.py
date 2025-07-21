#!/usr/bin/env python3

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat grothlab/glseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq_1,fastq_2,fastq_umi,okseq_part_file,replicate,exp_type,strandedness,antibody,control,control_replicate
    condition_1_H3K9me3,condition_1_bRep1_H3K9me3_R1.fastq.gz,condition_1_bRep1_H3K9me3_R3.fastq.gz,condition_1_bRep1_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,condition_1_INPUT,1
    condition_1_H3K9me3,condition_1_bRep2_H3K9me3_R1.fastq.gz,condition_1_bRep2_H3K9me3_R3.fastq.gz,condition_1_bRep2_H3K9me3_R2.fastq.gz,,2,chipseq,,H3K9me3,condition_1_INPUT,2
    condition_1_H3K27ac,condition_1_bRep1_H3K27ac_R1.fastq.gz,condition_1_bRep1_H3K27ac_R3.fastq.gz,condition_1_bRep1_H3K27ac_R2.fastq.gz,,1,chipseq,,H3K27ac,condition_1_INPUT,1
    condition_1_H3K27ac,condition_1_bRep2_H3K27ac_R1.fastq.gz,condition_1_bRep2_H3K27ac_R3.fastq.gz,condition_1_bRep2_H3K27ac_R2.fastq.gz,,2,chipseq,,H3K27ac,condition_1_INPUT,2
    condition_1_INPUT,condition_1_bRep1_INPUT_R1.fastq.gz,condition_1_bRep1_INPUT_R3.fastq.gz,condition_1_bRep1_INPUT_R2.fastq.gz,,1,chipseq,,,,
    condition_1_INPUT,condition_1_bRep2_INPUT_R1.fastq.gz,condition_1_bRep2_INPUT_R3.fastq.gz,condition_1_bRep2_INPUT_R2.fastq.gz,,2,chipseq,,,,
    For an example see:
    https://github.com/grothlab/glseq/blob/main/assets/samplesheets/ex1_multiBioRep_samplesheet.csv
    """
    file_out_tmp = file_out + ".tmp"
    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:

        ## Check header
        MIN_COLS = 3
        HEADER = ["sample", "fastq_1", "fastq_2", "fastq_umi", "okseq_part_file", "replicate", "exp_type", "strandedness", "antibody", "control", "control_replicate"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, fastq_1, fastq_2, fastq_umi, okseq_part_file, replicate, exp_type, strandedness, antibody, control, control_replicate = lspl[: len(HEADER)]
            if sample.find(" ") != -1:
                print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
                sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2, fastq_umi]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )
            ## Check replicate column is integer
            if not replicate.isdecimal():
                print_error("Replicate id not an integer!", "Line", line)
                sys.exit(1)

            ## Check exp_type
            if exp_type not in ["chipseq", "atacseq", "scarseq", "chorseq", "ChIP-exo", "OK-seq"]:
                print_error("Experiment type not 'chipseq', 'atacseq', 'scarseq', 'chorseq', 'ChIP-exo', or 'OK-seq'!", "Line", line)
                sys.exit(1)

            # strandedness needs to  be specified for scarseq or OK-seq samples
            if exp_type in ["scarseq", "OK-seq"]:
                if not strandedness:
                    print_error("Strandedness must be specified (only) for scarseq or OK-seq samples!", "Line", line)
                    sys.exit(1)
                if strandedness not in ["forward", "reverse"]:
                    print_error("Strandedness should be either 'forward' or 'reverse'!", "Line", line)
                    sys.exit(1)

            ## Check antibody and control columns have valid values
            if antibody:
                if antibody.find(" ") != -1:
                    print(f"WARNING: Spaces have been replaced by underscores for antibody: {antibody}")
                    antibody = antibody.replace(" ", "_")
            ## a control sample is no longer mandatory
                # if not control:
                #     print_error(
                #         "Both antibody and control columns must be specified!",
                #         "Line",
                #         line,
                #     )

            if control:
                if control.find(" ") != -1:
                    print(f"WARNING: Spaces have been replaced by underscores for control: {control}")
                    control = control.replace(" ", "_")
                if not control_replicate.isdecimal():
                    print_error("Control replicate id not an integer!", "Line", line)
                    sys.exit(1)
                control = "{}_REP{}".format(control, control_replicate)
                if not antibody and exp_type == "chipseq":
                    print_error(
                        "Both antibody and control columns must be specified for ChIP-seq samples!",
                        "Line",
                        line,
                    )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2, fastq_umi, okseq_part_file, replicate, antibody, control]
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", fastq_1, fastq_2, fastq_umi, okseq_part_file, replicate, exp_type, strandedness, antibody, control]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", fastq_1, fastq_2, fastq_umi, okseq_part_file, replicate, exp_type, strandedness, antibody, control]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Check that all ATAC-seq samples are paired-end, otherwise the alignmentsieve step will fail
            if exp_type == "atacseq" and sample_info[0] == "1":
                print_error("ATAC-seq samples must be paired-end for the alignmentsieve step to work!", "Line", line)
                
            ## Auto-detect UMI fastq file
            if sample and fastq_umi:
                sample_info.insert(1, "1")
            else:
                sample_info.insert(1, "0")

            ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, fastq_umi, okseq_part_file, replicate, antibody, control ]]}
            replicate = int(replicate)
            sample_info = sample_info + lspl[len(HEADER) :]
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = {}
            if replicate not in sample_mapping_dict[sample]:
                sample_mapping_dict[sample][replicate] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample][replicate]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample][replicate].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out_tmp)
        make_dir(out_dir)
        with open(file_out_tmp, "w") as fout:
            fout.write(
                ",".join(
                    [
                        "sample",
                        "single_end",
                        "sep_umi_fq",
                        "fastq_1",
                        "fastq_2",
                        "fastq_umi",
                        "okseq_part_file",
                        "replicate",
                        "exp_type",
                        "strandedness",
                        "antibody",
                        "control",
                    ]
                )
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that replicate ids are in format 1..<num_replicates>
                uniq_rep_ids = sorted(list(set(sample_mapping_dict[sample].keys())))
                if len(uniq_rep_ids) != max(uniq_rep_ids) or 1 != min(uniq_rep_ids):
                    print_error(
                        "Replicate ids must start with 1..<num_replicates>!",
                        "Sample",
                        "{}, replicate ids: {}".format(sample, ",".join([str(x) for x in uniq_rep_ids])),
                    )
                    sys.exit(1)

                ## Check that multiple replicates are of the same datatype i.e. single-end / paired-end
                # if not all(
                #     x[0][0] == sample_mapping_dict[sample][1][0][0] for x in sample_mapping_dict[sample].values()
                # ):
                #     print_error(
                #         f"Multiple replicates of a sample must be of the same datatype i.e. single-end or paired-end!",
                #         "Sample",
                #         sample,
                #     )

                for replicate in sorted(sample_mapping_dict[sample].keys()):
                    ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                    if not all(
                        x[0] == sample_mapping_dict[sample][replicate][0][0]
                        for x in sample_mapping_dict[sample][replicate]
                    ):
                        print_error(
                            f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                            "Sample",
                            sample,
                        )

                    for idx, val in enumerate(sample_mapping_dict[sample][replicate]):
                        control = "_REP".join(val[-1].split("_REP")[:-1])
                        control_replicate = val[-1].split("_REP")[-1]
                        if control and (
                            control not in sample_mapping_dict.keys()
                            or int(control_replicate) not in sample_mapping_dict[control].keys()
                        ):
                            print_error(
                                f"Control identifier and replicate has to match a provided sample identifier and replicate!",
                                "Control",
                                val[4],
                            )
                        # TODO: improve this maybe? prepend exp_type to control if not empty
                        # val[7] must be changed if a column is added
                        if control:
                            sample_mapping_dict[sample][replicate][idx][-1] = "{}_{}".format(val[7], val[-1])

                    ## Write to file
                    for idx in range(len(sample_mapping_dict[sample][replicate])):
                        fastq_files = sample_mapping_dict[sample][replicate][idx]
                        exp_type = fastq_files[7]
                        sample_id = "{}_{}_REP{}_T{}".format(exp_type, sample, replicate, idx + 1)
                        if len(fastq_files) == 1:
                            fout.write(",".join([sample_id] + fastq_files) + ",\n")
                        else:
                            fout.write(",".join([sample_id] + fastq_files) + "\n")

    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")

    ########### TODO: TEMPORARY FIX TO TRACK CONTROL SAMPLES
    with open(file_out_tmp, "r", encoding="utf-8-sig") as fin:
        lines = fin.readlines()

    headers = lines[0].strip().split(",")
    rows = [x.strip().split(",") for x in lines[1:]]
    control_index, sample_index = headers.index("control"), headers.index("sample")
    control_samples = {x[control_index] for x in rows}

    headers.append("is_control")

    with open(file_out, "w") as fout:
        fout.write(",".join(headers) + "\n")
        for row in rows:
            # sample_mod is row[sample_index] with trailing "_T<digit>" removed
            sample_mod = row[sample_index].rsplit("_T", 1)[0]
            row.append("1" if sample_mod in control_samples else "0")
            fout.write(",".join(row) + "\n")
    ########### TODO: TEMPORARY FIX TO TRACK CONTROL SAMPLES

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
