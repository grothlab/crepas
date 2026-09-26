#!/usr/bin/env python3

"""
===============================================================================
macs3_merged_expand.py

Originally created June 29th 2018 by:
    - Harshil Patel <https://github.com/drpatelh>

With contributions from:
    - Jose Espinosa-Carrasco <https://github.com/JoseEspinosa>

Source:
    https://github.com/nf-core/chipseq/blob/76e2382b6d443db4dc2396e6831d1243256d80b0/bin/macs3_merged_expand.py

Adapted for the grothlab/crepas pipeline by:
    - Samuel Ruiz-Pérez <samper@cancer.dk>
    https://github.com/grothlab/crepas/

Description:
    Add sample boolean files and aggregate columns from merged MACS narrow or broad peak file.

===============================================================================
"""

import os
import re
import errno
import argparse

############################################
############################################
## HELPER FUNCTIONS
############################################
############################################



def parse_min_replicates_per_sample(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("min_replicates_per_sample must be an integer (use -1 for all replicates)")
    if ivalue == 0 or ivalue < -1:
        raise argparse.ArgumentTypeError("min_replicates_per_sample must be >= 1, or -1 for all replicates")
    return ivalue



def makedir(path):
    if path:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def get_sample_id(name):
    """
    Remove trailing MACS peak suffix: _peak_<number>

    Example:
        ChIP-seq_bRep_2.foo_peak_8
        -> ChIP-seq_bRep_2.foo
    """
    return re.sub(r"_peak_\d+$", "", name)


def get_group_id(sample_id):
    """
    Remove everything from and including '_bRep_' onward.

    Example:
        ChIP-seq_bRep_2.mLb.mmr.flT1.T2T_CHM13.flT2.flT3.flTbl.rmO
        -> ChIP-seq
    """
    return re.sub(r"_bRep_.*$", "", sample_id)


############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = "Add sample boolean files and aggregate columns from merged MACS narrow or broad peak file."
Epilog = (
    "Example usage: python macs3_merged_expand.py <MERGED_INTERVAL_FILE> "
    "<SAMPLE_NAME_LIST> <OUTFILE> --is_narrow_peak --min_replicates_per_sample 1\n"
    "Or require all replicates per group: --min_replicates_per_sample -1"
)

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument(
    "MERGED_INTERVAL_FILE",
    help="Merged MACS3 interval file created using linux sort and mergeBed."
)
argParser.add_argument(
    "SAMPLE_NAME_LIST",
    help=(
        "Comma-separated list of sample names as named in individual MACS3 broadPeak/"
        "narrowPeak output file, but without the trailing _peak_<N>. "
        "Example: SAMPLE_bRep_1.mLb.mmr,SAMPLE_bRep_2.mLb.mmr"
    ),
)
argParser.add_argument("OUTFILE", help="Full path to output file.")

## OPTIONAL PARAMETERS
argParser.add_argument(
    "-in",
    "--is_narrow_peak",
    dest="IS_NARROW_PEAK",
    help="Whether merged interval file was generated from narrow or broad peak files (default: False).",
    action="store_true",
)
argParser.add_argument(
    "-mr",
    "--min_replicates_per_sample",
    type=parse_min_replicates_per_sample,
    dest="MIN_REPLICATES_PER_SAMPLE",
    default=1,
    help="Minimum number of replicates per sample required to contribute to merged peak. Use -1 to require all replicates for each sample group (default: 1).",
)

args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

## MergedIntervalTxtFile is file created using commands below:
## 1) broadPeak
## sort -k1,1 -k2,2n <MACS_BROADPEAK_FILES_LIST> | \
## mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse \
## > merged_peaks.txt
##
## 2) narrowPeak
## sort -k1,1 -k2,2n <MACS_NARROWPEAK_FILE_LIST> | \
## mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse \
## > merged_peaks.txt


def macs3_merged_expand(MergedIntervalTxtFile, SampleNameList, OutFile, isNarrow=False, minReplicatesPerSample=1):

    makedir(os.path.dirname(OutFile))

    combFreqDict = {}
    totalOutIntervals = 0
    SampleNameList = sorted(SampleNameList)

    ## TOTAL REPLICATES AVAILABLE PER GROUP
    ## Example:
    ##   ChIP-seq_bRep_1.xxx -> group ChIP-seq
    ##   ChIP-seq_bRep_2.xxx -> group ChIP-seq
    sampleReplicateCount = {}
    for sample in SampleNameList:
        gID = get_group_id(sample)
        if gID not in sampleReplicateCount:
            sampleReplicateCount[gID] = set()
        sampleReplicateCount[gID].add(sample)

    sampleReplicateCount = {gID: len(sIDs) for gID, sIDs in sampleReplicateCount.items()}

    fin = open(MergedIntervalTxtFile, "r")
    fout = open(OutFile, "w")

    oFields = (
        ["chr", "start", "end", "interval_id", "num_peaks", "num_samples"]
        + [x + ".bool" for x in SampleNameList]
        + [x + ".fc" for x in SampleNameList]
        + [x + ".qval" for x in SampleNameList]
        + [x + ".pval" for x in SampleNameList]
        + [x + ".start" for x in SampleNameList]
        + [x + ".end" for x in SampleNameList]
    )
    if isNarrow:
        oFields += [x + ".summit" for x in SampleNameList]

    fout.write("\t".join(oFields) + "\n")

    while True:
        line = fin.readline()
        if not line:
            fin.close()
            fout.close()
            break

        lspl = line.strip().split("\t")

        chromID = lspl[0]
        mstart = int(lspl[1])
        mend = int(lspl[2])
        starts = [int(x) for x in lspl[3].split(",")]
        ends = [int(x) for x in lspl[4].split(",")]
        names = lspl[5].split(",")
        fcs = [float(x) for x in lspl[8].split(",")]
        pvals = [float(x) for x in lspl[9].split(",")]
        qvals = [float(x) for x in lspl[10].split(",")]

        summits = []
        if isNarrow:
            summits = [int(x) for x in lspl[11].split(",")]

        ## GROUP SAMPLES BY REMOVING EVERYTHING FROM '_bRep_' ONWARD
        groupDict = {}
        for name in names:
            sID = get_sample_id(name)
            gID = get_group_id(sID)
            if gID not in groupDict:
                groupDict[gID] = []
            if sID not in groupDict[gID]:
                groupDict[gID].append(sID)

        ## GET SAMPLES THAT PASS REPLICATE THRESHOLD
        passRepThreshList = []
        for gID, sIDs in groupDict.items():
            if minReplicatesPerSample == -1:
                required_reps = sampleReplicateCount.get(gID, len(sIDs))
            else:
                required_reps = minReplicatesPerSample

            if len(sIDs) >= required_reps:
                passRepThreshList += sIDs

        ## GET VALUES FROM INDIVIDUAL PEAK SETS
        fcDict = {}
        qvalDict = {}
        pvalDict = {}
        startDict = {}
        endDict = {}
        summitDict = {}

        for idx in range(len(names)):
            sample = get_sample_id(names[idx])
            if sample in passRepThreshList:
                if sample not in fcDict:
                    fcDict[sample] = []
                fcDict[sample].append(str(fcs[idx]))

                if sample not in qvalDict:
                    qvalDict[sample] = []
                qvalDict[sample].append(str(qvals[idx]))

                if sample not in pvalDict:
                    pvalDict[sample] = []
                pvalDict[sample].append(str(pvals[idx]))

                if sample not in startDict:
                    startDict[sample] = []
                startDict[sample].append(str(starts[idx]))

                if sample not in endDict:
                    endDict[sample] = []
                endDict[sample].append(str(ends[idx]))

                if isNarrow:
                    if sample not in summitDict:
                        summitDict[sample] = []
                    summitDict[sample].append(str(summits[idx]))

        samples = sorted(fcDict.keys())

        if samples:
            numSamples = len(samples)
            boolList = ["TRUE" if x in samples else "FALSE" for x in SampleNameList]
            fcList = [";".join(fcDict[x]) if x in samples else "NA" for x in SampleNameList]
            qvalList = [";".join(qvalDict[x]) if x in samples else "NA" for x in SampleNameList]
            pvalList = [";".join(pvalDict[x]) if x in samples else "NA" for x in SampleNameList]
            startList = [";".join(startDict[x]) if x in samples else "NA" for x in SampleNameList]
            endList = [";".join(endDict[x]) if x in samples else "NA" for x in SampleNameList]

            oList = [
                str(x)
                for x in [
                    chromID,
                    mstart,
                    mend,
                    "Interval_" + str(totalOutIntervals + 1),
                    len(names),
                    numSamples,
                ]
                + boolList
                + fcList
                + qvalList
                + pvalList
                + startList
                + endList
            ]

            if isNarrow:
                oList += [";".join(summitDict[x]) if x in samples else "NA" for x in SampleNameList]

            fout.write("\t".join(oList) + "\n")

            tsamples = tuple(sorted(samples))
            if tsamples not in combFreqDict:
                combFreqDict[tsamples] = 0
            combFreqDict[tsamples] += 1
            totalOutIntervals += 1

    ## WRITE FILE FOR INTERVAL INTERSECT ACROSS SAMPLES
    ## COMPATIBLE WITH UPSETR PACKAGE
    fout = open(OutFile[:-4] + ".intersect.txt", "w")
    combFreqItems = sorted([(combFreqDict[x], x) for x in combFreqDict.keys()], reverse=True)
    for k, v in combFreqItems:
        fout.write("%s\t%s\n" % ("&".join(v), k))
    fout.close()


############################################
############################################
## RUN FUNCTION
############################################
############################################

macs3_merged_expand(
    MergedIntervalTxtFile=args.MERGED_INTERVAL_FILE,
    SampleNameList=args.SAMPLE_NAME_LIST.split(","),
    OutFile=args.OUTFILE,
    isNarrow=args.IS_NARROW_PEAK,
    minReplicatesPerSample=args.MIN_REPLICATES_PER_SAMPLE,
)
