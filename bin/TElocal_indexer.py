#!/usr/bin/env python

"""
===============================================================================

TElocal_indexer.py

Downloaded on 2025-07-26 from:
https://www.dropbox.com/scl/fo/tuy2s2tghlnxlnpjjq047/AJeMmkt1E7sztPfyFE1PEvc?rlkey=uoig7akt7iu1w2h5i00w1fd9o&dl=0

See for more information:
https://github.com/mhammell-laboratory/TEtranscripts/issues/221#issuecomment-2642814628


This script includes source code derived from the TEtranscripts project:
https://github.com/mhammell-laboratory/TEtranscripts

Original authors: Hammell Laboratory
License: GNU General Public License v3.0 or later (GPL-3.0+)
If you use or modify this code, please cite the original TEtranscripts publication:

Jin, Y., Tam, O.H., Paniagua, E., & Hammell, M. (2015).
TEtranscripts: A package for including transposable elements in differential
expression analysis of RNA-seq datasets. Bioinformatics, 31(22), 3593-3599.
https://doi.org/10.1093/bioinformatics/btv422

For license details, see:
    https://github.com/mhammell-laboratory/TEtranscripts/blob/master/LICENSE

===============================================================================
"""

import sys
import logging
import subprocess

try:
    import cPickle as pickle
except ImportError:
    import pickle

import argparse
import os
import time

# from TEToolkit.TEindex import *
# from TEToolkit.IntervalTree import *
# from TEToolkit.GeneFeatures import *

from TElocal_Toolkit.TEindex import *
from TElocal_Toolkit.GeneFeatures import *

TEindex_BINSIZE = 500


def read_opts4(parser):
    args = parser.parse_args()
    if not os.path.isfile(args.afile):
        logging.error("No such file: %s ! \n" % args.afile)
        sys.exit(1)

    if args.itype.lower() not in ['gene', 'te']:
        logging.error("indexing mode %s not supported \n" % args.itype)
        sys.exit(1)

    # Level of logging for tool
    logging.basicConfig(
        level=(4 - args.verbose) * 10,
        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        stream=sys.stderr,
        filemode="w")

    args.error = logging.critical  # function alias
    args.warn = logging.warning
    args.debug = logging.debug
    args.info = logging.info

#    args.argtxt = "# ARGUMENTS LIST: \n# prefix = %s \n# file to index = \
#    %s \n# index type = %s" % (args.prefix, args.afile, args.itype)
    args.argtxt = "# ARGUMENTS LIST: \n# file to index = \
    %s \n# index type = %s" % (args.afile, args.itype)
    return args


def prepare_parser():
    desc = "Building an index for the genome or transposable element annotations file."

    exmp = "Example: TElocal_indexer --afile gene_annotation.gtf --itype gene "

    parser = argparse.ArgumentParser(prog='TElocal_indexer', description=desc, epilog=exmp)

    parser.add_argument('--afile', metavar='annotation-file', dest='afile', type=str, required=True,
                        help='file for indexing of annotations')
    parser.add_argument('--itype', metavar='index-type', dest='itype', type=str, required=True,
                        help='index type to build for this gtf (gene or TE)')
#    parser.add_argument('--project', metavar='name', dest='prefix', default='',
#                        help='Prefix for the index file. Default is empty string')
    parser.add_argument('--verbose', metavar='verbose', dest='verbose', type=int, default=2,
                        help='Set verbose level. 0: only show critical message, 1: show additional warning message, \
                        2: show process information, 3: show debug messages. DEFAULT:2')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

    return parser


def main():
    """Start TElocal_indexer......parse options......"""

    args = read_opts4(prepare_parser())

    info = args.info

    # Output arguments used for program
    info("\n" + args.argtxt + "\n")

    info("Processing %s annotation file ... \n" % args.itype)

    if args.itype.lower() == 'gene':
        try:
            info("Building gene index ....... \n")
            geneIdx = GeneFeatures(args.afile, "exon", "gene_id")
            info("Done building gene index ...... \n")
            filename = args.afile.split('/')[-1] + '.ind'

            with open(filename,'wb') as filehandle:
                pickle.dump(geneIdx, filehandle, protocol=2)
            info("Done saving gene index")
            info("Gene index can be used by TEtranscripts and TElocal\n")
            info("Gene index saved to %s\n" % filename)

        except:
            sys.stderr.write("Error in building gene index \n")
            sys.exit(1)

    elif args.itype.lower() == 'te':
        try:
            teIdx = TEfeatures()
            cur_time = time.time()
            te_tmpfile = '.' + str(cur_time) + '.te.gtf'
            subprocess.call(['sort -k 1,1 -k 4,4g ' + args.afile + ' >' + te_tmpfile], shell=True)
            info("\nBuilding TE index ....... \n")
            teIdx.build(te_tmpfile)
            subprocess.call(['rm -f ' + te_tmpfile], shell=True)
            info("Done building TE index ...... \n")
            filename = args.afile.split('/')[-1] + '.locInd'

            with open(filename, 'wb') as filehandle:
                pickle.dump(teIdx, filehandle, protocol=2)

            info("Done saving local TE index")
            info("TE index can be only used by TElocal\n")
            info("TE index saved to %s\n" % filename)

        except:
            sys.stderr.write("Error in building TE index \n")
            sys.exit(1)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt !\n")
        sys.exit(0)
