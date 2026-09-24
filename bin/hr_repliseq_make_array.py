#!/usr/bin/env python3
"""
hr_repliseq_make_array.py

Build a Gaussian-smoothed, column-scaled high-resolution Repli-seq array from
per-S-phase-fraction bin counts.

--------------------------------------------------------------------------------
Adapted for the grothlab/crepas pipeline (https://github.com/grothlab/crepas)
from High-Resolution-RepliSeq:

    https://github.com/oliviacamel/High-Resolution-RepliSeq
    Zhao PA, Sasaki T, Gilbert DM (2020). High-resolution Repli-Seq defines the
    temporal choreography of initiation, elongation and termination of
    replication in mammalian cells. Genome Biology 21:76.
    https://doi.org/10.1186/s13059-020-01983-8

The original work is licensed under the Apache License, Version 2.0:

    http://www.apache.org/licenses/LICENSE-2.0

Changes made to the original (Apache-2.0 section 4b):
  - `gaussian_smoothing()` and `scale()` are taken unchanged from
    `HighResRepliSeq/make_array.py`; the rest of this file is new.
  - Added a command line interface, deepTools counts-table input alongside the
    original's one bedGraph per fraction, the bin-size and sequence filtering its
    README applied by hand, and RPM/G1 normalization of the input fractions.

The G1 mappability correction has no counterpart in the original, whose released
bedGraphs are plain RPM; see `g1_filter` and `g1_ratio` below.
--------------------------------------------------------------------------------
"""

import argparse
import sys

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------------
# From the original HighResRepliSeq/make_array.py, unchanged
# --------------------------------------------------------------------------------

def gaussian_smoothing(M, kernel_shape=(3, 3), sigma=1):
    """Gaussian smoothing of a (fractions x bins) Repli-seq array."""

    def _gaussian_kernel(shape=kernel_shape, sigma=sigma):
        m, n = [(edge - 1) / 2 for edge in shape]
        y, x = np.ogrid[-m:m + 1, -n:n + 1]
        kernel = np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2))
        kernel[kernel < np.finfo(kernel.dtype).eps * kernel.max()] = 0
        kernel /= kernel.sum()
        return kernel

    newM = np.zeros_like(M)
    _M = np.concatenate((
        np.array([M[0, :] for _ in range((kernel_shape[0] - 1) // 2)]),
        M,
        np.array([M[-1, :] for _ in range((kernel_shape[0] - 1) // 2)]),
    ))
    _M = np.pad(_M, ((0, 0), ((kernel_shape[0] - 1) // 2, (kernel_shape[0] - 1) // 2)),
                'constant', constant_values=np.nan)

    for i in range((kernel_shape[0] - 1) // 2, _M.shape[0] - (kernel_shape[0] - 1) // 2):
        for j in range((kernel_shape[1] - 1) // 2, _M.shape[1] - (kernel_shape[0] - 1) // 2):
            box = np.ma.masked_invalid(
                _M[i - (kernel_shape[0] - 1) // 2:i + (kernel_shape[0] - 1) // 2 + 1,
                   j - (kernel_shape[1] - 1) // 2:j + (kernel_shape[1] - 1) // 2 + 1])
            newM[i - (kernel_shape[0] - 1) // 2, j - (kernel_shape[1] - 1) // 2] = \
                np.nansum(np.multiply(box, _gaussian_kernel()))

    return newM


def scale(M):
    """Convert each bin's fractions to percentages summing to 100."""
    col_sums = np.nansum(M, axis=0, keepdims=True)
    col_sums[col_sums == 0] = np.nan
    scaled_M = (M / col_sums) * 100
    return np.nan_to_num(scaled_M, nan=0.0)


# --------------------------------------------------------------------------------
# Added for the pipeline
# --------------------------------------------------------------------------------

def to_rpm(df):
    """Reads per million, per fraction."""
    totals = df.sum(axis=0)
    if (totals == 0).any():
        empty = list(totals.index[totals == 0])
        sys.exit(f"Error: fraction(s) {', '.join(map(str, empty))} have no reads")
    return df / totals * 1e6


def g1_filter(df, g1_label, min_rpm):
    """Drop bins the G1 control does not cover, keeping fractions as RPM.

    This is the mappability correction of Zhao, Sasaki & Gilbert (2020):
    bins "filtered out as a result of mappability normalisation".
    A region the G1 control cannot sequence cannot be measured in the pull-down libraries
    either, so it is removed; the fractions that remain keep their own depth-normalized
    values.

    Prefer this to `g1_ratio` (see below), which empties the array.
    """
    rpm = to_rpm(df)
    keep = rpm[g1_label] > min_rpm
    dropped = int((~keep).sum())
    return rpm.loc[keep].drop(columns=[g1_label]), dropped


def g1_ratio(df, g1_label, min_rpm):
    """Replace each fraction with log2 of its ratio to the G1 control, floored at zero.

    A stricter correction: "The log2 ratio
    between RPM of BrdU pulldown and that of G1 WGS was calculated for each S phase
    fraction ... Bins with values below zero ... were converted to zero."
    """
    rpm = to_rpm(df)
    keep = rpm[g1_label] > min_rpm
    dropped = int((~keep).sum())
    rpm = rpm.loc[keep]
    control = rpm[g1_label]
    fractions = rpm.drop(columns=[g1_label])
    with np.errstate(divide='ignore'):
        ratio = np.log2(fractions.div(control, axis=0))
    return ratio.where(ratio > 0, 0.0), dropped


def read_bedgraphs(paths, labels):
    """One bedGraph per fraction, joined on (chrom, start, end)."""
    frames = [pd.read_table(p, header=None, index_col=[0, 1, 2]) for p in paths]
    df = pd.concat(frames, axis=1)
    df.columns = labels
    df.index.names = ['chrom', 'start', 'end']
    return df


def read_counts_table(path, labels):
    """deepTools multiBamSummary --outRawCounts table; its header line is skipped."""
    df = pd.read_table(path, header=None, skiprows=1)
    if df.shape[1] != 3 + len(labels):
        sys.exit(f"Error: {path} has {df.shape[1] - 3} count columns but {len(labels)} "
                 f"fraction labels were given")
    df.columns = ['chrom', 'start', 'end'] + labels
    return df.set_index(['chrom', 'start', 'end'])


def main():
    parser = argparse.ArgumentParser(
        description="Build a Gaussian-smoothed, column-scaled high-resolution Repli-seq array.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--counts", help="deepTools multiBamSummary --outRawCounts table, one column per fraction")
    source.add_argument("--bedgraphs", nargs="+", help="One bedGraph per fraction, earliest first")
    parser.add_argument("--labels", nargs="+", default=None,
                        help="Fraction names, earliest first [default: S1..SN]")
    parser.add_argument("--bin_size", type=int, default=None,
                        help="Keep only bins of exactly this width, as the original did for 50 kb bins [default: keep all]")
    parser.add_argument("--exclude_pattern", default=r"^chrUn|_random$|_alt$|_fix$|\.",
                        help="Regular expression matching sequence names to drop [default: %(default)s]")
    parser.add_argument("--chromosomes", nargs="+", default=None,
                        help="Restrict to these sequences [default: every sequence that survives --exclude_pattern]")
    parser.add_argument("--normalization",
                        choices=["auto", "none", "rpm", "g1_filter", "g1_ratio"], default="auto",
                        help="How to normalize the input fractions before the array is built. "
                             "'g1_filter' drops bins the G1 whole-genome control does not cover "
                             "and keeps the fractions as RPM. 'g1_ratio' additionally replaces "
                             "each fraction with log2 of its ratio to the control, floored at "
                             "zero, which is the stricter correction but "
                             "zeroes about half of every fraction and should not be used to call "
                             "replication features. 'auto' uses 'g1_filter' when a G1 column is "
                             "present, 'rpm' for --counts and 'none' for --bedgraphs, which the "
                             "original expects to be normalized already [default: %(default)s]")
    parser.add_argument("--g1_label", default="G1",
                        help="Label of the G1 whole-genome-sequencing control column [default: %(default)s]")
    parser.add_argument("--g1_min_rpm", type=float, default=0.0,
                        help="Drop bins whose G1 control RPM is at or below this, since their "
                             "ratio is undefined [default: %(default)s]")
    parser.add_argument("--exclude_chromosomes", default="chrX,chrY,chrM",
                        help="Comma-separated sequences to drop, matched by name. The default "
                             "removes the sex chromosomes, whose copy number differs from the "
                             "autosomes, and the mitochondrial genome, which does not replicate "
                             "on the nuclear schedule. Pass an empty string to keep everything "
                             "[default: %(default)s]")
    parser.add_argument("-o", "--output", required=True, help="Output array as CSV")
    parser.add_argument("--qc", default=None, help="Optional normalization QC report")
    args = parser.parse_args()

    n_fractions = len(args.bedgraphs) if args.bedgraphs else None
    if args.labels:
        labels = args.labels
    elif n_fractions:
        labels = [f"S{i}" for i in range(1, n_fractions + 1)]
    else:
        sys.exit("Error: --labels is required with --counts, so that the count columns can be named")

    df = (read_bedgraphs(args.bedgraphs, labels) if args.bedgraphs
          else read_counts_table(args.counts, labels))

    df = df.reset_index()
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    if args.bin_size:
        df = df[(df['end'] - df['start']) == args.bin_size]
    if args.chromosomes:
        df = df[df['chrom'].isin(args.chromosomes)]
    elif args.exclude_pattern:
        df = df[~df['chrom'].astype(str).str.contains(args.exclude_pattern, regex=True)]
    if args.exclude_chromosomes:
        drop = {c.strip() for c in args.exclude_chromosomes.split(',') if c.strip()}
        present = drop & set(df['chrom'].astype(str))
        df = df[~df['chrom'].astype(str).isin(drop)]
        if present:
            print(f"[hr_repliseq_make_array] dropped {', '.join(sorted(present))}", file=sys.stderr)
    if df.empty:
        sys.exit("Error: no bins left after filtering")
    df = df.set_index(['chrom', 'start', 'end'])

    # Replicates of one fraction arrive as repeated labels. The array has one column
    # per fraction, so their counts are summed before normalization.
    if len(set(labels)) != len(labels):
        pooled = df.T.groupby(level=0, sort=False).sum().T
        labels = list(dict.fromkeys(labels))
        df = pooled[labels]
        print(f"[hr_repliseq_make_array] pooled replicates into {len(labels)} fractions",
              file=sys.stderr)

    has_g1 = args.g1_label in df.columns
    normalization = args.normalization
    if normalization == 'auto':
        normalization = 'g1_filter' if has_g1 else ('rpm' if args.counts else 'none')
    if normalization.startswith('g1') and not has_g1:
        sys.exit(f"Error: --normalization {normalization} needs a '{args.g1_label}' column "
                 f"among the fractions")
    if not normalization.startswith('g1') and has_g1:
        sys.exit(f"Error: a '{args.g1_label}' column was given but --normalization is "
                 f"'{normalization}', which would treat the control as a timing fraction")

    n_bins_in = df.shape[0]
    dropped = 0
    if normalization.startswith('g1'):
        correct = g1_filter if normalization == 'g1_filter' else g1_ratio
        df, dropped = correct(df, args.g1_label, args.g1_min_rpm)
        labels = [label for label in labels if label != args.g1_label]
    elif normalization == 'rpm':
        df = to_rpm(df)
    print(f"[hr_repliseq_make_array] normalization: {normalization}", file=sys.stderr)
    if dropped:
        print(f"[hr_repliseq_make_array] dropped {dropped} of {n_bins_in} bins not covered by "
              f"{args.g1_label}", file=sys.stderr)
    if df.empty:
        sys.exit("Error: no bins left after normalization")

    if args.qc:
        with open(args.qc, 'w') as handle:
            handle.write(f"normalization\t{normalization}\n")
            handle.write(f"bins_before_normalization\t{n_bins_in}\n")
            handle.write(f"bins_dropped_no_g1_coverage\t{dropped}\n")
            handle.write(f"bins_analysed\t{df.shape[0]}\n")
            handle.write(f"fractions\t{len(labels)}\n")
            handle.write("fraction\tmean\tmedian\tzero_bins\n")
            for label in labels:
                column = df[label]
                handle.write(f"{label}\t{column.mean():.6f}\t{column.median():.6f}\t"
                             f"{int((column == 0).sum())}\n")

    # Follow the order --chromosomes was given in, so output can be genomic rather
    # than lexicographic.
    present = list(df.index.get_level_values(0).drop_duplicates())
    chrom_order = [c for c in args.chromosomes if c in set(present)] if args.chromosomes else present

    processed = []
    for chrom in chrom_order:
        _df = df.loc[[chrom]].droplevel(0).sort_index(level=0)
        # scale, smooth, then scale again, as in the original README
        M = np.nan_to_num(_df.values).T
        M = scale(M)
        M = gaussian_smoothing(M)
        M = scale(M)
        processed.append(pd.DataFrame(
            M.T,
            index=pd.MultiIndex.from_tuples([(chrom, s, e) for s, e in _df.index],
                                            names=['chrom', 'start', 'end']),
            columns=labels))
        print(f"[hr_repliseq_make_array] {chrom}: {_df.shape[0]} bins", file=sys.stderr)

    pd.concat(processed, axis=0).to_csv(args.output)
    print(f"[hr_repliseq_make_array] wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
