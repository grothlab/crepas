#!/usr/bin/env python3
"""
hr_repliseq_call_features.py

Call replication features from a Gaussian-smoothed, column-scaled high-resolution
Repli-seq array: initiation zones, timing transition regions, breakages,
termination sites and late constant timing regions.

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
  - `cluster_by_birch()`, `find_peaks()`, `find_valleys()`, `find_slopes()` and
    `get_time_label()` are taken unchanged from
    `HighResRepliSeq/findFeatures.py`; the rest of this file is new.
  - Each chromosome is clustered once and every feature type read off that one
    profile, instead of re-clustering per feature type, and detection runs per
    contiguous run of bins so no feature spans bins that were dropped.
  - Added breakage calling; added the size split of valleys into
    termination sites and late CTRs; added the per-region fork-speed estimate.
  - Output is one BED per feature and per partition with a QC summary, rather
    than one CSV per invocation.
--------------------------------------------------------------------------------
"""

import argparse
import sys

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------------
# From the original HighResRepliSeq/findFeatures.py, unchanged
# --------------------------------------------------------------------------------

def cluster_by_birch(array, n_clusters=None, threshold=0.5):
    """Birch cluster a repli-seq array of shape (fractions, genomic bins)."""
    from sklearn.cluster import Birch

    if not isinstance(array, np.ndarray) or array.ndim != 2:
        raise ValueError("Input must be a 2D NumPy array.")
    if not np.isfinite(array).all():
        raise ValueError("Input array contains NaN or infinite values.")
    birch = Birch(n_clusters=n_clusters, threshold=threshold)
    return birch.fit(array.T)


def find_peaks(ArgMaxArr):
    """Plateaus higher than both neighbours: initiation zones."""
    ArgMaxArr = np.array(ArgMaxArr)
    plateaus = []
    start = None
    for i in range(1, len(ArgMaxArr)):
        if ArgMaxArr[i] == ArgMaxArr[i - 1]:
            if start is None:
                start = i - 1
        else:
            if start is not None:
                end = i - 1
                left_lower = start == 0 or ArgMaxArr[start - 1] < ArgMaxArr[start]
                right_lower = i == len(ArgMaxArr) or ArgMaxArr[i] < ArgMaxArr[start]
                if left_lower and right_lower:
                    plateaus.append([start, end])
                start = None
        if (i == 0 or ArgMaxArr[i - 1] < ArgMaxArr[i]) and \
           (i == len(ArgMaxArr) - 1 or ArgMaxArr[i + 1] < ArgMaxArr[i]):
            plateaus.append([i, i])
    if start is not None:
        left_lower = start == 0 or ArgMaxArr[start - 1] < ArgMaxArr[start]
        right_lower = ArgMaxArr[-1] < ArgMaxArr[start]
        if left_lower and right_lower:
            plateaus.append([start, len(ArgMaxArr) - 1])
    return plateaus


def find_valleys(ArgMaxArr):
    """Plateaus lower than both neighbours: termination sites and late CTRs."""
    ArgMaxArr = np.array(ArgMaxArr)
    plateaus = []
    start = None
    for i in range(1, len(ArgMaxArr)):
        if ArgMaxArr[i] == ArgMaxArr[i - 1]:
            if start is None:
                start = i - 1
        else:
            if start is not None:
                end = i - 1
                left_upper = start == 0 or ArgMaxArr[start - 1] > ArgMaxArr[start]
                right_upper = i == len(ArgMaxArr) or ArgMaxArr[i] > ArgMaxArr[start]
                if left_upper and right_upper:
                    plateaus.append([start, end])
                start = None
        if (i == 0 or ArgMaxArr[i - 1] > ArgMaxArr[i]) and \
           (i == len(ArgMaxArr) - 1 or ArgMaxArr[i + 1] > ArgMaxArr[i]):
            plateaus.append([i, i])
    if start is not None:
        left_upper = start == 0 or ArgMaxArr[start - 1] < ArgMaxArr[start]
        right_upper = ArgMaxArr[-1] < ArgMaxArr[start]
        if left_upper and right_upper:
            plateaus.append([start, len(ArgMaxArr) - 1])
    return plateaus


def find_slopes(ArgMaxArr, direction='right'):
    """Segments of constant slope: timing transition regions."""
    ArgMaxArr = np.array(ArgMaxArr)
    slopes = []
    n = len(ArgMaxArr)
    i = 0
    while i < n - 1:
        while i < n - 1 and ArgMaxArr[i] == ArgMaxArr[i + 1]:
            i += 1
        left_slope_start = i
        slope_start = i + 1
        while slope_start < n - 1:
            if (direction == 'left' and ArgMaxArr[slope_start] > ArgMaxArr[slope_start - 1]) or \
               (direction == 'right' and ArgMaxArr[slope_start] < ArgMaxArr[slope_start - 1]):
                slope_start += 1
            elif ArgMaxArr[slope_start] == ArgMaxArr[slope_start - 1]:
                temp = slope_start
                while temp < n - 1 and ArgMaxArr[temp] == ArgMaxArr[temp + 1]:
                    temp += 1
                if temp < n - 1:
                    if (direction == 'left' and ArgMaxArr[slope_start] > ArgMaxArr[slope_start - 1]) or \
                       (direction == 'right' and ArgMaxArr[slope_start] < ArgMaxArr[slope_start - 1]):
                        slope_start = temp + 1
                    else:
                        break
                else:
                    break
            else:
                break
        slope_end = slope_start - 1
        if slope_end > left_slope_start:
            slopes.append([left_slope_start, slope_end])
        i = slope_end + 1
    return slopes


def get_time_label(S_frac):
    if S_frac <= 2:
        return 'early'
    elif 2 < S_frac <= 5:
        return 'earlymid'
    elif 5 < S_frac <= 8:
        return 'latemid'
    return 'late'


# --------------------------------------------------------------------------------
# Added for the pipeline
# --------------------------------------------------------------------------------

def find_breakages(ArgMaxArr, ttr_intervals=()):
    """Breakages: plateaus that step part-way down a slope.

    A plateau qualifies when one neighbour is earlier and the other later. Pass
    `ttr_intervals` to skip plateaus a TTR already covers.
    """
    ArgMaxArr = np.array(ArgMaxArr)
    n = len(ArgMaxArr)
    in_ttr = np.zeros(n, dtype=bool)
    for a, b in ttr_intervals:
        in_ttr[a:b + 1] = True

    breakages = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and ArgMaxArr[j + 1] == ArgMaxArr[i]:
            j += 1
        # interior plateaus only: a breakage needs a neighbour on each side
        if i > 0 and j < n - 1 and not in_ttr[i:j + 1].any():
            left_earlier = ArgMaxArr[i - 1] > ArgMaxArr[i]
            right_later = ArgMaxArr[j + 1] < ArgMaxArr[i]
            left_later = ArgMaxArr[i - 1] < ArgMaxArr[i]
            right_earlier = ArgMaxArr[j + 1] > ArgMaxArr[i]
            if (left_earlier and right_later) or (left_later and right_earlier):
                breakages.append([i, j])
        i = j + 1
    return breakages


def find_breakages_flanked(ArgMaxArr, ttr_intervals):
    """Breakages that sit between two TTRs.

    A plateau qualifies when both neighbours are in a TTR and no single TTR holds
    both, so it is the seam between two rather than a step inside one.
    """
    ArgMaxArr = np.array(ArgMaxArr)
    n = len(ArgMaxArr)
    # which TTRs cover each bin
    covering = [set() for _ in range(n)]
    for k, (a, b) in enumerate(ttr_intervals):
        for pos in range(a, min(b, n - 1) + 1):
            covering[pos].add(k)

    breakages = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and ArgMaxArr[j + 1] == ArgMaxArr[i]:
            j += 1
        if i > 0 and j < n - 1:
            left_earlier = ArgMaxArr[i - 1] > ArgMaxArr[i]
            right_later = ArgMaxArr[j + 1] < ArgMaxArr[i]
            left_later = ArgMaxArr[i - 1] < ArgMaxArr[i]
            right_earlier = ArgMaxArr[j + 1] > ArgMaxArr[i]
            steps_down = (left_earlier and right_later) or (left_later and right_earlier)
            left_ttrs, right_ttrs = covering[i - 1], covering[j + 1]
            if steps_down and left_ttrs and right_ttrs and not (left_ttrs & right_ttrs):
                breakages.append([i, j])
        i = j + 1
    return breakages


def features_to_rows(features, chrom, index, Arr):
    """
    Turn bin-index intervals into rows, dropping any bin with an empty fraction.

    Also returns the surviving bin-index intervals, so the partition is built from
    exactly the features that were emitted.
    """
    rows, kept = [], []
    for f in features:
        if len(np.where(Arr[:, f[0]] == 0)[0]) != 0:
            continue
        rows.append([chrom, index[f[0]][0], index[f[1]][1],
                     get_time_label(np.argmax(Arr[:, f[0]]))])
        kept.append((f[0], f[1]))
    return rows, kept


# Order the partition resolves overlaps in; later entries win. The raw calls
# overlap, so genome percentages only add up once each bin is counted once.
PARTITION_PRIORITY = ['TTR', 'breakage', 'CTR', 'termination', 'IZ']

# Shortest run of adjacent bins that can hold a feature
MIN_BLOCK_BINS = 3

# 'overlap' calls every qualifying plateau, so a breakage may sit inside a TTR.
# 'disjoint' keeps breakages out of bins that are part of a called TTR
# 'flanked': the plateau is the seam between two different TTRs
PARTITION_VIEWS = ('overlap', 'disjoint', 'flanked')


def merge_across_breakages(spans, breakages):
    """Join TTRs back together across the breakages that interrupt them.

    `find_slopes` only traces strictly monotone runs, so a flat stretch inside a
    transition is split off as a breakage and the transition comes back as two
    short slopes. Those flat stretches are where the fork covers the most ground,
    so measuring the slopes alone keeps the steep parts and drops the fast ones.

    Two consecutive same-direction slopes are merged when everything between them
    is breakage.
    """
    by_chrom = {}
    for row in breakages:
        by_chrom.setdefault(row[0], []).append((int(row[1]), int(row[2])))
    for chrom in by_chrom:
        by_chrom[chrom].sort()

    def gap_is_breakage(chrom, left, right):
        if left >= right:
            return True
        reach = left
        for b_start, b_end in by_chrom.get(chrom, ()):
            if b_end <= reach:
                continue
            if b_start > reach:
                return False
            reach = b_end
            if reach >= right:
                return True
        return reach >= right

    merged = []
    current = None
    for row in sorted(spans, key=lambda r: (r[0], int(r[1]))):
        chrom, start, end, timing, direction, first, last, _n = row
        if (current and chrom == current['chrom'] and direction == current['direction']
                and gap_is_breakage(chrom, current['end'], start)):
            if start > current['end']:
                current['breakages_spanned'] += 1
            current['end'] = max(current['end'], end)
            current['last_fraction'] = last
            current['segments'] += 1
        else:
            if current:
                merged.append(current)
            current = {'chrom': chrom, 'start': start, 'end': end, 'timing': timing,
                       'direction': direction, 'first_fraction': first, 'last_fraction': last,
                       'segments': 1, 'breakages_spanned': 0}
    if current:
        merged.append(current)
    return merged


def ttr_fork_speeds(spans, breakages, n_fractions, s_phase_hours):
    """Fork speed within each timing transition region, in kb/min.

    One fork crosses a transition region, so its length divided by the stretch of
    S phase it spans gives that fork's speed.

    Slopes are first joined across the breakages interrupting them, since a breakage
    is an initiation event inside the transition, not the end of one. Measuring raw
    slopes instead is degenerate: any region advancing one fraction per bin returns
    `bin size * fractions / S phase` whatever its length. `breakages_spanned` counts
    how many were joined over.
    """
    columns = ['chrom', 'start', 'end', 'timing', 'direction', 'first_fraction',
               'last_fraction', 'fractions_traversed', 'segments', 'breakages_spanned',
               'breakage_free', 'speed_kb_per_min']
    if not spans:
        return pd.DataFrame(columns=columns)

    merged = merge_across_breakages(spans, breakages)
    df = pd.DataFrame(merged)
    df['fractions_traversed'] = (df['last_fraction'] - df['first_fraction']).abs() + 1
    df['breakage_free'] = df['breakages_spanned'] == 0

    minutes = s_phase_hours * 60.0
    s_phase_traversed = df['fractions_traversed'] / float(n_fractions)
    df['speed_kb_per_min'] = ((df['end'] - df['start']) / 1000.0) / (s_phase_traversed * minutes)
    return df[columns].sort_values(['chrom', 'start'])


def contiguous_blocks(df):
    """Yield (chrom, block) for each run of genomically adjacent bins.

    A feature runs from its first bin's start to its last bin's end, so bins that are
    next to each other in the array but not in the genome would give one feature
    spanning the gap. Bins do go missing: the G1 correction drops uncovered ones.
    Working per contiguous block keeps every feature inside real sequence. With no
    gaps this gives one block per chromosome.
    """
    for chrom in df.index.get_level_values(0).drop_duplicates():
        _df = df.loc[chrom]
        starts = np.asarray(_df.index.get_level_values(0), dtype=np.int64)
        ends = np.asarray(_df.index.get_level_values(1), dtype=np.int64)
        # a break wherever the next bin does not begin where this one ends
        cuts = np.flatnonzero(starts[1:] != ends[:-1]) + 1
        for block in np.split(np.arange(len(_df)), cuts):
            if len(block):
                yield chrom, _df.iloc[block]


def merge_runs(labels, index, chrom):
    """Collapse a per-bin label array into contiguous intervals."""
    rows = []
    i = 0
    n = len(labels)
    while i < n:
        j = i
        while j + 1 < n and labels[j + 1] == labels[i]:
            j += 1
        if labels[i] is not None:
            rows.append([chrom, index[i][0], index[j][1], labels[i]])
        i = j + 1
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Call replication features from a high-resolution Repli-seq array.")
    parser.add_argument("-a", "--array", required=True,
                        help="Processed array CSV with chrom, start, end and one column per fraction")
    parser.add_argument("--n_clusters", type=int, default=None,
                        help="Number of Birch clusters [default: let Birch decide]")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Birch radius threshold [default: %(default)s]")
    parser.add_argument("--termination_max_size", type=int, default=100000,
                        help="Valleys up to this size are termination sites, larger ones are late constant timing regions [default: %(default)s]")
    parser.add_argument("--s_phase_hours", type=float, default=10.0,
                        help="Length of S phase in hours, used to turn the stretch of S phase a "
                             "timing transition region spans into a time and so its length into "
                             "a fork speed. 10 h is the value assumed by Zhao, Sasaki & Gilbert (2020), "
                             "lines it profiles [default: %(default)s]")
    parser.add_argument("--partition_priority", default=",".join(PARTITION_PRIORITY),
                        help="Class order for the one-class-per-bin partition, lowest priority first, comma separated. Only affects the partition; the raw calls are written whatever this is [default: %(default)s]")
    parser.add_argument("--sample_name", default="sample", help="Name to report this sample under")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument("-p", "--prefix", default="", help="Prefix for every output filename")
    args = parser.parse_args()

    def out_path(suffix):
        name = f"{args.prefix}.{suffix}" if args.prefix else suffix
        return f"{args.outdir.rstrip('/')}/{name}"

    priority = [k.strip() for k in args.partition_priority.split(",") if k.strip()]
    if sorted(priority) != sorted(PARTITION_PRIORITY):
        sys.exit(f"Error: --partition_priority must list each of {PARTITION_PRIORITY} exactly once")

    df = pd.read_csv(args.array, index_col=['chrom', 'start', 'end'])
    calls = {k: [] for k in ('IZ', 'TTR', 'breakage', 'termination', 'CTR')}
    # Breakages inside TTRs make the two overlap, and the partition then hands most
    # TTR bins to `breakage`. Excluding them is equally defensible, so all readings
    # are carried through rather than one being picked.
    calls_breakage_disjoint = []
    calls_breakage_flanked = []
    ttr_spans = []
    cluster_rank = []
    partitions = {view: [] for view in PARTITION_VIEWS}
    skipped_blocks = skipped_bins = 0

    for chrom, _df in contiguous_blocks(df):
        Arr = _df.values.T
        # too few bins to hold a peak, a slope or a valley
        if _df.shape[0] < MIN_BLOCK_BINS:
            skipped_blocks += 1
            skipped_bins += _df.shape[0]
            continue
        model = cluster_by_birch(Arr, n_clusters=args.n_clusters, threshold=args.threshold)
        if model is None:
            sys.exit(f"Error: Birch clustering failed on {chrom}")
        # centroid rank, negated so earlier bins score higher
        ArgMaxArr = np.array([np.argmax(model.subcluster_centers_, axis=1)[i] * -1
                              for i in model.labels_[:]])
        index = list(_df.index)
        # The fraction each bin's centroid sits in. Every feature is read off this
        # profile, so it is written out for plotting.
        for (start, end), rank in zip(index, -ArgMaxArr):
            cluster_rank.append([chrom, start, end, int(rank) + 1])

        kept = {}
        rows, kept['IZ'] = features_to_rows(find_peaks(ArgMaxArr), chrom, index, Arr)
        calls['IZ'] += rows

        # TTRs first: breakages are the qualifying plateaus these did not absorb
        kept['TTR'] = []
        ttr_intervals = []
        for direction in ('right', 'left'):
            slopes = find_slopes(ArgMaxArr, direction=direction)
            ttr_intervals += slopes
            rows, k = features_to_rows(slopes, chrom, index, Arr)
            for row in rows:
                calls['TTR'].append(row + [direction])
            kept['TTR'] += k
            # S phase traversed, first bin's fraction to last bin's. With the
            # region's length this gives fork speed.
            for row, (a, b) in zip(rows, k):
                first = int(np.argmax(Arr[:, a]))
                last = int(np.argmax(Arr[:, b]))
                ttr_spans.append(row + [direction, first + 1, last + 1,
                                        abs(last - first) + 1])

        rows, kept['breakage'] = features_to_rows(
            find_breakages(ArgMaxArr), chrom, index, Arr)
        calls['breakage'] += rows

        rows, kept_breakage_disjoint = features_to_rows(
            find_breakages(ArgMaxArr, ttr_intervals=ttr_intervals), chrom, index, Arr)
        calls_breakage_disjoint += rows

        rows, kept_breakage_flanked = features_to_rows(
            find_breakages_flanked(ArgMaxArr, ttr_intervals), chrom, index, Arr)
        calls_breakage_flanked += rows

        # V-shaped valleys are termination sites; the large U-shaped ones are late CTRs
        kept['termination'], kept['CTR'] = [], []
        rows, k = features_to_rows(find_valleys(ArgMaxArr), chrom, index, Arr)
        for row, interval in zip(rows, k):
            key = 'termination' if (row[2] - row[1]) <= args.termination_max_size else 'CTR'
            calls[key].append(row)
            kept[key].append(interval)

        for view in PARTITION_VIEWS:
            kept_view = dict(kept)
            if view == 'disjoint':
                kept_view['breakage'] = kept_breakage_disjoint
            elif view == 'flanked':
                kept_view['breakage'] = kept_breakage_flanked
            labels = np.full(len(index), None, dtype=object)
            for key in priority:
                for a, b in kept_view[key]:
                    labels[a:b + 1] = key
            partitions[view] += merge_runs(labels, index, chrom)

        print(f"[hr_repliseq_call_features] {chrom}:{index[0][0]}-{index[-1][1]}: "
              f"{len(index)} bins, {len(model.subcluster_centers_)} subclusters", file=sys.stderr)

    starts = np.asarray(df.index.get_level_values(1), dtype=np.int64)
    ends = np.asarray(df.index.get_level_values(2), dtype=np.int64)
    genome_bp = int((ends - starts).sum())
    summary = []
    for key, rows in calls.items():
        cols = ['chrom', 'start', 'end', 'timing'] + (['direction'] if key == 'TTR' else [])
        out = pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)
        if not out.empty:
            out = out.sort_values(['chrom', 'start'])
        # only breakages differ between readings, so they carry the view in the name
        name = f"hr_{key}.overlap.bed" if key == 'breakage' else f"hr_{key}.bed"
        out.to_csv(out_path(name), sep="\t", header=False, index=False)
        bp = int((out['end'] - out['start']).sum()) if not out.empty else 0
        summary.append((key, len(out), bp, bp / genome_bp if genome_bp else float('nan')))
        print(f"[hr_repliseq_call_features] {key}: {len(out)} features, "
              f"{100 * bp / genome_bp if genome_bp else 0:.1f}% of bins", file=sys.stderr)

    alt_breakage = {}
    for view, rows in (('disjoint', calls_breakage_disjoint), ('flanked', calls_breakage_flanked)):
        out = pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'timing'])
        if not out.empty:
            out = out.sort_values(['chrom', 'start'])
        out.to_csv(out_path(f"hr_breakage.{view}.bed"), sep="\t", header=False, index=False)
        bp = int((out['end'] - out['start']).sum()) if not out.empty else 0
        alt_breakage[view] = (len(out), bp)
        print(f"[hr_repliseq_call_features] breakage ({view}): {len(out)} features, "
              f"{100 * bp / genome_bp if genome_bp else 0:.1f}% of bins", file=sys.stderr)

    ranks = pd.DataFrame(cluster_rank, columns=['chrom', 'start', 'end', 'cluster_fraction'])
    if not ranks.empty:
        ranks = ranks.sort_values(['chrom', 'start'])
    ranks.to_csv(out_path("hr_cluster_rank.tsv"), sep="\t", index=False)

    speeds = ttr_fork_speeds(ttr_spans, calls_breakage_flanked, df.shape[1], args.s_phase_hours)
    speeds.to_csv(out_path("hr_TTR.speed.tsv"), sep="\t", index=False)
    # `all` is every joined region; `uninterrupted` had no breakage to join over
    free = speeds[speeds['breakage_free']] if not speeds.empty else speeds
    for tag, sub in (('all', speeds), ('uninterrupted', free)):
        if sub.empty:
            print(f"[hr_repliseq_call_features] TTR speed ({tag}): none", file=sys.stderr)
            continue
        print(f"[hr_repliseq_call_features] TTR speed ({tag}): n={len(sub)} "
              f"median={sub['speed_kb_per_min'].median():.2f} kb/min", file=sys.stderr)

    part_summary = {}
    unassigned = {}
    for view in PARTITION_VIEWS:
        part_df = pd.DataFrame(partitions[view], columns=['chrom', 'start', 'end', 'feature'])
        if not part_df.empty:
            part_df = part_df.sort_values(['chrom', 'start'])
        part_df.to_csv(out_path(f"hr_partition.{view}.bed"), sep="\t", header=False, index=False)

        rows = []
        assigned_bp = 0
        for key in calls:
            sub = part_df[part_df['feature'] == key]
            bp = int((sub['end'] - sub['start']).sum()) if not sub.empty else 0
            assigned_bp += bp
            rows.append((key, len(sub), bp, bp / genome_bp if genome_bp else float('nan')))
            print(f"[hr_repliseq_call_features] partition[{view}] {key}: {len(sub)} regions, "
                  f"{100 * bp / genome_bp if genome_bp else 0:.1f}% of bins", file=sys.stderr)
        part_summary[view] = rows
        unassigned[view] = genome_bp - assigned_bp
        print(f"[hr_repliseq_call_features] partition[{view}] unassigned: "
              f"{100 * unassigned[view] / genome_bp if genome_bp else 0:.1f}% of bins",
              file=sys.stderr)

    with open(out_path("hr_features.qc.txt"), "w") as fh:
        fh.write(f"array:\t{args.array}\n")
        fh.write(f"sample:\t{args.sample_name}\n")
        fh.write(f"fractions:\t{df.shape[1]}\n")
        fh.write(f"bins:\t{df.shape[0]}\n")
        fh.write(f"birch_threshold:\t{args.threshold}\n")
        fh.write(f"birch_n_clusters:\t{args.n_clusters}\n")
        fh.write(f"termination_max_size:\t{args.termination_max_size}\n")
        fh.write(f"total_bp:\t{genome_bp}\n")
        fh.write(f"blocks_skipped_too_short:\t{skipped_blocks}\n")
        fh.write(f"bins_skipped_too_short:\t{skipped_bins}\n")
        for key, n, bp, frac in summary:
            fh.write(f"{key}_features:\t{n}\n")
            fh.write(f"{key}_bp:\t{bp}\n")
            fh.write(f"{key}_fraction_of_bp:\t{frac:.4f}\n")
        fh.write(f"partition_priority:\t{' > '.join(reversed(priority))}\n")
        fh.write(f"s_phase_hours:\t{args.s_phase_hours}\n")
        fh.write(f"ttr_speed_regions_after_merging:\t{len(speeds)}\n")
        fh.write(f"ttr_speed_breakages_stitched:\t"
                 f"{int(speeds['breakages_spanned'].sum()) if not speeds.empty else 0}\n")
        for tag, sub in (('all', speeds), ('uninterrupted', free)):
            fh.write(f"ttr_speed_{tag}_n:\t{len(sub)}\n")
            if sub.empty:
                continue
            values = sub['speed_kb_per_min']
            fh.write(f"ttr_speed_{tag}_median_kb_per_min:\t{values.median():.4f}\n")
            fh.write(f"ttr_speed_{tag}_q1_kb_per_min:\t{values.quantile(0.25):.4f}\n")
            fh.write(f"ttr_speed_{tag}_q3_kb_per_min:\t{values.quantile(0.75):.4f}\n")
            fh.write(f"ttr_speed_{tag}_min_kb_per_min:\t{values.min():.4f}\n")
            fh.write(f"ttr_speed_{tag}_max_kb_per_min:\t{values.max():.4f}\n")
        for view, (n, bp) in alt_breakage.items():
            fh.write(f"breakage_{view}_features:\t{n}\n")
            fh.write(f"breakage_{view}_bp:\t{bp}\n")
        for view in PARTITION_VIEWS:
            for key, n, bp, frac in part_summary[view]:
                fh.write(f"partition_{view}_{key}_regions:\t{n}\n")
                fh.write(f"partition_{view}_{key}_bp:\t{bp}\n")
                fh.write(f"partition_{view}_{key}_fraction_of_bp:\t{frac:.4f}\n")
            fh.write(f"partition_{view}_unassigned_bp:\t{unassigned[view]}\n")
            fh.write(f"partition_{view}_unassigned_fraction_of_bp:\t"
                     f"{unassigned[view] / genome_bp if genome_bp else float('nan'):.4f}\n")
            origin_free = sum(bp for key, _n, bp, _f in part_summary[view]
                              if key in ('TTR', 'termination'))
            fh.write(f"partition_{view}_origin_free_bp:\t{origin_free}\n")
            fh.write(f"partition_{view}_origin_free_fraction_of_bp:\t"
                     f"{origin_free / genome_bp if genome_bp else float('nan'):.4f}\n")

    header = ["Sample"] + [k for k, _, _, _ in summary]
    values = [args.sample_name] + [str(n) for _, n, _, _ in summary]
    with open(f"{args.outdir.rstrip('/')}/hr_repliseq_features_table.tsv", "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(values) + "\n")


if __name__ == "__main__":
    main()
