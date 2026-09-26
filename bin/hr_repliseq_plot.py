#!/usr/bin/env python3
"""
hr_repliseq_plot.py

Draw the high-resolution Repli-seq heatmap and summarize the replication features
called from it.

Written for the grothlab/crepas pipeline (https://github.com/grothlab/crepas). The
heatmap approximates the presented in:

    Zhao PA, Sasaki T, Gilbert DM (2020). High-resolution Repli-Seq defines the
    temporal choreography of initiation, elongation and termination of replication
    in mammalian cells. Genome Biology 21:76.
    https://doi.org/10.1186/s13059-020-01983-8

"""

import argparse
import json
import re
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch

# The early/late line: warm where the ratio is positive, cool where negative.
EL_EARLY_COLOR = '#b2182b'
EL_LATE_COLOR = '#2166ac'

# One colour per replication feature, kept consistent between the heatmap track, the
# counts bar chart and the size distributions.
FEATURE_COLORS = {
    'IZ': '#2c7fb8',
    'TTR': '#7fcdbb',
    'breakage': '#d95f02',
    'termination': '#756bb1',
    'CTR': '#c51b8a',
}


def parse_region(text):
    """`chrom:start-end` in bp."""
    match = re.fullmatch(r"([^:]+):([\d,_]+)-([\d,_]+)", text.strip())
    if not match:
        sys.exit(f"Error: --region must look like chrom:start-end, got '{text}'")
    chrom = match.group(1)
    start = int(match.group(2).replace(",", "").replace("_", ""))
    end = int(match.group(3).replace(",", "").replace("_", ""))
    if end <= start:
        sys.exit(f"Error: --region end must be greater than start, got '{text}'")
    return chrom, start, end


def subset_region(df, chrom, start, end):
    """Bins of the array overlapping the window, kept whole rather than cut."""
    sub = df[df.index.get_level_values(0) == chrom]
    starts = np.asarray(sub.index.get_level_values(1), dtype=np.int64)
    ends = np.asarray(sub.index.get_level_values(2), dtype=np.int64)
    return sub[(ends > start) & (starts < end)]


def clip_bed(bed, chrom, start, end):
    """Features overlapping the window, with their coordinates clipped to it."""
    if bed.empty:
        return bed
    sub = bed[(bed['chrom'] == chrom) & (bed['end'] > start) & (bed['start'] < end)].copy()
    sub['start'] = sub['start'].clip(lower=start)
    sub['end'] = sub['end'].clip(upper=end)
    return sub


def read_bedgraph(path):
    """A 4-column bedGraph of the early/late log2 ratio."""
    if path is None:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'value'])
    try:
        df = pd.read_table(path, header=None, comment='t')
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'value'])
    df = df.iloc[:, :4]
    df.columns = ['chrom', 'start', 'end', 'value']
    return df[pd.to_numeric(df['value'], errors='coerce').notna()].astype(
        {'start': int, 'end': int, 'value': float})


def read_calls(paths):
    """The per-feature BEDs; column four is the fraction the feature falls in."""
    frames = []
    for path in paths:
        # Breakage BEDs are `hr_<feature>.<view>.bed`. The view is not part of the
        # name; leaving it on drops the feature from every panel.
        name = path.split('hr_')[-1].split('.')[0]
        try:
            df = pd.read_table(path, header=None)
        except pd.errors.EmptyDataError:
            continue
        # The TTR bed has a fifth column, the direction. Dropping it silently mixes
        # leftward and rightward regions.
        has_direction = df.shape[1] >= 5
        df = df.iloc[:, :5] if has_direction else df.iloc[:, :4]
        df.columns = ['chrom', 'start', 'end', 'timing'] + (['direction'] if has_direction else [])
        if not has_direction:
            df = df.assign(direction=None)
        frames.append(df.assign(feature=name))
    if not frames:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'timing', 'direction', 'feature'])
    return pd.concat(frames, ignore_index=True)


def read_bed(path):
    """A feature BED: chrom/start/end/timing[/direction]."""
    if path is None:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'feature'])
    try:
        df = pd.read_table(path, header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'feature'])
    df = df.iloc[:, :4]
    df.columns = ['chrom', 'start', 'end', 'feature']
    return df


def plot_el_track(ax, el, extent):
    """The early/late log2 ratio over the same window, above the heatmap.
    """
    mid = ((el['start'] + el['end']) / 2e6).to_numpy()
    values = el['value'].to_numpy()
    if len(mid) > 1:
        points = np.array([mid, values]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # colour each segment by its sign, so the line changes colour at the crossing
        midpoints = (values[:-1] + values[1:]) / 2
        colours = np.where(midpoints >= 0, EL_EARLY_COLOR, EL_LATE_COLOR)
        ax.add_collection(LineCollection(segments, colors=colours, linewidths=0.7))
    span = float(np.nanmax(np.abs(values))) if len(values) else 1.0
    ax.set_ylim(-span * 1.1, span * 1.1)
    ax.axhline(0, color='0.4', linewidth=0.5, zorder=0)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylabel('log2(E/L)', fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    ax.tick_params(labelbottom=False)


def plot_chromosome(fig, chrom, sub, partitions, fractions, el=None):
    """Heatmap of one chromosome, with a partition track per reading below."""
    starts = sub.index.get_level_values(1).to_numpy()
    ends = sub.index.get_level_values(2).to_numpy()
    extent = [starts.min() / 1e6, ends.max() / 1e6, len(fractions) + 0.5, 0.5]

    views = list(partitions)
    has_el = el is not None and not el.empty
    heights = ([3] if has_el else []) + [10] + [1] * len(views)
    grid = fig.add_gridspec(len(heights), 2, height_ratios=heights,
                            width_ratios=[40, 1], hspace=0.25, wspace=0.02, bottom=0.18)
    row = 0
    eax = None
    if has_el:
        eax = fig.add_subplot(grid[0, 0])
        fig.add_subplot(grid[0, 1]).axis('off')
        row = 1
    ax = fig.add_subplot(grid[row, 0], sharex=eax) if eax is not None else fig.add_subplot(grid[row, 0])
    cax = fig.add_subplot(grid[row, 1])
    taxes = []
    for i, view in enumerate(views):
        taxes.append(fig.add_subplot(grid[row + 1 + i, 0], sharex=ax))
        fig.add_subplot(grid[row + 1 + i, 1]).axis('off')
    tax = taxes[-1]

    image = ax.imshow(sub.to_numpy().T, aspect='auto', origin='upper', extent=extent,
                      cmap='RdYlBu_r', interpolation='nearest')
    ax.set_yticks(range(1, len(fractions) + 1))
    ax.set_yticklabels(fractions, fontsize=6)
    ax.set_ylabel('S phase fraction')
    ax.tick_params(labelbottom=False)
    if has_el:
        plot_el_track(eax, el, extent)
        eax.set_title(f'{chrom}  ({len(sub)} bins)', fontsize=10)
    else:
        ax.set_title(f'{chrom}  ({len(sub)} bins)', fontsize=10)
    fig.colorbar(image, cax=cax, label='% of bin replicated')

    for axis, view in zip(taxes, views):
        for _, row in partitions[view].iterrows():
            axis.axvspan(row['start'] / 1e6, row['end'] / 1e6, ymin=0.15, ymax=0.85,
                         color=FEATURE_COLORS.get(row['feature'], '#999999'), linewidth=0)
        axis.set_xlim(extent[0], extent[1])
        axis.set_yticks([])
        axis.set_ylabel(view, fontsize=6, rotation=0, ha='right', va='center')
        if axis is not tax:
            axis.tick_params(labelbottom=False)
    tax.set_xlabel('Position (Mb)')
    present = [f for f in FEATURE_COLORS
               if any(f in set(p['feature']) for p in partitions.values())]
    if present:
        # On the figure rather than the track axis, so it sits below the x-axis label
        fig.legend(handles=[Patch(facecolor=FEATURE_COLORS[f], label=f) for f in present],
                   loc='lower center', ncol=len(present), frameon=False, fontsize=7)


def plot_counts(ax, calls):
    """How many features of each kind were called.
    """
    counts = calls['feature'].value_counts()
    keys = [f for f in FEATURE_COLORS if f in counts.index]
    ax.bar(keys, [counts[f] for f in keys], color=[FEATURE_COLORS[f] for f in keys])
    ax.set_ylabel('Features')
    ax.tick_params(axis='x', labelsize=7)
    ax.set_title('Replication features called (raw calls)', fontsize=10)
    for i, key in enumerate(keys):
        ax.text(i, counts[key], f'{counts[key]}', ha='center', va='bottom', fontsize=7)


def plot_coverage(axes, partitions, total_bp):
    """Share of the analysed sequence per feature, one panel per partition.
    """
    keys = [f for f in FEATURE_COLORS
            if any(f in set(p['feature']) for p in partitions.values())]
    if not keys:
        for ax in axes:
            ax.axis('off')
        return

    pcts = {}
    for view, part in partitions.items():
        sizes = part.assign(bp=part['end'] - part['start']).groupby('feature')['bp'].sum()
        pcts[view] = [100 * sizes.get(f, 0) / total_bp for f in keys]
    ceiling = max(max(v) for v in pcts.values()) * 1.18 or 1

    for ax, (view, pct) in zip(axes, pcts.items()):
        ax.bar(keys, pct, color=[FEATURE_COLORS[f] for f in keys])
        ax.set_ylim(0, ceiling)
        ax.set_ylabel('% of analysed bp')
        ax.set_title(f'Genome covered ({view})', fontsize=9)
        ax.tick_params(axis='x', labelsize=7)
        for i, value in enumerate(pct):
            ax.text(i, value, f'{value:.1f}', ha='center', va='bottom', fontsize=6)
    for ax in list(axes)[len(pcts):]:
        ax.axis('off')


def plot_sizes(ax, calls):
    """Size per feature, on a log axis: they span three orders of magnitude."""
    keys = [f for f in FEATURE_COLORS if f in set(calls['feature'])]
    data = [(calls.loc[calls['feature'] == f, 'end']
             - calls.loc[calls['feature'] == f, 'start']).to_numpy() / 1e3 for f in keys]
    if not data:
        ax.axis('off')
        return
    boxes = ax.boxplot(data, tick_labels=keys, patch_artist=True, showfliers=False)
    for patch, key in zip(boxes['boxes'], keys):
        patch.set_facecolor(FEATURE_COLORS[key])
    ax.set_yscale('log')
    ax.set_ylabel('Size (kb, log scale)')
    ax.tick_params(axis='x', labelsize=7)
    ax.set_title('Feature size distribution (raw calls)', fontsize=10)


def plot_timing(ax, partition_with_timing):
    """When in S phase each feature class is called, from the timing column."""
    if partition_with_timing.empty:
        ax.axis('off')
        return
    keys = [f for f in FEATURE_COLORS if f in set(partition_with_timing['feature'])]
    data = [pd.to_numeric(partition_with_timing.loc[partition_with_timing['feature'] == f, 'timing'],
                          errors='coerce').dropna().to_numpy() for f in keys]
    keys = [k for k, d in zip(keys, data) if len(d)]
    data = [d for d in data if len(d)]
    if not data:
        ax.axis('off')
        return
    boxes = ax.boxplot(data, tick_labels=keys, patch_artist=True, showfliers=False)
    for patch, key in zip(boxes['boxes'], keys):
        patch.set_facecolor(FEATURE_COLORS[key])
    ax.set_ylabel('S phase fraction')
    ax.set_title('Replication timing of each feature', fontsize=10)


def plot_ttr_speed(ax, speeds, sample):
    """Fork speed within timing transition regions.
    """
    if speeds.empty:
        ax.axis('off')
        ax.text(0.5, 0.5, 'No transition-region speeds were estimated',
                ha='center', va='center', transform=ax.transAxes)
        return
    groups = [('all', speeds['speed_kb_per_min'].to_numpy())]
    free = speeds[speeds['breakage_free']]['speed_kb_per_min'].to_numpy()
    if len(free):
        groups.append(('uninterrupted', free))

    data = [values for _, values in groups]
    parts = ax.violinplot(data, showmeans=False, showmedians=False, showextrema=False,
                          widths=0.8)
    for body in parts['bodies']:
        body.set_facecolor(FEATURE_COLORS['TTR'])
        body.set_edgecolor('black')
        body.set_linewidth(0.5)
        body.set_alpha(0.9)

    for i, (_, values) in enumerate(groups, start=1):
        q1, median, q3 = np.percentile(values, [25, 50, 75])
        ax.vlines(i, values.min(), values.max(), color='black', linewidth=1)
        ax.vlines(i, q1, q3, color='black', linewidth=5)
        ax.plot(i, median, marker='o', markersize=4, color='white',
                markeredgecolor='black', markeredgewidth=0.5, zorder=3)
        ax.text(i, q3, f'{median:.2f}', ha='center', va='bottom', fontsize=7)

    # A few one-fraction regions reach tens of kb/min and would flatten the rest,
    # so the axis stops at the 99th percentile and says so.
    hi = max(np.percentile(values, 99) for _, values in groups)
    top = max(hi * 1.15, 1e-9)
    clipped = sum(int((values > top).sum()) for _, values in groups)
    ax.set_ylim(0, top)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([name for name, _ in groups])
    ax.set_ylabel('Fork speed (kb/min)')
    ax.set_title(f'{sample}: fork speed within TTRs', fontsize=10)
    if clipped:
        ax.text(0.99, 0.97, f'{clipped} region(s) above {top:.1f} kb/min not shown',
                transform=ax.transAxes, ha='right', va='top', fontsize=6, color='0.35')


# `key` is the BED to draw from, `direction` narrows TTRs to one orientation.
EXAMPLE_PANELS = [
    ('IZ', 'IZ', None),
    ('leftward\nTTR', 'TTR', 'left'),
    ('rightward\nTTR', 'TTR', 'right'),
    ('breakage', 'breakage', None),
    ('termination\nsites\n(<=100 kb)', 'termination', None),
    ('late CTRs', 'CTR', None),
]


def feature_window(ranks, chrom, start, end, flank_bins=4):
    """One feature's bins plus a margin, and which of them are the feature."""
    sub = ranks[ranks['chrom'] == chrom]
    if sub.empty:
        return None, None
    inside = (sub['end'] > start) & (sub['start'] < end)
    if not inside.any():
        return None, None
    positions = np.flatnonzero(inside.to_numpy())
    lo = max(positions[0] - flank_bins, 0)
    hi = min(positions[-1] + flank_bins, len(sub) - 1)
    window = sub.iloc[lo:hi + 1]
    is_feature = ((window['end'] > start) & (window['start'] < end)).to_numpy()
    # a margin that spans a gap in the analysed bins would not be contiguous sequence
    starts = window['start'].to_numpy()
    ends = window['end'].to_numpy()
    if len(window) > 1 and np.any(starts[1:] != ends[:-1]):
        return None, None
    return window, is_feature


def plot_feature_example(ax_line, ax_map, window, is_feature, array, fractions, show_scale):
    """One example column: rank profile above, heatmap of the same bins below."""
    x = np.arange(len(window))
    y = window['cluster_fraction'].to_numpy()
    # black for the surrounding bins, red for the ones that make up the feature
    ax_line.plot(x, y, color='black', linewidth=1.2)
    # Marker as well as line: a breakage is often one bin, and a one-point line
    # draws nothing at all.
    masked = np.where(is_feature, y, np.nan)
    ax_line.plot(x, masked, color='#d62728', linewidth=1.8, marker='o', markersize=2.2,
                 markeredgewidth=0)
    ax_line.set_xlim(-0.5, len(window) - 0.5)
    # Fixed to the full range, not the data, so panels share a scale and a shallow
    # feature looks shallow. Inverted so S1 is up, as in the heatmap.
    ax_line.set_ylim(len(fractions) + 0.5, 0.5)
    ax_line.set_xticks([])
    if show_scale:
        ax_line.set_yticks([1, len(fractions)])
        ax_line.set_yticklabels([fractions[0], fractions[-1]], fontsize=5)
    else:
        ax_line.set_yticks([])
    for side in ('top', 'right'):
        ax_line.spines[side].set_visible(False)

    keys = list(zip(window['chrom'], window['start']))
    values = np.array([array.loc[k] for k in keys])
    image = ax_map.imshow(values.T, aspect='auto', cmap='Greys', origin='upper',
                          vmin=0, vmax=float(np.percentile(array.to_numpy(), 99.5)))
    ax_map.set_xticks([])
    ax_map.set_yticks([0, len(fractions) - 1])
    ax_map.set_yticklabels([fractions[0], fractions[-1]], fontsize=6)
    return image


def plot_feature_examples(fig, ranks, calls, array, fractions, rng, sample, flank_bins=4):
    """One randomly drawn instance of each feature shape."""
    grid = fig.add_gridspec(2, len(EXAMPLE_PANELS) + 1, height_ratios=[1, 2],
                            width_ratios=[1] * len(EXAMPLE_PANELS) + [0.08],
                            hspace=0.08, wspace=0.35, top=0.76, bottom=0.12)
    image = None
    for col, (title, key, direction) in enumerate(EXAMPLE_PANELS):
        pool = calls[calls['feature'] == key]
        if direction is not None and 'direction' in pool.columns:
            pool = pool[pool['direction'] == direction]
        ax_line = fig.add_subplot(grid[0, col])
        ax_map = fig.add_subplot(grid[1, col])
        ax_line.set_title(title, fontsize=7)
        drawn = False
        # margins that run into a gap are skipped, so try a few before giving up
        for idx in rng.permutation(len(pool))[:25] if len(pool) else []:
            row = pool.iloc[int(idx)]
            window, is_feature = feature_window(ranks, row['chrom'], row['start'], row['end'],
                                                flank_bins=flank_bins)
            if window is None:
                continue
            image = plot_feature_example(ax_line, ax_map, window, is_feature, array, fractions,
                                         show_scale=(col == 0))
            ax_map.set_xlabel(f"{row['chrom']}:{int(row['start']) / 1e6:.2f}-"
                              f"{int(row['end']) / 1e6:.2f} Mb", fontsize=4)
            drawn = True
            break
        if not drawn:
            ax_line.axis('off')
            ax_map.axis('off')
            ax_line.text(0.5, 0.5, 'none', ha='center', va='center', fontsize=7,
                         transform=ax_line.transAxes)
    if image is not None:
        cax = fig.add_subplot(grid[1, len(EXAMPLE_PANELS)])
        fig.colorbar(image, cax=cax).set_label('% replicated', fontsize=6)
        cax.tick_params(labelsize=5)
    fig.suptitle(f'{sample}: features in the high-resolution Repli-seq heatmap', fontsize=9,
                 y=0.97)


def write_speed_mqc_json(path, speeds, sample):
    """MultiQC custom content: fork speed per sample, for comparing samples.

    A box rather than a violin: MultiQC's violin plot expects one value per sample, while these
    are whole distributions. The violin with the median, quartiles and range is on the PDF.
    """
    data = {f'{sample}: all': sorted(speeds['speed_kb_per_min'].tolist())}
    free = speeds[speeds['breakage_free']]['speed_kb_per_min']
    if len(free):
        data[f'{sample}: uninterrupted'] = sorted(free.tolist())
    payload = {
        'id': 'hr_repliseq_ttr_speed',
        'section_name': 'High-resolution Repli-seq: fork speed in TTRs',
        'description': ('Fork speed within each timing transition region, in kb/min, from its '
                        'length and the stretch of S phase it spans. Regions whose slope is '
                        'interrupted by a breakage are shown separately, since the breakage marks '
                        'an initiation event rather than a travelling fork, so each region is '
                        'measured across the breakages interrupting it rather than being '
                        'discarded. <code>uninterrupted</code> is the subset that had none.'),
        'plot_type': 'box',
        'pconfig': {'id': 'hr_repliseq_ttr_speed_plot',
                    'title': 'High-resolution Repli-seq: fork speed within TTRs',
                    'xlab': 'Fork speed (kb/min)'},
        'data': data,
    }
    with open(path, 'w') as handle:
        json.dump(payload, handle, indent=2)


def write_mqc_json(path, calls, sample):
    """MultiQC custom content: an interactive box plot of feature sizes in kb.

    From the raw calls, so that it agrees with the feature counts table in the same report.
    """
    keys = [f for f in FEATURE_COLORS if f in set(calls['feature'])]
    data = {}
    for key in keys:
        sub = calls[calls['feature'] == key]
        data[f'{sample}: {key}'] = sorted(((sub['end'] - sub['start']) / 1e3).tolist())
    payload = {
        'id': 'hr_repliseq_feature_sizes',
        'section_name': 'High-resolution Repli-seq: feature sizes',
        'description': ('Size of each replication feature called from the 16-fraction Repli-seq '
                        'array, in kb. These are the raw calls, which may overlap each other, '
                        'matching the feature counts above; the partition BED assigns every bin '
                        'to exactly one feature instead.'),
        'plot_type': 'box',
        'pconfig': {'id': 'hr_repliseq_feature_sizes_plot', 'title': 'High-resolution Repli-seq: feature sizes',
                    'xlab': 'Size (kb)'},
        'data': data,
    }
    with open(path, 'w') as handle:
        json.dump(payload, handle, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description="Plot a high-resolution Repli-seq array and the features called from it.")
    parser.add_argument("-a", "--array", required=True, help="Array CSV from hr_repliseq_make_array.py")
    parser.add_argument("--partition", default=None,
                        help="hr_partition.overlap.bed, the one-class-per-bin partition in which a "
                             "breakage may sit inside a timing transition region")
    parser.add_argument("--partition_flanked", default=None,
                        help="hr_partition.flanked.bed, the partition in which a breakage is the "
                             "seam between two timing transition regions. Each partition given "
                             "gets its own panel and track")
    parser.add_argument("--partition_disjoint", default=None,
                        help="hr_partition.disjoint.bed, the same partition with breakages kept "
                             "out of bins a timing transition region already claims. Both are "
                             "drawn, since the choice between them moves most of the genome "
                             "between the two classes")
    parser.add_argument("--cluster_rank", default=None,
                        help="hr_cluster_rank.tsv, the S phase fraction of the cluster centroid "
                             "each bin was assigned to. With it, pages of example features are "
                             "drawn: the rank profile of one instance of each feature shape with "
                             "its bins in red, over the heatmap of the same bins")
    parser.add_argument("--feature_example_pages", type=int, default=5,
                        help="How many pages of randomly drawn example features to add, when "
                             "--cluster_rank is given [default: %(default)s]")
    parser.add_argument("--example_flank_bins", type=int, default=4,
                        help="How many bins to show either side of an example feature. Smaller "
                             "zooms in on the shape, larger gives it more context "
                             "[default: %(default)s]")
    parser.add_argument("--seed", type=int, default=1,
                        help="Seed for drawing those examples, so a rerun gives the same ones "
                             "[default: %(default)s]")
    parser.add_argument("--el_track", default=None,
                        help="Early/late log2 ratio bedGraph for the same sample. When given, it "
                             "is drawn as a line above each heatmap, so the two measurements of "
                             "the same sample can be read against each other")
    parser.add_argument("--speeds", default=None,
                        help="hr_TTR.speed.tsv, the per-transition-region fork speeds. Drawn as "
                             "a violin per group with the median, quartiles and range marked")
    parser.add_argument("--calls", nargs="*", default=[],
                        help="Per-feature BEDs (hr_IZ.bed, hr_TTR.bed, ...), read for their timing column")
    parser.add_argument("--chromosomes", nargs="+", default=None,
                        help="Chromosomes to draw a heatmap page for [default: all]")
    parser.add_argument("--region", default=None,
                        help="Draw a single window instead of whole chromosomes, as "
                             "chrom:start-end in bp (e.g. chr4:20000000-25000000). Bins are kept "
                             "when they overlap the window at all. Intended for looking at one "
                             "locus by hand; the pipeline does not set it")
    parser.add_argument("--mqc_chromosome", default=None,
                        help="Chromosome shown in the MultiQC heatmap image [default: the first drawn]")
    parser.add_argument("--sample_name", default="sample", help="Name to report this sample under")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument("-p", "--prefix", default="", help="Prefix for every output filename")
    args = parser.parse_args()

    def out_path(suffix):
        name = f"{args.prefix}.{suffix}" if args.prefix else suffix
        return f"{args.outdir.rstrip('/')}/{name}"

    df = pd.read_csv(args.array, index_col=['chrom', 'start', 'end'])
    if df.empty:
        sys.exit(f"Error: {args.array} has no bins")

    fractions = list(df.columns)

    partitions = {'overlap': read_bed(args.partition)}
    if args.partition_disjoint:
        partitions['disjoint'] = read_bed(args.partition_disjoint)
    if args.partition_flanked:
        partitions['flanked'] = read_bed(args.partition_flanked)
    calls = read_calls(args.calls)
    el = read_bedgraph(args.el_track)
    ranks = pd.read_table(args.cluster_rank) if args.cluster_rank else pd.DataFrame(
        columns=['chrom', 'start', 'end', 'cluster_fraction'])
    speeds = pd.read_table(args.speeds) if args.speeds else pd.DataFrame(
        columns=['chrom', 'start', 'end', 'breakage_free', 'speed_kb_per_min'])

    if args.region:
        chrom, start, end = parse_region(args.region)
        df = subset_region(df, chrom, start, end)
        if df.empty:
            sys.exit(f"Error: no bins of the array overlap {args.region}")
        partitions = {view: clip_bed(part, chrom, start, end)
                      for view, part in partitions.items()}
        calls = clip_bed(calls, chrom, start, end)
        if not speeds.empty:
            speeds = speeds[(speeds['chrom'] == chrom) & (speeds['end'] > start)
                            & (speeds['start'] < end)]
        if not el.empty:
            el = el[(el['chrom'] == chrom) & (el['end'] > start) & (el['start'] < end)]
        args.chromosomes = [chrom]
        print(f"[hr_repliseq_plot] {args.region}: {len(df)} bins", file=sys.stderr)

    present = list(df.index.get_level_values(0).drop_duplicates())
    chroms = [c for c in args.chromosomes if c in set(present)] if args.chromosomes else present

    with PdfPages(out_path("hr_repliseq.plots.pdf")) as pdf:
        bin_starts = np.asarray(df.index.get_level_values(1), dtype=np.int64)
        bin_ends = np.asarray(df.index.get_level_values(2), dtype=np.int64)
        total_bp = int((bin_ends - bin_starts).sum())
        nothing_called = calls.empty and all(p.empty for p in partitions.values())

        # Page 1: what was called, how big it is and when it replicates
        summary = plt.figure(figsize=(13, 4.6))
        axes = summary.subplots(1, 3)
        if nothing_called:
            for ax in axes:
                ax.axis('off')
            axes[0].text(0.5, 0.5, 'No features were called', ha='center', va='center')
        else:
            plot_counts(axes[0], calls)
            plot_sizes(axes[1], calls)
            plot_timing(axes[2], calls)
        summary.suptitle(f'{args.sample_name}: high-resolution Repli-seq features')
        summary.tight_layout()
        pdf.savefig(summary)
        plt.close(summary)

        # Page 2: one genome-coverage panel per partition, on a shared scale
        coverage = plt.figure(figsize=(13, 4.6))
        n_views = max(len(partitions), 1)
        cov_axes = coverage.subplots(1, n_views, squeeze=False)[0]
        if nothing_called:
            for ax in cov_axes:
                ax.axis('off')
            cov_axes[0].text(0.5, 0.5, 'No features were called', ha='center', va='center')
        else:
            plot_coverage(cov_axes, partitions, total_bp)
        coverage.suptitle(f'{args.sample_name}: genome covered by each feature, '
                          f'by reading of the breakage/TTR relationship')
        coverage.tight_layout()
        pdf.savefig(coverage)
        plt.close(coverage)

        # Page 3: fork speed within transition regions
        speed_fig = plt.figure(figsize=(7, 5.5))
        plot_ttr_speed(speed_fig.add_subplot(111), speeds, args.sample_name)
        speed_fig.tight_layout()
        pdf.savefig(speed_fig)
        plt.close(speed_fig)

        if not ranks.empty and not calls.empty and args.feature_example_pages > 0:
            indexed = df.set_index([df.index.get_level_values(0), df.index.get_level_values(1)])
            rng = np.random.default_rng(args.seed)
            for page in range(args.feature_example_pages):
                fig = plt.figure(figsize=(11, 4.2))
                plot_feature_examples(fig, ranks, calls, indexed, fractions, rng,
                                      args.sample_name, flank_bins=args.example_flank_bins)
                pdf.savefig(fig)
                plt.close(fig)
            print(f"[hr_repliseq_plot] drew {args.feature_example_pages} example-feature pages",
                  file=sys.stderr)

        for chrom in chroms:
            sub = df.loc[[chrom]]
            fig = plt.figure(figsize=(11, 5))
            plot_chromosome(fig, chrom, sub,
                            {v: p[p['chrom'] == chrom] for v, p in partitions.items()},
                            fractions, el=el[el['chrom'] == chrom] if not el.empty else None)
            pdf.savefig(fig)
            plt.close(fig)
            print(f"[hr_repliseq_plot] drew {chrom}", file=sys.stderr)

    mqc_chrom = args.mqc_chromosome or (chroms[0] if chroms else None)
    if mqc_chrom:
        sub = df.loc[[mqc_chrom]]
        fig = plt.figure(figsize=(11, 5))
        plot_chromosome(fig, mqc_chrom, sub,
                        {v: p[p['chrom'] == mqc_chrom] for v, p in partitions.items()},
                        fractions, el=el[el['chrom'] == mqc_chrom] if not el.empty else None)
        fig.suptitle(f'{args.sample_name}', fontsize=9)
        fig.savefig(out_path(f"hr_repliseq_heatmap_{mqc_chrom}_mqc.png"), dpi=150,
                    bbox_inches='tight')
        plt.close(fig)

    if not calls.empty:
        write_mqc_json(out_path("hr_repliseq_feature_sizes_mqc.json"), calls, args.sample_name)
    if not speeds.empty:
        write_speed_mqc_json(out_path("hr_repliseq_ttr_speed_mqc.json"), speeds, args.sample_name)

    print(f"[hr_repliseq_plot] wrote {out_path('hr_repliseq.plots.pdf')}", file=sys.stderr)



if __name__ == "__main__":
    main()
