#!/usr/bin/env python3
"""Segment SACD DAPI MIPs and export intensity-preserving ImageJ nucleus montages.

Requires cellpose >=4, numpy, scipy, matplotlib, tifffile, and roifile (the existing
``/opt/miniconda3/envs/cellpose/bin/python`` environment includes these).

Example::

    python cli_sacd_nuclei.py /path/to/SACD_dataset --trial
    python cli_sacd_nuclei.py /path/to/SACD_dataset
    python cli_sacd_nuclei.py /path/to/SACD_dataset --rebuild --columns 8

Inputs: <FOV>__DAPI-SACD-MIP-YX.tif and matching <FOV>__<channel>-SACD-MIP-YX.tif.
Outputs default to <input>/individual nucleus. TIFF channels are DAPI, other
SACD channels (alphabetical), then Nucleus mask (0/1). Fluorescence pixels are
never normalized on export. ImageJ overlays identify nuclei without changing
pixel values; Image > Overlay > Hide Overlay hides the labels. The CSV uses
zero-based coordinates with exclusive upper bounds and zero-based tile indices.
Masks smaller than --min-nucleus-area (default 10000 pixels) are excluded.
Edge-touching nuclei with bounding-box aspect ratio >1.2 are also excluded.
The expression-matched montage retains eligible SHA nuclei and eligible dRRM/
dIDR nuclei within the inclusive SHA masked mean-SPEN range.
A 2-pixel edge tolerance accounts for Cellpose masks that stop just inside the
image boundary; --edge-margin controls this tolerance (0 means exact contact).
Masks are cached with source fingerprints and inference settings; --rebuild
requires matching cached inputs but deliberately reuses previous inference.
"""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import tempfile
import time
from collections import Counter

import numpy as np
from roifile import ImagejRoi, ROI_TYPE, ROI_SUBTYPE, ROI_OPTIONS
from scipy import ndimage
import tifffile

SUFFIX = '-SACD-MIP-YX.tif'
FIELDS = ['nucleus_id', 'condition', 'fov', 'mask_label', 'area_pixels', 'area_um2', 'mean_spen_intensity',
          'passes_boundary', 'passes_size', 'expression_matched', 'expression_exclusion_reason',
          'y0', 'x0', 'y1', 'x1', 'crop_y0', 'crop_x0', 'crop_y1', 'crop_x1',
          'touches_edge', 'boundary_distance_pixels', 'aspect_ratio', 'included', 'exclusion_reason',
          'crop_file', 'montage_row', 'montage_column', 'montage_y', 'montage_x',
          'matched_montage_row', 'matched_montage_column', 'matched_montage_y', 'matched_montage_x']
FIELDS += ['sorted_montage_row', 'sorted_montage_column', 'sorted_montage_y',
           'sorted_montage_x', 'expression_rank_within_condition']


class PipelineProgress:
    """Render callback events with Rich, or concise text when Rich is unavailable."""

    def __init__(self, enabled=True):
        self.progress = None
        self.tasks = {}
        self.stages = {}
        self.last_text = {}
        if enabled:
            try:
                from rich.progress import (Progress, SpinnerColumn, TextColumn, BarColumn,
                                           MofNCompleteColumn, TimeElapsedColumn, TimeRemainingColumn)
                self.progress = Progress(SpinnerColumn(), TextColumn('{task.description}'),
                                         BarColumn(), MofNCompleteColumn(), TimeElapsedColumn(),
                                         TimeRemainingColumn(), auto_refresh=False)
            except ImportError:
                pass

    def __enter__(self):
        if self.progress:
            self.progress.start()
        return self

    def __exit__(self, kind, value, traceback):
        if kind is not None:
            self({'task': 'stage', 'stage': 'failed', 'description':
                  'Interrupted' if issubclass(kind, KeyboardInterrupt) else f'Failed: {value}',
                  'total': None, 'completed': 0})
        if self.progress:
            self.progress.stop()

    def __call__(self, event):
        key = event.get('task', 'stage')
        stage = event.get('stage', key)
        description = event['description']
        completed, total = event.get('completed', 0), event.get('total')
        if self.progress:
            if key not in self.tasks:
                self.tasks[key] = self.progress.add_task(description, total=total)
            elif self.stages.get(key) != stage:
                self.progress.reset(self.tasks[key], total=total, description=description)
            self.progress.update(self.tasks[key], description=description, completed=completed)
            self.stages[key] = stage
            self.progress.refresh()
        elif self.last_text.get(key) != description:
            print(description, flush=True)
            self.last_text[key] = description


def emit(callback, description, *, stage='stage', task='stage', completed=0, total=None):
    """Publish progress without tying processing to a particular display backend."""
    event = dict(description=description, stage=stage, task=task, completed=completed, total=total)
    if callback is not None:
        callback(event)
    elif total is None or completed in (0, total):
        print(description, flush=True)


def expression_sorted_rows(rows):
    """Return matched rows by condition, descending mean SPEN, then natural FOV/ID."""
    for row in rows:
        row['expression_rank_within_condition'] = ''
        for coordinate in ('row', 'column', 'y', 'x'):
            row[f'sorted_montage_{coordinate}'] = ''
    selected = sorted((r for r in rows if r['expression_matched']),
                      key=lambda r: (natural(r['condition']), -r['mean_spen_intensity'],
                                     natural(r['fov']), r['mask_label']))
    ranks = Counter()
    for row in selected:
        ranks[row['condition']] += 1
        row['expression_rank_within_condition'] = ranks[row['condition']]
    return selected


def natural(value):
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', str(value))]


def condition(fov):
    return fov.split('-DAPI-')[0]


def calibration(tif):
    tags = tif.pages[0].tags
    resolution = []
    for key in ('XResolution', 'YResolution'):
        value = tags[key].value if key in tags else (1, 1)
        resolution.append(float(value[0]) / value[1])
    return {'resolution': resolution,
            'unit': (tif.imagej_metadata or {}).get('unit', 'pixel')}


def discover(root):
    groups = []
    expected_channels = None
    expected_cal = None
    for dapi in root.glob('*__DAPI' + SUFFIX):
        fov = dapi.name.removesuffix('__DAPI' + SUFFIX)
        # Exact prefix matching avoids FOV-1/FOV-10 cross-pairing.
        pairs = {p.name[len(fov) + 2:-len(SUFFIX)]: p
                 for p in root.glob('*' + SUFFIX)
                 if p.name.startswith(fov + '__')}
        channels = ['DAPI'] + sorted(set(pairs) - {'DAPI'}, key=natural)
        if len(channels) < 2:
            raise ValueError(f'{fov}: no other SACD MIP channel found')
        if expected_channels is None:
            expected_channels = channels
        if channels != expected_channels:
            raise ValueError(f'{fov}: channels {channels} differ from {expected_channels}')
        shape = None
        for channel in channels:
            with tifffile.TiffFile(pairs[channel]) as tif:
                series = tif.series[0]
                if len(series.shape) != 2 or series.dtype != np.dtype('float32'):
                    raise ValueError(f'{pairs[channel]}: expected a float32 YX SACD MIP')
                if shape is not None and shape != series.shape:
                    raise ValueError(f'{fov}: channel dimensions differ')
                shape = series.shape
                cal = calibration(tif)
                if expected_cal is None:
                    expected_cal = cal
                if cal['unit'] != expected_cal['unit'] or not np.allclose(
                        cal['resolution'], expected_cal['resolution'], rtol=1e-7):
                    raise ValueError(f'{pairs[channel]}: incompatible spatial calibration')
        groups.append({'fov': fov, 'condition': condition(fov), 'paths': pairs,
                       'shape': shape})
    if not groups:
        raise ValueError(f'No *__DAPI{SUFFIX} inputs found in {root}')
    groups.sort(key=lambda g: (natural(g['condition']), natural(g['fov'])))
    return groups, expected_channels, expected_cal


def atomic_json(path, data):
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def text_roi(text, x, y, width, height=24):
    return ImagejRoi(
        version=228, roitype=ROI_TYPE.RECT, subtype=ROI_SUBTYPE.TEXT,
        options=ROI_OPTIONS.SUB_PIXEL_RESOLUTION,
        name=text, text=text, text_name='SansSerif', text_size=14,
        left=x, top=y, right=x + width, bottom=y + height,
        xd=float(x), yd=float(y), widthd=float(width), heightd=float(height),
        stroke_color=b'\xff\xff\xff\xff', c_position=0, z_position=0,
        t_position=0).tobytes()


def write_imagej(path, data, labels, cal, overlays=(), ranges=None):
    if ranges is None:
        ranges = [(float(np.min(c)), float(np.max(c))) for c in data[:-1]] + [(0., 1.)]
    metadata = {'axes': 'CYX', 'mode': 'composite', 'unit': cal['unit'],
                'Labels': labels, 'Ranges': tuple(np.asarray(ranges).ravel())}
    if overlays:
        metadata['Overlays'] = list(overlays)
    # ImageJ uses contiguous classic TIFF even above 4 GiB (not BigTIFF).
    temporary = path.with_suffix('.tif.tmp')
    tifffile.imwrite(temporary, data, imagej=True, photometric='minisblack',
                     resolution=cal['resolution'], metadata=metadata)
    temporary.replace(path)


def nucleus_rows(masks, fov, padding, threshold, edge_margin=2):
    rows = []
    height, width = masks.shape
    areas = np.bincount(masks.ravel())
    for label, box in enumerate(ndimage.find_objects(masks), start=1):
        if box is None:
            continue
        sy, sx = box
        h, w = sy.stop - sy.start, sx.stop - sx.start
        boundary_distance = min(sy.start, sx.start, height-sy.stop, width-sx.stop)
        edge = boundary_distance <= edge_margin
        ratio = max(h, w) / min(h, w)
        keep = not (edge and ratio > threshold)
        match = re.search(r'-FOV(?:-(\d+))?$', fov)
        fov_id = (match.group(1) or 'unnumbered') if match else fov
        nucleus_id = f'{condition(fov)}_FOV-{fov_id}_N{label:04d}'
        rows.append(dict(nucleus_id=nucleus_id, condition=condition(fov), fov=fov,
                         mask_label=label, area_pixels=int(areas[label]),
                         y0=sy.start, x0=sx.start, y1=sy.stop, x1=sx.stop,
                         crop_y0=max(0, sy.start-padding), crop_x0=max(0, sx.start-padding),
                         crop_y1=min(height, sy.stop+padding), crop_x1=min(width, sx.stop+padding),
                         touches_edge=edge, boundary_distance_pixels=boundary_distance,
                         aspect_ratio=ratio, included=keep,
                         exclusion_reason='' if keep else 'edge_aspect_ratio',
                         crop_file=f'nuclei/{nucleus_id}.tif' if keep else '',
                         montage_row='', montage_column='', montage_y='', montage_x=''))
    return rows


def physical_pixel_area(cal):
    factors = {'um': 1., 'µm': 1., 'μm': 1., 'micron': 1., 'microns': 1.,
               'nm': 0.001, 'mm': 1000.}
    factor = factors.get(cal['unit'])
    if factor is None or min(cal['resolution']) <= 0:
        raise ValueError('Physical calibration in um, nm, or mm is required for area_um2')
    return factor**2 / np.prod(cal['resolution'])


def measure_and_filter(rows, masks, spen, minimum_area, pixel_area_um2):
    for row in rows:
        region = np.s_[row['y0']:row['y1'], row['x0']:row['x1']]
        selected = masks[region] == row['mask_label']
        row['mean_spen_intensity'] = float(np.mean(spen[region][selected], dtype=np.float64))
        row['area_um2'] = float(row['area_pixels'] * pixel_area_um2)
        row['passes_boundary'] = row['included']
        row['passes_size'] = row['area_pixels'] >= minimum_area
        row['included'] = row['passes_boundary'] and row['passes_size']
        reasons = ([] if row['passes_boundary'] else ['edge_aspect_ratio'])
        if not row['passes_size']:
            reasons.append('below_min_nucleus_area')
        row['exclusion_reason'] = ';'.join(reasons)
        if not row['included']:
            row['crop_file'] = ''
        for coordinate in ('row', 'column', 'y', 'x'):
            row[f'matched_montage_{coordinate}'] = ''


def match_expression(rows):
    reference = [r['mean_spen_intensity'] for r in rows
                 if r['condition'] == 'SHA' and r['included']]
    if not reference:
        raise ValueError('No eligible SHA nuclei remain after boundary and size filtering')
    low, high = min(reference), max(reference)
    for row in rows:
        reason = ''
        if not row['included']:
            reason = 'failed_boundary_or_size_filter'
        elif row['condition'] not in ('SHA', 'dRRM', 'dIDR'):
            reason = 'not_expression_matching_condition'
        elif row['mean_spen_intensity'] < low:
            reason = 'below_SHA_mean_spen_min'
        elif row['mean_spen_intensity'] > high:
            reason = 'above_SHA_mean_spen_max'
        row['expression_matched'] = not reason
        row['expression_exclusion_reason'] = reason
    return [low, high]


def save_size_histogram(output, rows, minimum_area):
    os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'sacd-matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    conditions = sorted({r['condition'] for r in rows}, key=natural)
    maximum = max((r['area_pixels'] for r in rows), default=minimum_area)
    zoom_max = max(15000, minimum_area*1.5)
    bins = np.linspace(0, max(maximum, minimum_area)*1.05, 31)
    fig, axes = plt.subplots(len(conditions), 2, squeeze=False,
                             figsize=(11, 2.8*len(conditions)), constrained_layout=True)
    colors = ['#0072B2', '#D55E00', '#009E73']
    for i, c in enumerate(conditions):
        # Show the size distribution before size filtering, as in cutoff selection.
        areas = [r['area_pixels'] for r in rows if r['condition'] == c and r['passes_boundary']]
        for j, edges in enumerate((bins, np.linspace(0, zoom_max, 31))):
            ax = axes[i, j]
            ax.hist(areas, bins=edges, color=colors[i % len(colors)], edgecolor='white')
            ax.axvline(minimum_area, color='black', linestyle='--', label=f'Cutoff: {minimum_area:,}')
            ax.set(xlabel='Nuclear mask area (pixels)', ylabel='Nuclei', xlim=(0, edges[-1]),
                   title=f'{c}, n={len(areas)}' + (' — small-area zoom' if j else ''))
            ax.spines[['top', 'right']].set_visible(False)
            ax.legend(frameon=False)
    fig.suptitle('Mask sizes after boundary filtering, before size filtering')
    fig.savefig(output/'nucleus_size_histogram.png', dpi=160)
    plt.close(fig)


def fingerprint(group):
    return {c: {'path': str(p), 'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns}
            for c, p in group['paths'].items()}


def build_montage(output, rows, labels, cal, columns, filename='all_nuclei_montage.tif',
                  selection='included', prefix='montage', progress_callback=None):
    """Pack selected crops in input order and write coordinates under their own prefix."""
    included = [r for r in rows if r[selection]]
    if not included:
        # Remove only this script's previously generated montage when rebuilding.
        (output / filename).unlink(missing_ok=True)
        return None
    label_height, gutter = 26, 8
    tile_h = max(r['crop_y1'] - r['crop_y0'] for r in included) + label_height + gutter
    tile_w = max(max(r['crop_x1'] - r['crop_x0'] for r in included) + gutter,
                 max(len(r['nucleus_id']) for r in included) * 9 + gutter)
    current_row, column, previous_condition = -1, 0, None
    for row in included:
        if row['condition'] != previous_condition or column == columns:
            current_row += 1
            column = 0
        row[f'{prefix}_row'], row[f'{prefix}_column'] = current_row, column
        previous_condition = row['condition']
        column += 1
    shape = (len(labels), (current_row + 1)*tile_h, columns*tile_w)
    overlays = []
    ranges = [[math.inf, -math.inf] for _ in labels[:-1]] + [[0., 1.]]
    # Local disk backing prevents a large montage from exhausting RAM.
    with tempfile.TemporaryDirectory(prefix='sacd-montage-') as temp:
        montage = np.memmap(Path(temp)/'montage.dat', dtype=np.float32, mode='w+', shape=shape)
        montage[:] = 0
        emit(progress_callback, f'Assembling {filename}', stage=filename, total=len(included))
        for tile_index, row in enumerate(included, 1):
            crop = tifffile.imread(output / row['crop_file'])
            y = row[f'{prefix}_row']*tile_h + label_height
            x = row[f'{prefix}_column']*tile_w + (tile_w-crop.shape[2])//2
            row[f'{prefix}_y'], row[f'{prefix}_x'] = y, x
            montage[:, y:y+crop.shape[1], x:x+crop.shape[2]] = crop
            for c in range(len(labels)-1):
                ranges[c][0] = min(ranges[c][0], float(crop[c].min()))
                ranges[c][1] = max(ranges[c][1], float(crop[c].max()))
            overlays.append(text_roi(row['nucleus_id'], row[f'{prefix}_column']*tile_w+4,
                                     row[f'{prefix}_row']*tile_h, tile_w-8))
            emit(progress_callback, f'Assembling {filename}', stage=filename,
                 completed=tile_index, total=len(included))
        emit(progress_callback, f'Writing {filename}', stage=filename + ':write')
        write_imagej(output/filename, montage, labels, cal, overlays, ranges)
        del montage
        emit(progress_callback, f'Saved {filename}', stage=filename + ':done', completed=1, total=1)
    return {'shape_cyx': shape, 'tile_height': tile_h, 'tile_width': tile_w,
            'columns': columns, 'overlay_count': len(overlays)}


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input', type=Path)
    p.add_argument('--output', type=Path, help='Default: INPUT/individual nucleus')
    p.add_argument('--diameter', type=float, default=None, help='Optional nucleus diameter in input pixels')
    p.add_argument('--min-size', type=int, default=15, help='Cellpose minimum area in pixels (default: 15)')
    p.add_argument('--min-nucleus-area', type=int, default=10000,
                   help='Post-segmentation minimum mask area in pixels (default: 10000)')
    p.add_argument('--spen-channel', default='SPEN_JFX650',
                   help='Channel used for masked mean expression (default: SPEN_JFX650)')
    p.add_argument('--flow-threshold', type=float, default=0.4)
    p.add_argument('--cellprob-threshold', type=float, default=0.)
    p.add_argument('--padding', type=int, default=10)
    p.add_argument('--edge-aspect-ratio', type=float, default=1.2,
                   help='Exclude edge-touching masks above this bounding-box ratio (default: 1.2)')
    p.add_argument('--edge-margin', type=int, default=2,
                   help='Boundary tolerance in pixels for slightly inset Cellpose masks (default: 2)')
    p.add_argument('--columns', type=int, default=10)
    p.add_argument('--device', choices=['auto', 'cpu', 'mps', 'cuda'], default='auto')
    p.add_argument('--trial', action='store_true', help='One FOV per condition; subsequent full run reuses masks')
    p.add_argument('--limit', type=int, help='Process only the first N selected FOVs')
    p.add_argument('--no-progress', action='store_true', help='Use plain-text progress instead of Rich')
    p.add_argument('--rebuild', action='store_true', help='Require cached masks; do not run Cellpose')
    return p


def preflight(args):
    """Validate inputs and report cache reuse without writing files or loading Cellpose."""
    if (args.padding < 0 or args.edge_margin < 0 or args.columns < 1 or args.min_size < 1 or args.min_nucleus_area < 1 or args.edge_aspect_ratio < 1
            or (args.limit is not None and args.limit < 1)
            or (args.diameter is not None and args.diameter <= 0)):
        raise ValueError('Invalid padding, columns, minimum size, edge ratio, limit, or diameter')
    root = args.input.expanduser().resolve()
    groups, channels, cal = discover(root)
    if args.spen_channel not in channels:
        raise ValueError(f'Missing SPEN channel: {args.spen_channel}')
    pixel_area_um2 = physical_pixel_area(cal)
    if args.trial:
        seen = set()
        groups = [g for g in groups if g['condition'] not in seen and not seen.add(g['condition'])]
    if args.limit:
        groups = groups[:args.limit]
    if not any(g['condition'] == 'SHA' for g in groups):
        raise ValueError('No SHA reference FOVs selected; expression matching requires SHA')
    output = (args.output or root/'individual nucleus').expanduser().resolve()
    inference = {k: getattr(args, k) for k in ('diameter', 'min_size', 'flow_threshold', 'cellprob_threshold')}
    inference.update(model='cpsam', cellpose_version=importlib.metadata.version('cellpose'))
    for group in groups:
        mask_path = output/'masks'/f"{group['fov']}__nuclei.tif"
        cache_path = mask_path.with_suffix('.json')
        cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
        group['cached'] = bool(mask_path.exists() and cache.get('sources') == fingerprint(group)
                               and (args.rebuild or cache.get('inference') == inference))
        if args.rebuild and not group['cached']:
            raise ValueError(f"{group['fov']}: missing/stale cache; disable cache-only mode to segment")
    return dict(root=root, output=output, groups=groups, channels=channels, calibration=cal,
                pixel_area_um2=pixel_area_um2, inference=inference,
                fovs_by_condition=dict(Counter(g['condition'] for g in groups)),
                cached_masks=sum(g['cached'] for g in groups),
                requires_segmentation=sum(not g['cached'] for g in groups))


def run(args, progress_callback=None):
    """Run the complete pipeline; return settings and paths for notebook review."""
    emit(progress_callback, 'Validating inputs and cached masks', stage='preflight')
    plan = preflight(args)
    root, output, groups = plan['root'], plan['output'], plan['groups']
    channels, cal = plan['channels'], plan['calibration']
    pixel_area_um2, inference = plan['pixel_area_um2'], plan['inference']
    (output/'nuclei').mkdir(parents=True, exist_ok=True)
    (output/'masks').mkdir(exist_ok=True)
    labels = channels + ['Nucleus mask']
    old_files = set()
    index_path = output/'nucleus_index.csv'
    if index_path.exists():
        with index_path.open() as f:
            old_files = {r['crop_file'] for r in csv.DictReader(f) if r['crop_file']}
    settings = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    settings.update(input=str(root), output=str(output), channels=labels, calibration=cal,
                    inference=inference, selected_fovs=[g['fov'] for g in groups], status='running')
    atomic_json(output/'run_settings.json', settings)
    model = None
    rows, summaries, ids = [], [], set()
    emit(progress_callback, 'FOVs | cached 0', task='fovs', total=len(groups))
    for i, group in enumerate(groups, 1):
        start = time.monotonic()
        fov = group['fov']
        mask_path = output/'masks'/f'{fov}__nuclei.tif'
        cache_path = mask_path.with_suffix('.json')
        sources = fingerprint(group)
        cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
        cached = mask_path.exists() and cache.get('sources') == sources
        cached = cached and (args.rebuild or cache.get('inference') == inference)
        if args.rebuild and not cached:
            raise ValueError(f'{fov}: missing cached mask or changed source images; run without --rebuild')
        emit(progress_callback, f'{fov}: {"cached mask" if cached else "segmenting"}', stage=fov)
        images = [tifffile.imread(group['paths'][c]) for c in channels]
        if any(not np.isfinite(image).all() for image in images):
            raise ValueError(f'{fov}: non-finite source intensities')
        if cached:
            masks = tifffile.imread(mask_path)
        else:
            if model is None:
                os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'sacd-matplotlib'))
                import torch
                from cellpose import models
                device = args.device
                if device == 'auto':
                    device = 'cuda' if torch.cuda.is_available() else ('mps' if torch.backends.mps.is_available() else 'cpu')
                emit(progress_callback, f'Loading Cellpose-SAM on {device}', stage='model')
                model = models.CellposeModel(device=torch.device(device), gpu=device != 'cpu',
                                             pretrained_model='cpsam')
                settings['actual_device'] = device
            masks = model.eval(images[0], diameter=args.diameter, min_size=args.min_size,
                               flow_threshold=args.flow_threshold,
                               cellprob_threshold=args.cellprob_threshold, normalize=True)[0]
            masks = np.asarray(masks, dtype=np.uint32)
            temporary = mask_path.with_suffix('.tif.tmp')
            tifffile.imwrite(temporary, masks, photometric='minisblack',
                             resolution=cal['resolution'], metadata={'axes': 'YX', 'unit': cal['unit']})
            temporary.replace(mask_path)
            atomic_json(cache_path, {'sources': sources, 'inference': inference})
        if masks.shape != images[0].shape or masks.dtype.kind not in 'ui':
            raise ValueError(f'{fov}: invalid cached mask dimensions or dtype')
        fov_rows = nucleus_rows(masks, fov, args.padding, args.edge_aspect_ratio, args.edge_margin)
        measure_and_filter(fov_rows, masks, images[channels.index(args.spen_channel)],
                           args.min_nucleus_area, pixel_area_um2)
        for row in fov_rows:
            if row['nucleus_id'] in ids:
                raise ValueError(f'Duplicate nucleus ID: {row["nucleus_id"]}')
            ids.add(row['nucleus_id'])
            if not row['included']:
                continue
            region = np.s_[row['crop_y0']:row['crop_y1'], row['crop_x0']:row['crop_x1']]
            crop = np.stack([image[region] for image in images] +
                            [(masks[region] == row['mask_label']).astype(np.float32)])
            write_imagej(output/row['crop_file'], crop, labels, cal)
        rows.extend(fov_rows)
        kept = sum(r['included'] for r in fov_rows)
        summary = {'fov': fov, 'condition': group['condition'], 'retained': kept,
                   'excluded': len(fov_rows)-kept, 'cached': cached,
                   'inference': cache.get('inference', inference) if cached else inference,
                   'seconds': round(time.monotonic()-start, 2)}
        summaries.append(summary)
        emit(progress_callback, f'FOVs | cached {sum(x["cached"] for x in summaries)} | last retained {kept}',
             task='fovs', completed=i, total=len(groups))
        atomic_json(output/'progress.json', summaries)
    expression_range = match_expression(rows)
    sorted_rows = expression_sorted_rows(rows)
    layout = build_montage(output, rows, labels, cal, args.columns, progress_callback=progress_callback)
    matched_layout = build_montage(output, rows, labels, cal, args.columns,
                                   'exp_matched_nuclei_montage.tif',
                                   'expression_matched', 'matched_montage', progress_callback)
    sorted_layout = build_montage(output, sorted_rows, labels, cal, args.columns,
                                  'exp_matched_nuclei_montage-sorted.tif',
                                  'expression_matched', 'sorted_montage', progress_callback)
    emit(progress_callback, 'Saving histogram and index', stage='reports')
    save_size_histogram(output, rows, args.min_nucleus_area)
    temporary = index_path.with_suffix('.csv.tmp')
    with temporary.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(index_path)
    new_files = {r['crop_file'] for r in rows if r['crop_file']}
    for stale in old_files - new_files:
        path = output/stale
        if path.parent == output/'nuclei' and path.suffix == '.tif':
            path.unlink(missing_ok=True)
    counts = {c: {'retained': sum(r['included'] for r in rows if r['condition'] == c),
                  'expression_matched': sum(r['expression_matched'] for r in rows if r['condition'] == c),
                  'excluded': sum(not r['included'] for r in rows if r['condition'] == c)}
              for c in sorted({r['condition'] for r in rows}, key=natural)}
    settings.update(status='complete', montage=layout, matched_montage=matched_layout, sorted_montage=sorted_layout,
                    sha_mean_spen_range=expression_range, counts_by_condition=counts,
                    expression_matched=sum(r['expression_matched'] for r in rows), fovs=summaries,
                    retained=sum(r['included'] for r in rows), excluded=sum(not r['included'] for r in rows))
    atomic_json(output/'run_settings.json', settings)
    emit(progress_callback, f'Complete: {settings["retained"]} nuclei, {settings["expression_matched"]} expression-matched',
         stage='complete', completed=1, total=1)
    return {'settings': settings, 'paths': {name: output/name for name in (
        'all_nuclei_montage.tif', 'exp_matched_nuclei_montage.tif',
        'exp_matched_nuclei_montage-sorted.tif', 'nucleus_size_histogram.png',
        'nucleus_index.csv', 'run_settings.json', 'nuclei', 'masks')}}


if __name__ == '__main__':
    arguments = parser().parse_args()
    try:
        with PipelineProgress(enabled=not arguments.no_progress) as progress:
            run(arguments, progress_callback=progress)
    except (ValueError, OSError) as exc:
        raise SystemExit(f'Error: {exc}') from exc
