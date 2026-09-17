#!/usr/bin/env python3
"""Segment SACD DAPI MIPs and export intensity-preserving ImageJ nucleus montages.

Requires cellpose >=4, numpy, scipy, tifffile, and roifile (the existing
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
Only edge-touching nuclei with bounding-box aspect ratio >1.2 are excluded.
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

import numpy as np
from roifile import ImagejRoi, ROI_TYPE, ROI_SUBTYPE, ROI_OPTIONS
from scipy import ndimage
import tifffile

SUFFIX = '-SACD-MIP-YX.tif'
FIELDS = ['nucleus_id', 'condition', 'fov', 'mask_label', 'area_pixels',
          'y0', 'x0', 'y1', 'x1', 'crop_y0', 'crop_x0', 'crop_y1', 'crop_x1',
          'touches_edge', 'boundary_distance_pixels', 'aspect_ratio', 'included', 'exclusion_reason',
          'crop_file', 'montage_row', 'montage_column', 'montage_y', 'montage_x']


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


def fingerprint(group):
    return {c: {'path': str(p), 'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns}
            for c, p in group['paths'].items()}


def build_montage(output, rows, labels, cal, columns):
    included = [r for r in rows if r['included']]
    if not included:
        # Remove only this script's previously generated montage when rebuilding.
        (output / 'all_nuclei_montage.tif').unlink(missing_ok=True)
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
        row['montage_row'], row['montage_column'] = current_row, column
        previous_condition = row['condition']
        column += 1
    shape = (len(labels), (current_row + 1)*tile_h, columns*tile_w)
    overlays = []
    ranges = [[math.inf, -math.inf] for _ in labels[:-1]] + [[0., 1.]]
    # Local disk backing prevents a large montage from exhausting RAM.
    with tempfile.TemporaryDirectory(prefix='sacd-montage-') as temp:
        montage = np.memmap(Path(temp)/'montage.dat', dtype=np.float32, mode='w+', shape=shape)
        montage[:] = 0
        for row in included:
            crop = tifffile.imread(output / row['crop_file'])
            y = row['montage_row']*tile_h + label_height
            x = row['montage_column']*tile_w + (tile_w-crop.shape[2])//2
            row['montage_y'], row['montage_x'] = y, x
            montage[:, y:y+crop.shape[1], x:x+crop.shape[2]] = crop
            for c in range(len(labels)-1):
                ranges[c][0] = min(ranges[c][0], float(crop[c].min()))
                ranges[c][1] = max(ranges[c][1], float(crop[c].max()))
            overlays.append(text_roi(row['nucleus_id'], row['montage_column']*tile_w+4,
                                     row['montage_row']*tile_h, tile_w-8))
        write_imagej(output/'all_nuclei_montage.tif', montage, labels, cal, overlays, ranges)
        del montage
    return {'shape_cyx': shape, 'tile_height': tile_h, 'tile_width': tile_w,
            'columns': columns, 'overlay_count': len(overlays)}


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input', type=Path)
    p.add_argument('--output', type=Path, help='Default: INPUT/individual nucleus')
    p.add_argument('--diameter', type=float, default=None, help='Optional nucleus diameter in input pixels')
    p.add_argument('--min-size', type=int, default=15, help='Cellpose minimum area in pixels (default: 15)')
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
    p.add_argument('--rebuild', action='store_true', help='Require cached masks; do not run Cellpose')
    return p


def run(args):
    if (args.padding < 0 or args.edge_margin < 0 or args.columns < 1 or args.min_size < 1 or args.edge_aspect_ratio < 1
            or (args.limit is not None and args.limit < 1)
            or (args.diameter is not None and args.diameter <= 0)):
        raise ValueError('Invalid padding, columns, minimum size, edge ratio, limit, or diameter')
    root = args.input.expanduser().resolve()
    groups, channels, cal = discover(root)
    if args.trial:
        seen = set()
        groups = [g for g in groups if g['condition'] not in seen and not seen.add(g['condition'])]
    if args.limit:
        groups = groups[:args.limit]
    output = (args.output or root/'individual nucleus').expanduser().resolve()
    (output/'nuclei').mkdir(parents=True, exist_ok=True)
    (output/'masks').mkdir(exist_ok=True)
    labels = channels + ['Nucleus mask']
    inference = {k: getattr(args, k) for k in ('diameter', 'min_size', 'flow_threshold', 'cellprob_threshold')}
    inference.update(model='cpsam', cellpose_version=importlib.metadata.version('cellpose'))
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
        print(f'[{i}/{len(groups)}] {fov}: {"cached mask" if cached else "segmenting"}', flush=True)
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
                print(f'Loading Cellpose-SAM on {device}', flush=True)
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
        print(f'  retained {kept}, excluded {len(fov_rows)-kept}; {summary["seconds"]:.1f}s', flush=True)
        atomic_json(output/'progress.json', summaries)
    print('Building montage...', flush=True)
    layout = build_montage(output, rows, labels, cal, args.columns)
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
    settings.update(status='complete', montage=layout, fovs=summaries,
                    retained=sum(r['included'] for r in rows), excluded=sum(not r['included'] for r in rows))
    atomic_json(output/'run_settings.json', settings)
    print(f'Done: {settings["retained"]} retained, {settings["excluded"]} excluded. {output}', flush=True)


if __name__ == '__main__':
    arguments = parser().parse_args()
    try:
        run(arguments)
    except (ValueError, OSError) as exc:
        raise SystemExit(f'Error: {exc}') from exc
