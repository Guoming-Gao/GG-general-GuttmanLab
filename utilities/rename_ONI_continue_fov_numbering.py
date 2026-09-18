#!/usr/bin/env python3
"""Continue ONI FOV numbering across two acquisition folders (preview by default)."""

import argparse
import json
import os
from pathlib import Path
import re
import sys


JOURNAL = '.oni_continue_fov_numbering.jsonl'
FOV = re.compile(r'(.+)-FOV(?:-([0-9]+))?\Z')


def discover(root):
    conditions = {}
    for path in sorted(root.iterdir()):
        match = FOV.fullmatch(path.name)
        if path.is_symlink() or not path.is_dir() or not match:
            continue
        condition, number = match.groups()
        index = int(number or 0)
        group = conditions.setdefault(condition, {})
        if index in group:
            raise ValueError(f'Duplicate FOV index {index}: {group[index]} and {path}')
        group[index] = path
    return conditions


def ordered_moves(pairs, occupied):
    """Prefer supplied order, but vacate source names before using them."""
    pending = [(src, dst) for src, dst in pairs if src != dst]
    destinations = [dst for _, dst in pending]
    if len(set(destinations)) != len(destinations):
        raise ValueError('Multiple sources map to the same destination')
    occupied = set(occupied)
    result = []
    while pending:
        for i, (src, dst) in enumerate(pending):
            if dst not in occupied:
                result.append((src, dst))
                occupied.remove(src)
                occupied.add(dst)
                pending.pop(i)
                break
        else:
            raise ValueError(f'Conflicting destination or rename cycle: {pending[0][1]}')
    return result


def tree_moves(root, old, new):
    # ONI embeds identifiers after underscores and before underscores/extensions.
    token = re.compile(r'(?<![A-Za-z0-9-])' + re.escape(old)
                       + r'(?![A-Za-z0-9]|-[0-9])')
    operations = []

    def visit(directory):
        entries = sorted(directory.iterdir())
        pairs = []
        for entry in entries:
            # Symlinks are neither followed nor renamed.
            if entry.is_symlink():
                continue
            if entry.is_dir():
                visit(entry)
            pairs.append((entry, entry.with_name(token.sub(lambda _: new, entry.name))))
        operations.extend(ordered_moves(pairs, entries))

    visit(root)
    return operations


def build_plan(earlier, later):
    earlier, later = Path(earlier).resolve(), Path(later).resolve()
    if not earlier.is_dir() or not later.is_dir():
        raise ValueError('Both inputs must be existing directories')
    if earlier == later or earlier in later.parents or later in earlier.parents:
        raise ValueError('Input directories must be distinct and must not be nested')
    before, after = discover(earlier), discover(later)
    summaries, pairs = [], []
    for condition in sorted(after):
        group = after[condition]
        if condition not in before:
            summaries.append(f'{condition}: later-only condition; unchanged')
            continue
        maximum = max(before[condition])
        if min(group) > maximum:
            summaries.append(f'{condition}: already beyond earlier maximum {maximum}; unchanged')
            continue
        assigned = [(group[index], later / f'{condition}-FOV-{maximum + 1 + rank}')
                    for rank, index in enumerate(sorted(group))]
        summaries.append(f'{condition}: {len(group)} FOVs -> '
                         f'FOV-{maximum + 1} through FOV-{maximum + len(group)}')
        pairs.extend(reversed(assigned))
    # Usually descending index order. Gaps can require a lower source to move first.
    top_moves = ordered_moves(pairs, later.iterdir())
    operations = []
    for source, destination in top_moves:
        operations.extend(tree_moves(source, source.name, destination.name))
        operations.append((source, destination))
    return later, summaries, operations


def check_journal(path, allow_completed=False):
    if not path.exists():
        return
    try:
        with path.open() as handle:
            records = [json.loads(line) for line in handle]
        complete = records and records[-1].get('event') == 'complete'
    except (OSError, ValueError):
        complete = False
    if not complete:
        raise ValueError(f'Incomplete journal: {path}. Inspect and resolve the recorded '
                         'operations, then remove the journal before applying again.')
    if allow_completed:
        return
    raise ValueError(f'A completed journal already exists: {path}. '
                     'Archive it before applying a new plan to this folder.')


def apply_plan(root, operations):
    journal = root / JOURNAL
    check_journal(journal, allow_completed=not operations)
    if not operations:
        return
    # Exclusive creation also prevents simultaneous apply runs.
    with journal.open('x') as handle:
        def record(event, **data):
            handle.write(json.dumps(dict(event=event, **data)) + '\n')
            handle.flush()
            os.fsync(handle.fileno())

        record('plan', operations=[{'source': str(src.relative_to(root)),
                                   'destination': str(dst.relative_to(root))}
                                  for src, dst in operations])
        for index, (src, dst) in enumerate(operations):
            if not os.path.lexists(src) or os.path.lexists(dst):
                raise ValueError(f'Filesystem changed; cannot rename {src} -> {dst}')
            # Intent resolves the crash window between rename and completion log.
            record('intent', index=index)
            src.rename(dst)
            record('renamed', index=index)
        record('complete')


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog='Bare FOV is index 0. Later gaps are closed; conditions already entirely '
               'above the earlier maximum are skipped. Only names change, not file '
               'contents. Symlinks are skipped. Stop acquisition before applying. '
               f'Apply writes {JOURNAL}; incomplete journals block further applies.')
    parser.add_argument('earlier', type=Path)
    parser.add_argument('later', type=Path)
    parser.add_argument('--apply', action='store_true', help='perform renames (default: preview only)')
    args = parser.parse_args(argv)
    try:
        root, summaries, operations = build_plan(args.earlier, args.later)
        print('\n'.join(summaries) or 'No FOV directories found in the later folder.')
        for src, dst in operations:
            print(f'{src.relative_to(root)} -> {dst.relative_to(root)}')
        print(f'{len(operations)} renames; {"apply" if args.apply else "preview only"}.', flush=True)
        if args.apply:
            apply_plan(root, operations)
            print('Done.')
        return 0
    except (OSError, ValueError, KeyboardInterrupt) as error:
        print(f'Error: {error}', file=sys.stderr)
        if args.apply:
            print(f'If created, inspect journal: {args.later.resolve() / JOURNAL}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
