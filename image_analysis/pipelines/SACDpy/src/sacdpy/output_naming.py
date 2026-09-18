"""Read-only selection of recorded outputs when filename conventions change."""
import json
from pathlib import Path


COMPLETED_STATUSES = {"written", "skipped_existing", "migrated_order_root", "resumed_manifest"}


def manifest_record(root, relative_fov):
    path = Path(root) / "_processing" / "manifest.json"
    if not path.exists():
        return None
    records = [r for r in json.loads(path.read_text()).get("fovs", [])
               if r.get("relative_fov") == relative_fov]
    if len(records) > 1:
        raise ValueError(f"Ambiguous duplicate manifest records for {relative_fov}")
    return records[0] if records else None


def resolve_recorded_paths(canonical, legacy, record, recorded_paths, expected, *, dataset_name):
    """Use a complete, identified manifest set; never infer ownership from a name."""
    expected = json.loads(json.dumps(expected))  # Manifest JSON represents tuples as lists.
    canonical = tuple(Path(p) for p in canonical)
    root = canonical[0].parent.resolve()
    known = {Path(p).resolve() for p in legacy}
    # Also cover old aliases prefixed by the dataset name, with or without its date.
    # Do not use a wildcard suffix: another legitimate FOV could share that suffix.
    dataset_prefixes = {dataset_name}
    date, separator, remainder = dataset_name.partition("_")
    if separator and len(date) == 8 and date.isdigit():
        dataset_prefixes.add(remainder)
    known.update((root / f"{prefix}{separator}{p.name}").resolve()
                 for prefix in dataset_prefixes for separator in ("_", "-") for p in canonical)
    canonical_set = {p.resolve() for p in canonical}
    existing_legacy = {p for p in known - canonical_set if p.exists()}
    recorded = tuple(Path(p) for p in recorded_paths) if recorded_paths else ()
    recorded_set = {p.resolve() for p in recorded}
    existing_recorded = {p for p in recorded_set if p.exists()}
    complete = record is not None and record.get("status") in COMPLETED_STATUSES
    label = expected["relative_fov"]

    if not complete:
        if existing_legacy or existing_recorded - canonical_set:
            raise ValueError(f"Legacy outputs for {label} lack successful manifest provenance; "
                             "use a new output location or explicitly resolve the existing files")
        return canonical

    if not recorded or len(recorded) != len(canonical) or len(recorded_set) != len(recorded):
        raise ValueError(f"Incomplete or duplicate manifest output paths for {label}")
    if any(p.resolve().parent != root for p in recorded):
        raise ValueError(f"Manifest outputs for {label} are outside its configured output directory")
    for key, value in expected.items():
        observed = record.get(key)
        matches = (isinstance(observed, str) and Path(observed).resolve() == Path(value).resolve()
                   if key == "source_folder" else observed == value)
        if not matches:
            raise ValueError(f"Incompatible manifest {key} for {label}: {observed!r}; "
                             f"expected {value!r}. These outputs cannot be resumed.")
    alternatives = existing_legacy | {p for p in canonical_set if p.exists()}
    if alternatives - recorded_set:
        raise ValueError(f"Ambiguous legacy and current naming variants for {label}; "
                         "resolve the duplicate output sets before resuming")
    if existing_recorded and len(existing_recorded) != len(recorded):
        raise ValueError(f"Incomplete recorded output set for {label}; refusing to create duplicates")
    if not existing_recorded:
        # A missing legacy result is not silently reconstructed under a different name.
        if recorded_set != canonical_set:
            raise ValueError(f"Recorded legacy outputs are missing for {label}; "
                             "use a new output location for reconstruction")
        return canonical
    return recorded


def reserve_paths(owners, paths, owner):
    for path in paths:
        path = Path(path).resolve()
        if path in owners and owners[path] != owner:
            raise ValueError(f"Output path collision between {owners[path]} and {owner}: {path}")
        owners[path] = owner
