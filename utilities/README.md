# Utilities

## Continue ONI FOV numbering across midnight

`oni_continue_fov_numbering.py` requires Python 3 and no extra packages. Run from
this repository, quoting paths that contain spaces:

```bash
# Preview only; no files are written.
python3 utilities/oni_continue_fov_numbering.py "/path/to/earlier" "/path/to/later"

# Apply the displayed renames.
python3 utilities/oni_continue_fov_numbering.py "/path/to/earlier" "/path/to/later" --apply
```

Finish acquisition before applying, and keep both folders unchanged during the
run. The arguments specify chronological order; dates are not inferred.

The script matches immediate subdirectories named `<condition>-FOV` and
`<condition>-FOV-N`, using the entire condition string. Bare `FOV` counts as zero.
For shared conditions, later FOVs are sorted numerically and numbered consecutively
from the earlier maximum plus one, closing later gaps. For example, earlier
`FOV` through `FOV-13` and later `FOV` through `FOV-25` produce later `FOV-14`
through `FOV-39`. Earlier folders and later-only conditions stay unchanged.

Matching identifiers in nested directory and file names are updated as well,
including ONI filenames with experiment prefixes. Matching respects identifier
boundaries (`FOV-1` does not match `FOV-10`). File contents, including metadata,
are never edited. Symlinks are skipped and never followed; links that reference
renamed paths are not repaired.

A condition whose later indices are all greater than its earlier maximum is
considered already numbered and skipped, even if there are gaps. This makes a
completed run a no-op when repeated. Duplicate indices (such as `FOV` and
`FOV-0`), conflicting destinations, and identical/nested inputs are rejected.
The whole mapping is checked before any rename. Descendants move before parents;
higher-index FOVs move first unless a destination must first be vacated.

### Journal and interrupted runs

Applying creates `.oni_continue_fov_numbering.jsonl` in the later folder. It records
the ordered plan as paths relative to that folder, an intent before each rename,
a completion record after each rename, and a final `complete` event. Paths in the
plan refer to the location at the time of each operation; a subsequent parent
rename moves previously renamed descendants too.

An incomplete or malformed journal blocks another apply, including a run with
no planned changes. There is no automatic resume or rollback. Inspect the journal
and filesystem, resolve the recorded operations manually, then archive or remove
the journal before applying again. An intent without a matching `renamed` event
requires checking both source and destination: interruption may have occurred
just before or after the rename. Keep the journal until recovery is complete.
A completed journal permits unchanged repeat runs; archive it before applying a
new nonempty plan to the same later folder.

Tests use temporary directories only:

```bash
python3 -m unittest discover -s utilities/tests -p 'test_oni*.py'
```
