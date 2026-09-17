# Parallel reconstruction and notebook progress validation

Date: 2026-09-16. Host: 10 CPU cores, 16 GiB RAM.

The three reconstruction batch pipelines default to two workers. Their shared
spawn pool bounds outstanding work to twice the worker count, caps numerical
library threads to one per worker, and settles running jobs after cancellation.
FOVs are sequential. Single-image batch and Dox remain serial and use Rich
progress, as do the three parallel notebooks.

## Real-FOV benchmark limitation

The representative 0915 FOV was `SHA-DAPI-SPEN_JFX650-200ms-FOV` (seven z planes,
two simultaneous camera channels). Each benchmark used a separate temporary
output directory and the actual TIFF validation/publication path. Parent plus
descendant RSS and system-wide swap were sampled every 0.1 seconds.

The benchmark guard stopped a run if swap usage grew by more than 1 GiB, or
available memory stayed below 512 MiB for five seconds. Both recorded runs
stopped on the swap-growth criterion:

| Workers | Elapsed before stop | Sampled peak process-tree RSS | Swap used at start | Peak swap used |
| --- | ---: | ---: | ---: | ---: |
| 2 | 9.07 s | 2,797,682,688 bytes | 7,441,285,120 bytes | 8,531,148,800 bytes |
| 3 | 10.95 s | 3,741,581,312 bytes | 8,423,342,080 bytes | 9,663,676,416 bytes |

Minimum available system memory was 2,713,026,560 and 2,689,122,304 bytes,
respectively. The monitor reported swap-in deltas of 5,832,704 and 21,200,896
bytes, and swap-out deltas of 2,555,904 and 2,801,664 bytes. Swap is system-wide,
so these changes cannot be attributed exclusively to SACD. RSS sums may also
count shared pages more than once and omit compressed/swapped-out memory.

Neither run completed; these are **not full-run runtimes or full-run memory
peaks**, and no two-versus-three speedup or real-FOV numerical comparison can
be claimed from them. After reviewing these results, the user selected two
workers as the default and requested no further testing. A clean
comparison requires a later run with less system memory pressure. No completed
dataset was modified.

## Automated acceptance checks

Final complete suite: **136 passed, 2 subtests passed** in 59.57 seconds.
This suite ran before the subsequent default-only change from three to two;
no additional tests or benchmarks were run after that change.
The 27 warnings are TIFF fixture and plotting-library deprecation warnings.
All five notebooks passed schema and Python-cell syntax validation; nbformat
also reported the existing missing-cell-ID warning.

Deterministic fixtures compare serial and three-worker arrays for multicolor
single-plane, sequential multicolor z-stack, simultaneous dual-view z-stack,
and timelapse z-stack, including saved MIPs and raw references where applicable.
Additional checks cover bounded submissions, completion order, worker errors,
interruption cleanup, resume after worker-count changes, failed-entry retries,
and progress accounting. All five notebook schemas and code cells are validated
without executing historical dataset processing cells.
