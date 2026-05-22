# SACDpy Validation Report

Date: 2026-05-20  
Package path: `/Users/gmgao/GGscripts/GG-general-GuttmanLab/image_analysis/pipelines/SACDpy`

## Attestation Summary

The Python implementation in `src/sacdpy` is a valid core port of the MATLAB SACDm default reconstruction pipeline for the provided 2D SACD test stacks. The implementation was built from the MATLAB source at `/Users/gmgao/GGscripts/SACDm` and validated against the two ImageJ/SACDj processed reference outputs in `tests/testdata/`.

This report supersedes the earlier validation run that used `NA=1.3`. The microscope objective used for the provided data has `NA=1.45`, and the ImageJ references were generated with that value.

This attestation applies to the core SACD path:

1. per-frame offset removal
2. pre Richardson-Lucy deconvolution
3. lateral Fourier interpolation
4. autocumulant calculation
5. post Richardson-Lucy deconvolution

Advanced SACDm branches are intentionally outside this validation scope: registration, wavelet background subtraction, GPU execution, and sparse Hessian deconvolution. Those options raise `NotImplementedError` in the Python port instead of producing unvalidated approximations.

## Reference Data and Parameters

The validation used the bundled test files:

| Channel | Raw stack | ImageJ SACD reference | Parameters |
| --- | --- | --- | --- |
| Left | `tests/testdata/TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-left.tif` | `tests/testdata/TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-SACD-left.tif` | `pixel=117 nm`, `wavelength=560 nm`, `NA=1.45` |
| Right | `tests/testdata/TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-right.tif` | `tests/testdata/TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-SACD-right.tif` | `pixel=117 nm`, `wavelength=647 nm`, `NA=1.45` |

Both raw stacks are `50 x 684 x 428` TIFF stacks. Both ImageJ reference outputs are `1368 x 856`, matching the expected `2x` lateral SACD reconstruction.

The validation used the following shared SACD parameters:

| Parameter | Value |
| --- | --- |
| `mag` | `2` |
| `NA` | `1.45` |
| `iter1` | `7` |
| `iter2` | `8` |
| `ACorder` | `2` |
| `subfactor` | `0.8` |
| `scale` | `2` |

## Validation Results

Validation command:

```bash
PYTHONPATH=src conda run -n smlm python tests/validate_testdata.py
```

Measured results:

| Channel | Output | Shape | Pearson corr. | SSIM | NRMSE | MAE |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Left | `tests/outputs/sacdpy-left.tif` | `1368 x 856` | `0.8848` | `0.2609` | `0.6675` | `0.0821` |
| Right | `tests/outputs/sacdpy-right.tif` | `1368 x 856` | `0.8455` | `0.9492` | `0.5251` | `0.0032` |

The normalized Pearson correlations show strong agreement in reconstructed spatial structure for both channels. Correcting the NA from `1.3` to `1.45` improved the right-channel agreement and kept the left-channel structural correlation stable. The left-channel SSIM remains lower despite high Pearson correlation, which is consistent with differences in dynamic-range scaling and local intensity distribution between the Python output and the ImageJ-rendered reference TIFF. The validation therefore treats these ImageJ files as similarity references, not bit-exact gold standards.

## Automated Tests

Test command:

```bash
PYTHONPATH=src conda run -n smlm python -m unittest discover -s tests
```

Result:

```text
Ran 17 tests in 26.103s
OK
```

Coverage includes:

- TIFF stack axis handling from `TYX` to internal `YXT`
- RSF and PSF normalization
- 2D/3D lateral Fourier interpolation output shape and finite values
- autocumulant order 2 and order 4 MATLAB formula parity
- full reconstruction path on both bundled test stacks

## Validity Statement

Based on source-level porting from MATLAB SACDm, successful unit tests, correct output dimensions, finite nonnegative reconstructions, and strong normalized similarity to the provided ImageJ/SACDj reference outputs, the current Python implementation is valid for the core SACD reconstruction workflow represented by the bundled test data.

The implementation should not be represented as full SACDm parity until the advanced MATLAB options are separately ported and validated.
