from __future__ import annotations

from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from .params import SACDParams
from .reconstruction import reconstruct
from .tiffio import read_tiff_stack, write_tiff_image


@dataclass(frozen=True)
class BatchResult:
    input: Path
    output: Path
    status: str
    input_shape: tuple[int, ...] | None = None
    output_shape: tuple[int, ...] | None = None
    wavelength_nm: float | None = None
    frames_per_sacd: int | None = None


def find_input_files(input_path: str | Path, glob_pattern: str = "*frames_1-50.tif") -> list[Path]:
    path = Path(input_path)
    if path.is_file():
        return [path]
    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {path}")
    return sorted(f for f in path.glob(glob_pattern) if f.is_file())


def sacdpy_output_path(
    input_file: str | Path,
    output_dir: str | Path | None = None,
    frame_suffix: str = "_frames_1-50",
    output_label: str = "SACDpy",
    channels: tuple[str, ...] = ("left", "right"),
) -> Path:
    path = Path(input_file)
    stem = path.stem

    for channel in channels:
        suffix = f"-{channel}{frame_suffix}"
        if stem.endswith(suffix):
            output_name = f"{stem[: -len(suffix)]}-{output_label}-{channel}.tif"
            break
    else:
        output_name = f"{stem.removesuffix(frame_suffix)}-{output_label}.tif"

    parent = path.parent if output_dir is None else Path(output_dir)
    return parent / output_name


def wavelength_for_file(
    input_file: str | Path,
    default_wavelength_nm: float,
    wavelength_by_name: Mapping[str, float] | None = None,
) -> float:
    name = Path(input_file).name.lower()
    for token, wavelength in (wavelength_by_name or {}).items():
        if token.lower() in name:
            return float(wavelength)
    return float(default_wavelength_nm)


def params_for_file(
    input_file: str | Path,
    *,
    pixel_nm: float,
    na: float,
    default_wavelength_nm: float,
    wavelength_by_name: Mapping[str, float] | None = None,
    mag: int = 2,
    iter1: int = 7,
    iter2: int = 8,
    ac_order: int = 2,
    subfactor: float = 0.8,
    frames_per_sacd: int | None = None,
    ifregistration: bool = False,
    ifbackground: bool = False,
    backgroundfactor: float = 2.0,
    ifsparsedecon: bool = False,
    fidelity: float = 100.0,
    tcontinuity: float = 0.1,
    sparsity: float = 1.0,
    sparse_iterations: int = 100,
) -> SACDParams:
    return SACDParams(
        pixel_nm=pixel_nm,
        wavelength_nm=wavelength_for_file(input_file, default_wavelength_nm, wavelength_by_name),
        na=na,
        mag=mag,
        iter1=iter1,
        iter2=iter2,
        ac_order=ac_order,
        subfactor=subfactor,
        frames_per_sacd=frames_per_sacd,
        ifregistration=ifregistration,
        ifbackground=ifbackground,
        backgroundfactor=backgroundfactor,
        ifsparsedecon=ifsparsedecon,
        fidelity=fidelity,
        tcontinuity=tcontinuity,
        sparsity=sparsity,
        sparse_iterations=sparse_iterations,
    )


def run_batch_reconstruction(
    input_files: list[Path],
    *,
    output_dir: str | Path | None = None,
    overwrite_outputs: bool = False,
    pixel_nm: float,
    na: float,
    default_wavelength_nm: float,
    wavelength_by_name: Mapping[str, float] | None = None,
    mag: int = 2,
    iter1: int = 7,
    iter2: int = 8,
    ac_order: int = 2,
    subfactor: float = 0.8,
    frames_per_sacd: int | None = None,
    ifregistration: bool = False,
    ifbackground: bool = False,
    backgroundfactor: float = 2.0,
    ifsparsedecon: bool = False,
    fidelity: float = 100.0,
    tcontinuity: float = 0.1,
    sparsity: float = 1.0,
    sparse_iterations: int = 100,
    show_progress: bool = False,
) -> list[BatchResult]:
    results: list[BatchResult] = []

    progress_context = nullcontext()
    if show_progress:
        from rich.progress import (
            BarColumn,
            MofNCompleteColumn,
            Progress,
            SpinnerColumn,
            TextColumn,
            TimeElapsedColumn,
            TimeRemainingColumn,
        )

        progress_context = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
        )

    with progress_context as progress:
        task_id = None
        if progress is not None:
            task_id = progress.add_task("SACD batch", total=len(input_files))

        for input_file in input_files:
            output_file = sacdpy_output_path(input_file, output_dir=output_dir)
            output_file.parent.mkdir(parents=True, exist_ok=True)

            if progress is not None and task_id is not None:
                progress.update(task_id, description=f"Processing {input_file.name}")

            if output_file.exists() and not overwrite_outputs:
                results.append(BatchResult(input=input_file, output=output_file, status="skipped_existing"))
                if progress is not None and task_id is not None:
                    progress.advance(task_id)
                continue

            params = params_for_file(
                input_file,
                pixel_nm=pixel_nm,
                na=na,
                default_wavelength_nm=default_wavelength_nm,
                wavelength_by_name=wavelength_by_name,
                mag=mag,
                iter1=iter1,
                iter2=iter2,
                ac_order=ac_order,
                subfactor=subfactor,
                frames_per_sacd=frames_per_sacd,
                ifregistration=ifregistration,
                ifbackground=ifbackground,
                backgroundfactor=backgroundfactor,
                ifsparsedecon=ifsparsedecon,
                fidelity=fidelity,
                tcontinuity=tcontinuity,
                sparsity=sparsity,
                sparse_iterations=sparse_iterations,
            )
            stack = read_tiff_stack(input_file)
            sacd = reconstruct(stack, params)
            write_tiff_image(output_file, sacd)
            results.append(
                BatchResult(
                    input=input_file,
                    output=output_file,
                    status="written",
                    input_shape=stack.shape,
                    output_shape=sacd.shape,
                    wavelength_nm=params.wavelength_nm,
                    frames_per_sacd=params.frames_per_sacd,
                )
            )
            if progress is not None and task_id is not None:
                progress.advance(task_id)

    return results
