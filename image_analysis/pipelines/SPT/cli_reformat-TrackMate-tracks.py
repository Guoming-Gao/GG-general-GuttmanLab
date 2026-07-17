"""Run the preserved historical TrackMate reformatter."""

from pathlib import Path
import runpy

runpy.run_path(str(Path(__file__).parent / "legacy" / "cli_reformat-TrackMate-tracks.py"), run_name="__main__")
