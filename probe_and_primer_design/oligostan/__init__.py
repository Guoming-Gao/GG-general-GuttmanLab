"""Python Oligostan port; see PARITY_PLAN.md before treating output as R-equivalent."""

from .config import DEFAULT_SETTINGS, FLAP_SEQUENCES
from .oligostan_core import get_probes_from_rna_dg37, process_probes_for_output

__all__ = [
    "DEFAULT_SETTINGS",
    "FLAP_SEQUENCES",
    "get_probes_from_rna_dg37",
    "process_probes_for_output",
]
