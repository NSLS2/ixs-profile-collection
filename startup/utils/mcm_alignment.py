"""Hardware-free logic for the MCM alignment plan.

Implements the decisions described in docs/mcm_alignment_scope.md:
HRM-state classification, analyzer-height preparation, HRM-dependent
settings, correction calculations with limit checks, and convergence.

Nothing here touches devices, EPICS or the RunEngine, so it can be tested
offline. Corrections that exceed a limit raise CorrectionLimitError; they are
never clipped.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

_EPS = 1e-12  # guards limit comparisons against float rounding


class CorrectionLimitError(ValueError):
    """A correction exceeds a configured limit or motor travel limit."""


@dataclass
class McmAlignConfig:
    """Configurable parameters. Defaults are the confirmed scope values."""

    restore_positions: bool = False
    interactive: bool = True
    align_horizontal: bool = True
    align_vertical: bool = True

    mcm_xpos: float = -0.941  # mm, fixed horizontal reference

    # analyzer height
    analyzer_lowering: float = -4.0  # mm, relative
    analyzer_min: float = 0.0  # allowed interval, strictly between
    analyzer_max: float = 1.0

    # scans (either axis)
    scan_half_range: float = 0.2  # mm
    scan_points: int = 41
    scan_time: float = 1.0  # s per point

    # convergence
    horizontal_tolerance: float = 0.002  # mm, abs error <
    vertical_tolerance: float = 0.002  # mm, abs error <
    max_scans: int = 4

    low_intensity_threshold: float = 1e5

    # calculate_max_value arguments (resampling improves resolution only)
    fit_delta: int = 1
    fit_sampling: int = 200

    # correction limits (magnitudes)
    max_tx_correction: float = 0.02  # mm
    max_dy_correction: float = 0.2  # mm
    max_angular_term: float = 0.008  # degrees

    # correction factors and two-step offsets
    angular_divisor: float = 10.0  # angular_term = dxmax / divisor
    dy_factor: float = -5.0  # dy correction = factor * dymax
    dch_first_step: float = -0.01  # degrees
    dch_second_offset: float = 0.009  # degrees

    # HRM state: both axes within `hrm_tolerance` of the position
    hrm_in_position: float = 0.0  # mm
    hrm_out_position: float = -20.0  # mm
    hrm_tolerance: float = 2.0  # mm

    # preparation targets (mm unless noted)
    slit_top: float = 2.0
    slit_bottom: float = -2.0
    slit_outboard: float = 2.0
    slit_inboard: float = -2.0
    anpd_position: float = 0.0
    airpad_on: int = 1
    spec_tth_position: float = 0.0

    # HRM-dependent measurement settings
    whl_hrm_in: float = 0
    whl_hrm_out: float = 2
    em_range_hrm_in: int = 1
    em_range_hrm_out: int = 0


# ---------------------------------------------------------------------------
# HRM state
# ---------------------------------------------------------------------------

HRM_IN = "in"
HRM_OUT = "out"
HRM_INCONSISTENT = "inconsistent"


def classify_hrm_state(ux, dx, cfg: Optional[McmAlignConfig] = None) -> str:
    """Return 'in', 'out' or 'inconsistent' from hrm2.ux and hrm2.dx.

    in  : both axes within in_position +/- tolerance (0 +/- 2 mm)
    out : both axes within out_position +/- tolerance (-20 +/- 2 mm)
    anything else, including an unreadable value, is inconsistent.
    """
    cfg = cfg or McmAlignConfig()
    try:
        ux = float(ux)
        dx = float(dx)
    except (TypeError, ValueError):
        return HRM_INCONSISTENT
    if ux != ux or dx != dx:  # NaN
        return HRM_INCONSISTENT

    def near(value, center):
        return abs(value - center) <= cfg.hrm_tolerance

    if near(ux, cfg.hrm_in_position) and near(dx, cfg.hrm_in_position):
        return HRM_IN
    if near(ux, cfg.hrm_out_position) and near(dx, cfg.hrm_out_position):
        return HRM_OUT
    return HRM_INCONSISTENT


def hrm_measurement_settings(hrm_state: str, cfg: Optional[McmAlignConfig] = None) -> dict:
    """Return {'whl': ..., 'em_range': ...} for an 'in' or 'out' HRM state."""
    cfg = cfg or McmAlignConfig()
    if hrm_state == HRM_IN:
        return {"whl": cfg.whl_hrm_in, "em_range": cfg.em_range_hrm_in}
    if hrm_state == HRM_OUT:
        return {"whl": cfg.whl_hrm_out, "em_range": cfg.em_range_hrm_out}
    raise ValueError(f"No measurement settings for HRM state {hrm_state!r}.")


# ---------------------------------------------------------------------------
# Analyzer height
# ---------------------------------------------------------------------------

@dataclass
class AnalyzerDecision:
    ok: bool
    move: bool  # True if the analyzer must be moved
    y_initial: float
    y_target: Optional[float]  # None when no move is needed
    message: str


def _strictly_inside(value, cfg: McmAlignConfig) -> bool:
    return cfg.analyzer_min < value < cfg.analyzer_max


def decide_analyzer_move(y_initial, cfg: Optional[McmAlignConfig] = None) -> AnalyzerDecision:
    """Decide how to prepare the analyzer height.

    If 0 < y_initial < 1 mm the analyzer is already lowered and is left
    unchanged. Otherwise y_target = y_initial - 4 mm must satisfy
    0 < y_target < 1 mm (both boundaries rejected); if not, preparation stops
    and the message reports the initial and calculated positions.
    """
    cfg = cfg or McmAlignConfig()
    y_initial = float(y_initial)
    if _strictly_inside(y_initial, cfg):
        return AnalyzerDecision(
            True, False, y_initial, None,
            f"Analyzer already lowered (y = {y_initial:.4f} mm); unchanged.",
        )
    y_target = y_initial + cfg.analyzer_lowering
    if _strictly_inside(y_target, cfg):
        return AnalyzerDecision(
            True, True, y_initial, y_target,
            f"Lower analyzer from {y_initial:.4f} mm to {y_target:.4f} mm.",
        )
    return AnalyzerDecision(
        False, False, y_initial, y_target,
        f"Analyzer target outside ({cfg.analyzer_min}, {cfg.analyzer_max}) mm: "
        f"initial = {y_initial:.4f} mm, calculated target = {y_target:.4f} mm. "
        "Preparation stopped.",
    )


def analyzer_readback_ok(y_readback, cfg: Optional[McmAlignConfig] = None) -> bool:
    """True if the analyzer readback is strictly inside the allowed interval."""
    cfg = cfg or McmAlignConfig()
    try:
        return _strictly_inside(float(y_readback), cfg)
    except (TypeError, ValueError):
        return False


# ---------------------------------------------------------------------------
# Corrections
# ---------------------------------------------------------------------------

def _check_limit(value, limit, what):
    if abs(value) > limit + _EPS:
        raise CorrectionLimitError(
            f"{what} correction {value:+.6g} exceeds limit {limit:g}; not applied."
        )


def _check_travel(target, travel, what):
    """travel = (low, high) absolute motor limits, or None to skip."""
    if travel is None:
        return
    low, high = travel
    if not (low <= target <= high):
        raise CorrectionLimitError(
            f"{what} target {target:+.6g} outside motor limits [{low:g}, {high:g}]."
        )


def tx_correction(dxmax, cfg: Optional[McmAlignConfig] = None,
                  position=None, travel=None) -> float:
    """Relative sample_stage.tx move (mm) for HRM out: dxmax.

    position and travel (low, high) are optional; if both are given the
    absolute target is checked against the motor limits.
    """
    cfg = cfg or McmAlignConfig()
    move = float(dxmax)
    _check_limit(move, cfg.max_tx_correction, "sample_stage.tx")
    if position is not None:
        _check_travel(position + move, travel, "sample_stage.tx")
    return move


def dch_moves(dxmax, cfg: Optional[McmAlignConfig] = None,
              position=None, travel=None) -> list:
    """Relative hrm2.dch moves (degrees) for HRM in.

    dxmax < 0 : [-dxmax / 10]
    dxmax > 0 : [-0.01, 0.009 - dxmax / 10]
    dxmax == 0: []
    The limit applies to abs(dxmax / 10) only. Both absolute targets are
    validated against travel (if given) before any move is returned.
    """
    cfg = cfg or McmAlignConfig()
    dxmax = float(dxmax)
    term = dxmax / cfg.angular_divisor
    _check_limit(term, cfg.max_angular_term, "hrm2.dch angular term")
    if dxmax < 0:
        moves = [-term]
    elif dxmax > 0:
        moves = [cfg.dch_first_step, cfg.dch_second_offset - term]
    else:
        moves = []
    if position is not None:
        pos = position
        for step in moves:
            pos += step
            _check_travel(pos, travel, "hrm2.dch")
    return moves


def dy_correction(dymax, cfg: Optional[McmAlignConfig] = None,
                  position=None, travel=None) -> float:
    """Relative hrm2.dy move (mm) for HRM in: dy_factor * dymax (-5 * dymax)."""
    cfg = cfg or McmAlignConfig()
    move = cfg.dy_factor * float(dymax)
    _check_limit(move, cfg.max_dy_correction, "hrm2.dy")
    if position is not None:
        _check_travel(position + move, travel, "hrm2.dy")
    return move


# ---------------------------------------------------------------------------
# Convergence
# ---------------------------------------------------------------------------

CONVERGED = "converged"
CORRECT = "correct"  # apply a correction, then rescan
NOT_CONVERGED = "not_converged"  # final scan outside tolerance


def convergence_step(scan_number, error, tolerance, max_scans=4) -> str:
    """Decide what to do after scan number `scan_number` (1-based).

    abs(error) < tolerance           -> CONVERGED (no correction applied)
    outside, scan_number < max_scans -> CORRECT (apply correction, rescan)
    outside, scan_number >= max_scans-> NOT_CONVERGED (no further correction)
    """
    if abs(error) < tolerance:
        return CONVERGED
    if scan_number < max_scans:
        return CORRECT
    return NOT_CONVERGED
