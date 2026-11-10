"""Offline tests for startup/utils/mcm_alignment.py (no hardware)."""

import importlib.util
from pathlib import Path

import pytest

_PATH = (
    Path(__file__).resolve().parent.parent
    / "startup" / "utils" / "mcm_alignment.py"
)
_spec = importlib.util.spec_from_file_location("mcm_alignment", _PATH)
m = importlib.util.module_from_spec(_spec)
import sys
sys.modules["mcm_alignment"] = m  # dataclasses needs the module registered
_spec.loader.exec_module(m)


# --- HRM state -------------------------------------------------------------

@pytest.mark.parametrize("ux,dx,state", [
    (0, 0, "in"),
    (2, -2, "in"),
    (-20, -20, "out"),
    (-18, -22, "out"),
    (0, -20, "inconsistent"),
    (-20, 0, "inconsistent"),
    (-10, -10, "inconsistent"),
    (2.01, 0, "inconsistent"),
    (-17.9, -20, "inconsistent"),
    (None, 0, "inconsistent"),
    ("bad", 0, "inconsistent"),
    (float("nan"), 0, "inconsistent"),
])
def test_classify_hrm_state(ux, dx, state):
    assert m.classify_hrm_state(ux, dx) == state


def test_hrm_settings():
    assert m.hrm_measurement_settings("in") == {"whl": 0, "em_range": 1}
    assert m.hrm_measurement_settings("out") == {"whl": 2, "em_range": 0}
    with pytest.raises(ValueError):
        m.hrm_measurement_settings("inconsistent")


# --- analyzer height -------------------------------------------------------

def test_analyzer_already_lowered_unchanged():
    d = m.decide_analyzer_move(0.5)
    assert d.ok and not d.move and d.y_target is None


def test_analyzer_lowering_in_range():
    d = m.decide_analyzer_move(4.5)
    assert d.ok and d.move and d.y_target == pytest.approx(0.5)


@pytest.mark.parametrize("y", [4.0, 5.0])
def test_analyzer_target_at_boundary_rejected(y):
    d = m.decide_analyzer_move(y)
    assert not d.ok and not d.move
    assert f"{y:.4f}" in d.message and "calculated target" in d.message


@pytest.mark.parametrize("y", [0.0, 1.0])
def test_analyzer_initial_on_boundary_is_not_lowered(y):
    # 0 and 1 are not strictly inside, so the -4 mm target is evaluated
    d = m.decide_analyzer_move(y)
    assert not d.ok


@pytest.mark.parametrize("y", [-3.0, 3.0, 10.0])
def test_analyzer_target_out_of_range_rejected(y):
    assert not m.decide_analyzer_move(y).ok


def test_analyzer_readback():
    assert m.analyzer_readback_ok(0.5)
    for bad in (0.0, 1.0, -0.1, 1.2, None, "x"):
        assert not m.analyzer_readback_ok(bad)


# --- corrections -----------------------------------------------------------

def test_tx_correction_sign_and_limit():
    assert m.tx_correction(0.01) == pytest.approx(0.01)
    assert m.tx_correction(-0.02) == pytest.approx(-0.02)  # at limit passes
    with pytest.raises(m.CorrectionLimitError):
        m.tx_correction(0.0201)


def test_tx_travel_limit():
    with pytest.raises(m.CorrectionLimitError):
        m.tx_correction(0.01, position=0.995, travel=(-1, 1))


def test_dch_negative_dxmax_single_move():
    assert m.dch_moves(-0.04) == pytest.approx([0.004])


def test_dch_positive_dxmax_two_steps():
    assert m.dch_moves(0.04) == pytest.approx([-0.01, 0.009 - 0.004])


def test_dch_zero_no_moves():
    assert m.dch_moves(0.0) == []


def test_dch_limit_applies_to_term_only():
    # dxmax = 0.08 -> term 0.008 (at limit): allowed
    assert m.dch_moves(0.08) == pytest.approx([-0.01, 0.001])
    assert m.dch_moves(-0.08) == pytest.approx([0.008])
    with pytest.raises(m.CorrectionLimitError):
        m.dch_moves(0.081)
    with pytest.raises(m.CorrectionLimitError):
        m.dch_moves(-0.081)


def test_dch_validates_both_targets_before_returning():
    # first target (-0.01) is fine, second (-0.01 + 0.005) exceeds high limit
    with pytest.raises(m.CorrectionLimitError):
        m.dch_moves(0.04, position=0.0, travel=(-0.005, 0.0))
    # first target outside
    with pytest.raises(m.CorrectionLimitError):
        m.dch_moves(0.04, position=0.0, travel=(-0.001, 1.0))
    assert m.dch_moves(0.04, position=0.0, travel=(-1, 1))


def test_dy_correction_sign_and_limit():
    assert m.dy_correction(0.01) == pytest.approx(-0.05)
    assert m.dy_correction(-0.04) == pytest.approx(0.2)  # at limit passes
    with pytest.raises(m.CorrectionLimitError):
        m.dy_correction(0.041)


def test_configurable_factors():
    cfg = m.McmAlignConfig(dy_factor=-2.0, max_dy_correction=0.1)
    assert m.dy_correction(0.01, cfg) == pytest.approx(-0.02)
    with pytest.raises(m.CorrectionLimitError):
        m.dy_correction(0.06, cfg)


# --- convergence -----------------------------------------------------------

def test_converged_when_inside_tolerance():
    assert m.convergence_step(1, 0.0019, 0.002) == m.CONVERGED
    assert m.convergence_step(1, -0.0019, 0.002) == m.CONVERGED


def test_tolerance_is_strict():
    assert m.convergence_step(1, 0.002, 0.002) == m.CORRECT


@pytest.mark.parametrize("n", [1, 2, 3])
def test_correct_before_last_scan(n):
    assert m.convergence_step(n, 0.01, 0.002) == m.CORRECT


def test_fourth_scan_outside_is_not_converged():
    assert m.convergence_step(4, 0.01, 0.002) == m.NOT_CONVERGED


def test_fourth_scan_inside_converges():
    assert m.convergence_step(4, 0.001, 0.002) == m.CONVERGED


# --- defaults match scope --------------------------------------------------

def test_defaults():
    c = m.McmAlignConfig()
    assert (c.restore_positions, c.interactive) == (False, True)
    assert c.mcm_xpos == -0.941
    assert (c.scan_half_range, c.scan_points, c.scan_time) == (0.2, 41, 1.0)
    assert c.max_scans == 4 and c.low_intensity_threshold == 1e5
    assert (c.max_tx_correction, c.max_dy_correction, c.max_angular_term) == (0.02, 0.2, 0.008)
    assert not hasattr(c, "min_contrast")
    assert (c.fit_delta, c.fit_sampling) == (1, 200)
