"""Simulated-device tests for the MCM alignment plan in startup/97-setup_plans.py.

The plan functions are extracted with `ast` (the startup file is not
importable: numeric prefix, injected globals, hardware at import) and exec'd
into a namespace of simulated devices. dscan, calculate_max_value and
flush_tiled are scripted stubs. No EPICS or hardware is involved.
"""

import ast
import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
import pytest
from bluesky import RunEngine
from bluesky.utils import RequestAbort, RequestStop
from ophyd import Signal
from ophyd.sim import SynAxis
from ophyd.status import Status

ROOT = Path(__file__).resolve().parent.parent
PLANS = ROOT / "startup" / "97-setup_plans.py"
_spec = importlib.util.spec_from_file_location(
    "mcm_alignment", ROOT / "startup" / "utils" / "mcm_alignment.py")
m = importlib.util.module_from_spec(_spec)
sys.modules["mcm_alignment"] = m
_spec.loader.exec_module(m)

XPOS = -0.941


class LogAxis(SynAxis):
    """Simulated motor that logs sets and can fail on chosen target values."""

    def __init__(self, *args, log=None, fail_values=(), **kwargs):
        super().__init__(*args, precision=8, **kwargs)
        self._log = log if log is not None else []
        self._fail_values = tuple(fail_values)

    def set(self, value):
        self._log.append((self.name, float(value)))
        if any(abs(float(value) - f) < 1e-9 for f in self._fail_values):
            st = Status()
            st.set_exception(RuntimeError(f"{self.name}: motor fault"))
            return st
        return super().set(value)


def _load_plan_namespace(ns):
    tree = ast.parse(PLANS.read_text(), filename=str(PLANS))
    keep = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name in ("AlignmentFailure", "OperatorAbort"):
            keep.append(node)
        elif isinstance(node, ast.FunctionDef) and (
                node.name.startswith("_mcm_") or node.name.startswith("mcm_setup")):
            keep.append(node)
    assert any(n.name == "mcm_setup" for n in keep)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(PLANS), "exec"), ns)


class Env:
    """Simulated beamline plus plan namespace."""

    def __init__(self, hrm="out", y=4.5, intensity=2e5, script=(), fail=None,
                 hrm_pos=None):
        self.log = []          # every motor set, in order
        self.events = []       # dscan / fit / flush order
        self.script = list(script)
        self.fit_calls = []
        self.flush_error = None
        fail = fail or {}

        def ax(name, value, **kw):
            return LogAxis(name=name, value=value, log=self.log,
                           fail_values=fail.get(name, ()), **kw)

        ux, dx = hrm_pos if hrm_pos else ((0.0, 0.0) if hrm == "in" else (-20.0, -20.0))
        self.d = SimpleNamespace(
            hrm2=SimpleNamespace(ux=ax("hrm2_ux", ux), dx=ax("hrm2_dx", dx),
                                 dch=ax("hrm2_dch", 0.0), dy=ax("hrm2_dy", 0.0)),
            anc_xtal=SimpleNamespace(y=ax("anc_xtal_y", y)),
            analyzer_slits=SimpleNamespace(
                top=ax("slit_top", 1.0), bottom=ax("slit_bottom", -1.0),
                outboard=ax("slit_outboard", 1.5), inboard=ax("slit_inboard", -1.5)),
            anpd=ax("anpd", -90.0), whl=ax("whl", 1.0),
            spec=SimpleNamespace(tth=ax("spec_tth", 30.0)),
            airpad=Signal(name="airpad", value=0),
            det2=SimpleNamespace(
                current1=SimpleNamespace(mean_value=Signal(
                    name="det2_current1_mean_value", value=intensity)),
                em_range=Signal(name="det2_em_range", value=7)),
            mcm=SimpleNamespace(x=ax("mcm_x", -0.9), y=ax("mcm_y", 0.1)),
            sample_stage=SimpleNamespace(tx=ax("s_tx", 0.0)),
        )
        env = self

        def dscan(mot, start, stop, steps, det, ct, det_ch=None):
            env.events.append(("dscan", mot.name))
            yield from bpp.run_wrapper(bps.null())
            return []

        def calculate_max_value(uid=-1, x="", y="", delta=1, sampling=200):
            env.events.append(("fit", uid))
            env.fit_calls.append(dict(uid=uid, x=x, y=y, delta=delta, sampling=sampling))
            item = env.script.pop(0)
            if isinstance(item, Exception):
                raise item
            return item, 1e6

        def flush_tiled(timeout=30):
            env.events.append(("flush",))
            if env.flush_error:
                raise env.flush_error

        ns = dict(
            bps=bps, bpp=bpp, RequestAbort=RequestAbort, RequestStop=RequestStop,
            dscan=dscan, calculate_max_value=calculate_max_value,
            flush_tiled=flush_tiled, **vars(self.d))
        for name in ("McmAlignConfig", "CorrectionLimitError", "HRM_IN", "HRM_OUT",
                     "CONVERGED", "CORRECT", "classify_hrm_state",
                     "hrm_measurement_settings", "decide_analyzer_move",
                     "analyzer_readback_ok", "tx_correction", "dch_moves",
                     "dy_correction", "convergence_step"):
            ns[name] = getattr(m, name)
        _load_plan_namespace(ns)
        self.ns = ns
        self.RE = RunEngine({}, call_returns_result=True)
        self.start_uids = []
        self.RE.subscribe(lambda n, d: self.start_uids.append(d["uid"]), "start")

    def run(self, ask=None, **kw):
        kw.setdefault("interactive", False)
        if ask is None:
            def ask(prompt):
                raise AssertionError(f"unexpected prompt: {prompt}")
        return self.RE(self.ns["mcm_setup"](ask=ask, **kw)).plan_result

    def moved(self, name):
        return [v for n, v in self.log if n == name]

    @property
    def pos(self):
        d = self.d
        return dict(
            y=d.anc_xtal.y.position, whl=d.whl.position, anpd=d.anpd.position,
            tth=d.spec.tth.position, em=d.det2.em_range.get(), air=d.airpad.get(),
            top=d.analyzer_slits.top.position, bottom=d.analyzer_slits.bottom.position,
            outboard=d.analyzer_slits.outboard.position,
            inboard=d.analyzer_slits.inboard.position,
            mx=d.mcm.x.position, my=d.mcm.y.position, tx=d.sample_stage.tx.position,
            dch=d.hrm2.dch.position, dy=d.hrm2.dy.position)


ENTRY = dict(whl=1.0, anpd=-90.0, top=1.0, bottom=-1.0, outboard=1.5, inboard=-1.5, em=7)
PREP = dict(y=0.5, anpd=0.0, top=2.0, bottom=-2.0, outboard=2.0, inboard=-2.0, tth=0.0, air=1)


# --- success paths ---------------------------------------------------------

def test_hrm_out_success_retains_measurement_configuration():
    e = Env("out", script=[XPOS + 0.01, XPOS + 0.0005, 0.123])
    r = e.run()
    assert r["status"] == "success"
    assert r["hrm_state"] == "out"
    assert r["axes"] == {"horizontal": "aligned", "vertical": "aligned"}
    p = e.pos
    for k, v in PREP.items():
        assert p[k] == pytest.approx(v, abs=1e-6), k
    assert (p["whl"], p["em"]) == (2.0, 0)           # HRM out
    assert p["tx"] == pytest.approx(0.01, abs=1e-6)   # tx correction = dxmax
    assert p["my"] == pytest.approx(0.123, abs=1e-6)  # moved to ymax
    assert p["mx"] == pytest.approx(XPOS, abs=1e-6)
    assert [c["motor"] for c in r["corrections"]] == ["s_tx"]
    assert r["preparation_restored"] is False


def test_hrm_out_success_restore_positions_true():
    e = Env("out", script=[XPOS + 0.0005, 0.2])
    r = e.run(restore_positions=True)
    assert r["status"] == "success"
    p = e.pos
    assert p["y"] == pytest.approx(4.5, abs=1e-6)
    for k, v in ENTRY.items():
        assert p[k] == pytest.approx(v, abs=1e-6), k
    assert p["tth"] == pytest.approx(0.0, abs=1e-6)  # never restored
    assert p["air"] == 1                              # never restored
    assert r["preparation_restored"] is True


def test_hrm_in_success_two_axes():
    e = Env("in", script=[XPOS - 0.04, XPOS + 0.001, 0.1 + 0.01, 0.1 + 0.0005])
    r = e.run()
    assert r["status"] == "success"
    p = e.pos
    assert (p["whl"], p["em"]) == (0.0, 1)           # HRM in
    assert p["dch"] == pytest.approx(0.004, abs=1e-6)  # -dxmax/10
    assert p["dy"] == pytest.approx(-0.05, abs=1e-6)   # -5 * dymax
    assert p["tx"] == 0.0                              # HRM-in uses no tower move
    assert p["my"] == pytest.approx(0.1, abs=1e-6)     # back at MCM_YPOS


def test_hrm_in_positive_dxmax_two_step_order():
    e = Env("in", script=[XPOS + 0.04, XPOS + 0.001], )
    r = e.run(align_vertical=False)
    assert r["status"] == "success"
    assert e.moved("hrm2_dch") == pytest.approx([-0.01, -0.01 + 0.005])
    assert e.pos["dch"] == pytest.approx(-0.005, abs=1e-6)


def test_already_lowered_analyzer_is_not_moved():
    e = Env("out", y=0.5, script=[XPOS, 0.1])
    r = e.run(restore_positions=True)
    assert r["status"] == "success"
    assert e.moved("anc_xtal_y") == []


def test_exact_uid_and_flush_before_fit():
    e = Env("out", script=[XPOS, 0.1])
    e.run()
    uids = [c["uid"] for c in e.fit_calls]
    assert uids == e.start_uids and all(isinstance(u, str) for u in uids)
    assert e.fit_calls[0]["sampling"] == 200 and e.fit_calls[0]["x"] == "mcm_x"
    for i, ev in enumerate(e.events):
        if ev[0] == "fit":
            assert e.events[i - 1] == ("flush",) and e.events[i - 2][0] == "dscan"


# --- preparation failures --------------------------------------------------

@pytest.mark.parametrize("pos", [(0.0, -20.0), (-20.0, 0.0), (-10.0, -10.0)])
def test_inconsistent_hrm_stops_before_any_move(pos):
    e = Env(hrm_pos=pos)
    r = e.run()
    assert r["status"] == "failed" and "inconsistent" in r["reason"]
    assert e.log == []


@pytest.mark.parametrize("y", [4.0, 5.0, 9.0])
def test_analyzer_target_out_of_interval_stops_before_any_move(y):
    e = Env(y=y)
    r = e.run()
    assert r["status"] == "failed" and f"{y:.4f}" in r["reason"]
    assert e.log == []


def test_low_intensity_fails_and_restores():
    e = Env("out", intensity=5e4)
    r = e.run()
    assert r["status"] == "failed" and "Low intensity" in r["reason"]
    p = e.pos
    assert p["y"] == pytest.approx(4.5, abs=1e-6)
    for k, v in ENTRY.items():
        assert p[k] == pytest.approx(v, abs=1e-6), k
    assert e.fit_calls == []


def test_prep_move_failure_aborts_and_restores():
    e = Env("out", fail={"spec_tth": [0.0]})
    r = e.run()
    assert r["status"] == "failed" and "motor fault" in r["reason"]
    assert e.pos["y"] == pytest.approx(4.5, abs=1e-6)
    assert e.pos["whl"] == pytest.approx(1.0, abs=1e-6)
    assert e.fit_calls == []


def test_restoration_error_does_not_hide_original_failure():
    # restoring whl to its entry value (1.0) fails
    e = Env("out", intensity=5e4, fail={"whl": [1.0]})
    r = e.run()
    assert r["status"] == "failed" and "Low intensity" in r["reason"]
    assert any("whl" in s for s in r["restore_errors"])
    assert r["preparation_restored"] is False
    assert e.pos["y"] == pytest.approx(4.5, abs=1e-6)  # other items still restored


# --- alignment failures ----------------------------------------------------

def test_non_convergence_after_four_scans():
    e = Env("out", script=[XPOS + 0.01] * 4)
    r = e.run()
    assert r["status"] == "failed" and "did not converge" in r["reason"]
    assert len(r["history"]) == 4
    assert len(r["corrections"]) == 3            # none after the fourth scan
    p = e.pos
    assert p["tx"] == pytest.approx(0.03, abs=1e-6)  # applied corrections kept
    assert p["mx"] == pytest.approx(XPOS, abs=1e-6)
    assert p["y"] == pytest.approx(4.5, abs=1e-6)    # failure restores


def test_fit_failure_stops_and_restores():
    e = Env("out", script=[ValueError("Delta value is greater than the lower limit")])
    r = e.run()
    assert r["status"] == "failed" and "Fit failed" in r["reason"]
    assert e.pos["mx"] == pytest.approx(XPOS, abs=1e-6)
    assert e.pos["y"] == pytest.approx(4.5, abs=1e-6)
    assert e.pos["tx"] == 0.0


def test_flush_failure_is_a_failure():
    e = Env("out", script=[XPOS])
    e.flush_error = RuntimeError("Tiled raw buffer was not flushed")
    r = e.run()
    assert r["status"] == "failed" and "not flushed" in r["reason"]


@pytest.mark.parametrize("hrm,err,motor", [("out", 0.0201, "s_tx"), ("in", 0.081, "hrm2_dch")])
def test_correction_over_limit_fails_without_moving(hrm, err, motor):
    e = Env(hrm, script=[XPOS + err])
    r = e.run()
    assert r["status"] == "failed" and "exceeds limit" in r["reason"]
    assert e.moved(motor) == []
    assert e.pos["y"] == pytest.approx(4.5, abs=1e-6)


def test_hrm_in_vertical_failure_returns_to_fixed_ypos_and_restores():
    e = Env("in", script=[XPOS + 0.001, 0.1 + 0.05])  # dy = -0.25, over 0.2
    r = e.run()
    assert r["status"] == "failed" and "hrm2.dy" in r["reason"] or "hrm2_dy" in r["reason"]
    assert e.pos["my"] == pytest.approx(0.1, abs=1e-6)
    assert e.pos["y"] == pytest.approx(4.5, abs=1e-6)


def test_hrm_out_vertical_failure_returns_to_pre_scan_position():
    e = Env("out", script=[XPOS, ValueError("edge")])
    e.d.mcm.y.set(0.37).wait()
    r = e.run()
    assert r["status"] == "failed"
    assert e.pos["my"] == pytest.approx(0.37, abs=1e-6)
    assert e.pos["y"] == pytest.approx(4.5, abs=1e-6)


# --- interaction, interruption, skipping ----------------------------------

def _scripted_ask(answers):
    it = iter(answers)
    prompts = []

    def ask(prompt):
        prompts.append(prompt)
        return next(it)
    ask.prompts = prompts
    return ask


def test_declined_correction_is_operator_abort_without_restoration():
    e = Env("out", script=[XPOS + 0.01])
    ask = _scripted_ask(["yes", "no"])  # horizontal axis yes, correction no
    r = e.run(ask=ask, interactive=True)
    assert r["status"] == "aborted" and "declined" in r["reason"]
    p = e.pos
    assert p["y"] == pytest.approx(0.5, abs=1e-6)     # stays lowered
    assert p["whl"] == pytest.approx(2.0, abs=1e-6)   # stays at measurement value
    assert p["tx"] == 0.0
    assert r["preparation_restored"] is False


def test_interactive_confirmation_shows_motor_position_displacement_target(capsys):
    e = Env("out", script=[XPOS + 0.01, XPOS, 0.2])
    ask = _scripted_ask(["yes", "yes", "yes", "yes"])  # h axis, tx, v axis, mcm.y
    r = e.run(ask=ask, interactive=True)
    assert r["status"] == "success"
    out = capsys.readouterr().out
    assert "Proposed: s_tx at 0" in out and "+0.01" in out and "target 0.01" in out
    assert "Proposed: mcm_y" in out


def test_two_step_correction_needs_one_confirmation():
    e = Env("in", script=[XPOS + 0.04, XPOS])
    ask = _scripted_ask(["yes", "yes"])  # h axis, one confirmation for both dch moves
    r = e.run(ask=ask, interactive=True, align_vertical=False)
    assert r["status"] == "success"
    assert len(ask.prompts) == 2
    assert len(e.moved("hrm2_dch")) == 2


def test_operator_interruption_does_not_restore():
    e = Env("out", script=[XPOS + 0.01])

    def ask(prompt):
        if "Perform" in prompt:
            return "yes"
        raise RequestAbort()  # operator aborts at the correction prompt

    # RunEngine handles RequestAbort itself and ends the run as 'abort'
    result = e.RE(e.ns["mcm_setup"](ask=ask, interactive=True))
    assert result.exit_status == "abort"
    assert e.pos["y"] == pytest.approx(0.5, abs=1e-6)   # no restoration
    assert e.pos["whl"] == pytest.approx(2.0, abs=1e-6)


def test_both_axes_skipped_by_operator_counts_as_success():
    e = Env("out")
    ask = _scripted_ask(["no", "no"])
    r = e.run(ask=ask, interactive=True)
    assert r["status"] == "success"
    assert r["axes"] == {"horizontal": "skipped", "vertical": "skipped"}
    assert e.fit_calls == []
    assert e.pos["y"] == pytest.approx(0.5, abs=1e-6)  # restore_positions False


def test_silent_mode_axis_selection_and_no_prompts():
    e = Env("out", script=[XPOS + 0.01, XPOS])
    r = e.run(align_vertical=False)  # default ask raises if a prompt appears
    assert r["status"] == "success"
    assert r["axes"] == {"horizontal": "aligned", "vertical": "skipped"}
    assert e.pos["tx"] == pytest.approx(0.01, abs=1e-6)


def test_skipped_axes_use_normal_success_restoration():
    e = Env("out")
    r = e.run(align_horizontal=False, align_vertical=False, restore_positions=True)
    assert r["status"] == "success" and r["preparation_restored"] is True
    assert e.pos["y"] == pytest.approx(4.5, abs=1e-6)


def test_unknown_option_rejected():
    e = Env("out")
    with pytest.raises(TypeError):
        e.RE(e.ns["mcm_setup"](not_an_option=1))
