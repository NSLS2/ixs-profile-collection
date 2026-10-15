"""Offline tests for alignment scan diagnostic helpers."""

import ast
import io
from contextlib import redirect_stdout
from pathlib import Path
from time import monotonic
from types import SimpleNamespace

import bluesky.preprocessors as bpp
import bluesky.plan_stubs as bps
from bluesky import RunEngine
from bluesky.utils import Msg
from ophyd import Signal
from ophyd.sim import SynAxis


STARTUP_FILE = (
    Path(__file__).resolve().parent.parent / "startup" / "98-alignment.py"
)


def _load_helpers():
    """Load helpers without executing startup-time hardware and file access."""
    source = STARTUP_FILE.read_text()
    module = ast.parse(source, filename=str(STARTUP_FILE))
    names = {
        "_diagnostic_message_tracer",
        "_diagnostic_ah_preflight",
        "_capture_start_uid",
        "_ccr_temporary_settings_wrapper",
        "_confirm_tiled_run",
        "_check_bpm1_y_at_zero",
        "alignment_scan",
    }
    functions = [
        node
        for node in module.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    namespace = {
        "bpp": bpp,
        "bps": bps,
        "monotonic": monotonic,
        "warnings": __import__("warnings"),
    }
    helper_module = ast.Module(body=functions, type_ignores=[])
    exec(compile(helper_module, str(STARTUP_FILE), "exec"), namespace)
    return namespace


def test_message_tracer_preserves_messages_and_reports_wait_objects():
    helpers = _load_helpers()
    reports = []
    motor = SimpleNamespace(name="dcm_p1")

    def plan():
        yield Msg("set", motor, 12.5, group="move-group")
        yield Msg("wait", None, group="move-group")

    messages = list(
        helpers["_diagnostic_message_tracer"](plan(), reports.append)
    )

    assert [msg.command for msg in messages] == ["set", "wait"]
    assert messages[0].obj is motor
    assert "object=dcm_p1 target=12.5" in reports[0]
    assert "for=dcm_p1" in reports[1]


def test_ah_preflight_reports_each_synchronous_read():
    helpers = _load_helpers()
    reports = []

    class Signal:
        def __init__(self, pvname, value):
            self.pvname = pvname
            self.value = value

        def get(self):
            return self.value

    detector = SimpleNamespace(
        name="tm1",
        values_per_read=Signal("TM:ValuesPerRead", 8),
        averaging_time=Signal("TM:AveragingTime", 1.0),
    )
    helpers["_get_ah_scan_dets"] = lambda primary: [primary]

    helpers["_diagnostic_ah_preflight"](detector, reports.append)

    assert any("tm1.values_per_read (TM:ValuesPerRead)" in line for line in reports)
    assert any("tm1.values_per_read = 8" in line for line in reports)
    assert any("tm1.averaging_time = 1.0" in line for line in reports)


def test_alignment_diagnostics_are_opt_in():
    helpers = _load_helpers()

    class Signal:
        pvname = "TEST:PV"

        def get(self):
            return 1

    motor = SimpleNamespace(name="dcm_p1")
    detector = SimpleNamespace(
        name="tm1",
        values_per_read=Signal(),
        averaging_time=Signal(),
    )
    helpers["_SCAN_REGISTRY"] = {
        "dcm": {
            "motor": lambda: motor,
            "start": -1,
            "stop": 1,
            "steps": 3,
            "detector": lambda: detector,
            "ct": 1,
            "det_ch": [0],
            "instrument": "dcm",
            "step": "p1_scan",
            "det_label": "tm1",
            "mot_label": "dcm.p1",
            "checks": [],
        }
    }
    helpers["_get_ah_scan_dets"] = lambda primary: [primary]

    def fake_dscan(*args, **kwargs):
        uid = yield Msg("open_run", None)
        yield Msg("close_run", None)
        return []

    helpers["dscan"] = fake_dscan
    RE = RunEngine({}, call_returns_result=True)

    normal_output = io.StringIO()
    with redirect_stdout(normal_output):
        RE(helpers["alignment_scan"]("dcm"))

    diagnostic_output = io.StringIO()
    with redirect_stdout(diagnostic_output):
        RE(helpers["alignment_scan"]("dcm", diagnostics=True))

    assert "[alignment_scan:dcm" not in normal_output.getvalue()
    assert "plan entered" in diagnostic_output.getvalue()
    assert "AH501D read starting" in diagnostic_output.getvalue()
    assert "RE message: open_run" in diagnostic_output.getvalue()


def test_start_uid_callback_captures_exact_run():
    helpers = _load_helpers()
    holder = {}

    callback = helpers["_capture_start_uid"](holder)
    callback("start", {"uid": "exact-run-uid"})

    assert holder == {"uid": "exact-run-uid"}


def test_tiled_confirmation_retries_exact_uid():
    helpers = _load_helpers()
    reports = []

    class Catalog:
        def __init__(self):
            self.calls = []

        def __getitem__(self, uid):
            self.calls.append(uid)
            if len(self.calls) < 3:
                raise KeyError(uid)
            return SimpleNamespace(
                start={"uid": uid},
                metadata={"start": {"uid": uid}, "stop": {}},
            )

    catalog = Catalog()
    helpers["trc"] = catalog
    plan = helpers["_confirm_tiled_run"](
        "exact-run-uid", reports.append, attempts=3, delay=0
    )

    messages = list(plan)

    assert [msg.command for msg in messages] == ["sleep", "sleep"]
    assert catalog.calls == ["exact-run-uid"] * 3
    assert reports == ["Tiled run confirmed: exact-run-uid"]


def test_tiled_confirmation_warns_and_preserves_uid():
    helpers = _load_helpers()
    reports = []

    class Catalog:
        def __getitem__(self, uid):
            raise KeyError(uid)

    helpers["trc"] = Catalog()
    plan = helpers["_confirm_tiled_run"](
        "captured-run-uid", reports.append, attempts=2, delay=0
    )

    with __import__("warnings").catch_warnings(record=True) as caught:
        __import__("warnings").simplefilter("always")
        messages = list(plan)

    assert [msg.command for msg in messages] == ["sleep"]
    assert len(caught) == 1
    assert "captured-run-uid" in str(caught[0].message)


def test_alignment_records_exact_uid_confirmed_by_tiled():
    helpers = _load_helpers()
    motor = SimpleNamespace(name="dcm_p1")
    detector = SimpleNamespace(name="tm1")
    peak = SimpleNamespace(
        max=(1.25, 42.0),
        fwhm=0.5,
        cen=1.2,
        com=1.22,
    )
    recorded = {}

    helpers["_SCAN_REGISTRY"] = {
        "dcm": {
            "motor": lambda: motor,
            "start": -1,
            "stop": 1,
            "steps": 3,
            "detector": lambda: detector,
            "ct": 0.1,
            "det_ch": [0],
            "instrument": "dcm",
            "step": "p1_scan",
            "det_label": "tm1",
            "mot_label": "dcm.p1",
            "checks": [],
        }
    }

    def fake_dscan(*args, **kwargs):
        yield Msg("open_run", None)
        yield Msg("close_run", None)
        return [peak]

    class Catalog:
        def __init__(self):
            self.requested = []

        def __getitem__(self, uid):
            self.requested.append(uid)
            return SimpleNamespace(
                start={"uid": uid},
                metadata={"start": {"uid": uid}, "stop": {}},
            )

    def fake_record_alignment(store, **kwargs):
        recorded.update(kwargs)
        yield from bps.null()

    catalog = Catalog()
    helpers.update(
        {
            "dscan": fake_dscan,
            "trc": catalog,
            "store": object(),
            "record_alignment": fake_record_alignment,
            "_read_crl_state": lambda: "out",
            "_read_hrm_state": lambda: "out",
            "_read_ring_current": lambda: 500.0,
        }
    )

    RE = RunEngine({}, call_returns_result=True)
    result = RE(helpers["alignment_scan"]("dcm"))
    uid = result.run_start_uids[0]

    assert catalog.requested == [uid]
    assert recorded["scan_id"] == uid
    assert recorded["max_intensity"] == 42.0
    assert recorded["x_at_max"] == 1.25


def test_alignment_production_has_no_legacy_or_latest_lookup():
    source = STARTUP_FILE.read_text()

    assert "db[" not in source
    assert "trc[-1]" not in source


def test_bpm1_foil_precondition_uses_bpm1_diag_and_passes_near_zero():
    helpers = _load_helpers()
    helpers["bpm1_diag"] = SimpleNamespace(position=0.05)

    ok, message = helpers["_check_bpm1_y_at_zero"]()

    assert ok is True
    assert message == ""


def test_bpm1_foil_precondition_rejects_inserted_foil_position():
    helpers = _load_helpers()
    helpers["bpm1_diag"] = SimpleNamespace(position=25.0003)

    ok, message = helpers["_check_bpm1_y_at_zero"]()

    assert ok is False
    assert "bpm1_diag = 25.000 mm" in message
    assert "BPM1 Y-foil" in message


def _install_ccr_devices(helpers, anpd_position, em_range):
    motor = SynAxis(name="anpd")
    motor.set(anpd_position).wait()
    range_signal = Signal(name="det2_em_range", value=em_range)
    helpers["anpd"] = motor
    helpers["det2"] = SimpleNamespace(em_range=range_signal)
    return motor, range_signal


def test_ccr_temporary_settings_are_applied_and_restored():
    helpers = _load_helpers()
    motor, range_signal = _install_ccr_devices(helpers, -90.0, 2)
    values_during_scan = []

    def scan_plan():
        values_during_scan.append((motor.position, range_signal.get()))
        yield from bps.null()

    wrapped = helpers["_ccr_temporary_settings_wrapper"](
        scan_plan(), lambda message: None
    )
    RunEngine({})(wrapped)

    assert values_during_scan == [(40.0, 0)]
    assert motor.position == -90.0
    assert range_signal.get() == 2


def test_ccr_temporary_settings_skip_values_already_at_target():
    helpers = _load_helpers()
    motor, range_signal = _install_ccr_devices(helpers, 40.05, 0)
    set_messages = []

    def scan_plan():
        yield from bps.null()

    wrapped = helpers["_ccr_temporary_settings_wrapper"](
        scan_plan(), lambda message: None
    )

    def capture_set(msg):
        if msg.command == "set":
            set_messages.append(msg)
        return msg

    wrapped = bpp.msg_mutator(
        wrapped,
        capture_set,
    )
    RunEngine({})(wrapped)

    assert set_messages == []
    assert motor.position == 40.05
    assert range_signal.get() == 0


def test_ccr_temporary_settings_change_only_value_outside_target():
    helpers = _load_helpers()
    motor, range_signal = _install_ccr_devices(helpers, 40.0, 3)
    set_messages = []

    def scan_plan():
        assert motor.position == 40.0
        assert range_signal.get() == 0
        yield from bps.null()

    wrapped = helpers["_ccr_temporary_settings_wrapper"](
        scan_plan(), lambda message: None
    )

    def capture_set(msg):
        if msg.command == "set":
            set_messages.append((msg.obj.name, msg.args[0]))
        return msg

    RunEngine({})(bpp.msg_mutator(wrapped, capture_set))

    assert set_messages == [("det2_em_range", 0), ("det2_em_range", 3)]
    assert motor.position == 40.0
    assert range_signal.get() == 3


def test_ccr_temporary_settings_restore_after_scan_failure():
    helpers = _load_helpers()
    motor, range_signal = _install_ccr_devices(helpers, -145.0, 1)
    values_during_scan = []

    def failing_scan_plan():
        values_during_scan.append((motor.position, range_signal.get()))
        yield from bps.null()
        raise RuntimeError("synthetic scan failure")

    wrapped = helpers["_ccr_temporary_settings_wrapper"](
        failing_scan_plan(), lambda message: None
    )
    try:
        RunEngine({})(wrapped)
    except RuntimeError as exc:
        assert str(exc) == "synthetic scan failure"
    else:
        raise AssertionError("The synthetic scan failure was not raised")

    assert values_during_scan == [(40.0, 0)]
    assert motor.position == -145.0
    assert range_signal.get() == 1


def test_ccr_failed_precondition_does_not_apply_temporary_settings():
    helpers = _load_helpers()
    motor, range_signal = _install_ccr_devices(helpers, -90.0, 2)
    scan_called = []
    detector = SimpleNamespace(name="det2")

    def fake_dscan(*args, **kwargs):
        scan_called.append(True)
        yield from bps.null()

    helpers["_SCAN_REGISTRY"] = {
        "ccr": {
            "motor": lambda: SimpleNamespace(name="analyzer_cfth"),
            "start": -150,
            "stop": 150,
            "steps": 31,
            "detector": lambda: detector,
            "ct": 1,
            "det_ch": [0],
            "instrument": "ccr",
            "step": "cfth_scan",
            "det_label": "det2",
            "mot_label": "analyzer.cfth",
            "checks": [lambda: (False, "synthetic precondition failure")],
        }
    }
    helpers["dscan"] = fake_dscan

    RunEngine({})(helpers["alignment_scan"]("ccr"))

    assert scan_called == []
    assert motor.position == -90.0
    assert range_signal.get() == 2
