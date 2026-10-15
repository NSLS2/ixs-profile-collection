"""Offline tests for electrometer settings in custom scans."""

import ast
from pathlib import Path
from types import SimpleNamespace

import bluesky.plan_stubs as bps


ROOT = Path(__file__).resolve().parent.parent
CUSTOM_SCANS_FILE = ROOT / "startup" / "97-custom_scans.py"
ALIGNMENT_FILE = ROOT / "startup" / "98-alignment.py"


class _Signal:
    def __init__(self, name, value):
        self.name = name
        self.parent = None
        self._value = value

    def get(self):
        return self._value


def _load_custom_scan_functions():
    source = CUSTOM_SCANS_FILE.read_text()
    module = ast.parse(source, filename=str(CUSTOM_SCANS_FILE))
    names = {"_scan_values_per_read", "dscan", "ascan"}
    functions = [
        node
        for node in module.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    namespace = {"bps": bps}
    exec(
        compile(
            ast.Module(body=functions, type_ignores=[]),
            str(CUSTOM_SCANS_FILE),
            "exec",
        ),
        namespace,
    )
    return namespace


def _detector(name, values_per_read, averaging_time):
    return SimpleNamespace(
        name=name,
        values_per_read=_Signal(f"{name}_values_per_read", values_per_read),
        averaging_time=_Signal(f"{name}_averaging_time", averaging_time),
    )


def _configuration_messages(scan_name):
    namespace = _load_custom_scan_functions()
    det1 = _detector("det1", 10, 0.02)
    tm1 = _detector("tm1", 50, 0.1)
    namespace["_get_ah_scan_dets"] = lambda primary: [det1, primary]
    namespace["det134"] = SimpleNamespace(name="det134")
    namespace["sclr"] = SimpleNamespace(
        channels=SimpleNamespace(chan13=SimpleNamespace(name="sclr_chan13"))
    )
    namespace["sr_curr"] = SimpleNamespace(name="sr_curr")

    md = {"test": True}
    scan = namespace[scan_name](
        SimpleNamespace(name="motor"),
        -1,
        1,
        3,
        tm1,
        0.1,
        det_ch=[4],
        md=md,
    )

    messages = []
    for _ in range(6):
        messages.append(next(scan) if not messages else scan.send(None))
    scan.close()
    return messages, md


def test_tm1_values_per_read_default():
    namespace = _load_custom_scan_functions()

    assert namespace["_scan_values_per_read"](SimpleNamespace(name="tm1")) == 50
    assert namespace["_scan_values_per_read"](SimpleNamespace(name="det1")) == 1
    assert namespace["_scan_values_per_read"](SimpleNamespace(name="det5")) == 1


def test_dscan_configures_tm1_with_50_and_count_time():
    messages, md = _configuration_messages("dscan")

    assert [(msg.command, msg.args) for msg in messages] == [
        ("set", (1,)),
        ("set", (0.1,)),
        ("wait", ()),
        ("set", (50,)),
        ("set", (0.1,)),
        ("wait", ()),
    ]
    assert md["ah501d_config"]["values_per_read"] == {"det1": 1, "tm1": 50}
    assert md["ah501d_config"]["averaging_time"] == 0.1


def test_ascan_configures_tm1_with_50_and_count_time():
    messages, md = _configuration_messages("ascan")

    assert messages[0].args == (1,)
    assert messages[1].args == (0.1,)
    assert messages[3].args == (50,)
    assert messages[4].args == (0.1,)
    assert md["ah501d_config"]["values_per_read"] == {"det1": 1, "tm1": 50}


def test_tm1_alignment_scans_use_point_one_second():
    module = ast.parse(ALIGNMENT_FILE.read_text(), filename=str(ALIGNMENT_FILE))
    registry = next(
        node.value
        for node in module.body
        if isinstance(node, ast.AnnAssign)
        and isinstance(node.target, ast.Name)
        and node.target.id == "_SCAN_REGISTRY"
    )

    count_times = {}
    for key, value in zip(registry.keys, registry.values):
        keyword = ast.literal_eval(key)
        count_times[keyword] = next(
            ast.literal_eval(item.value)
            for item in value.keywords
            if item.arg == "ct"
        )

    assert {name: count_times[name] for name in ("dcm", "crl", "ugap")} == {
        "dcm": 0.1,
        "crl": 0.1,
        "ugap": 0.1,
    }
    assert {name: count_times[name] for name in ("mcm", "hrm", "ccr", "wcr")} == {
        "mcm": 1,
        "hrm": 1,
        "ccr": 1,
        "wcr": 1,
    }
