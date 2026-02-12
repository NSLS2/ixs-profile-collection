import os
import time
import socket
import getpass
from datetime import datetime
from bluesky.callbacks.core import CallbackBase

SPEC_TIME_FORMAT = "%a %b %d %H:%M:%S %Y"

def _format_spec_command(start):
    pn = start.get("plan_name", "unknown")
    dets = list(start.get("detectors", []))  # these are strings like ['det1']
    motors = list(start.get("motors", []))   # strings like ['hrmE'] or ('mcm_x',)

    pa = start.get("plan_args", {}) or {}

    # COUNT / CT
    if pn in ("count", "trigger_and_read"):
        # count() uses plan_args: num, delay (and detectors but those are verbose objects)
        num = pa.get("num", start.get("num_points", 1))
        delay = pa.get("delay", 0)
        if dets:
            return f"count({dets}, num={num}, delay={delay})"
        return f"count(num={num}, delay={delay})"

    # SCANS (rel_scan / scan / ascan / dscan etc.)
    # We aim for something SPEC-ish: scan motor start stop num time
    if pn in ("rel_scan", "scan", "rel_grid_scan", "grid_scan"):
        # For rel_scan/scan, bluesky stores arguments in plan_args['args'] and plan_args['num']
        args = pa.get("args", [])
        num = pa.get("num", start.get("num_points", None))

        # Common case: args = [motor, start, stop]  (rel_scan) OR [motor, start, stop] (scan)
        # Sometimes args is stringified; handle both.
        if len(motors) == 1 and len(args) >= 3:
            m = motors[0]
            start_val = args[-2]
            stop_val = args[-1]
            # count_time if you track it
            ct = start.get("count_time", pa.get("per_step", None))
            if num is None:
                num = start.get("num_points", "")
            # If ct is None/ugly, omit it
            if isinstance(ct, (int, float)) and ct >= 0:
                return f"{pn} {m} {start_val} {stop_val} {num} {ct}"
            return f"{pn} {m} {start_val} {stop_val} {num}"

    # Fallback: short summary, no giant repr
    pieces = [pn]
    if motors:
        pieces.append(f"motors={motors}")
    if dets:
        pieces.append(f"detectors={dets}")
    # include only scalar-ish args
    scalar_keys = ("num", "delay", "step", "per_step")
    scalars = {k: pa[k] for k in scalar_keys if k in pa and isinstance(pa[k], (int, float, str, bool, type(None)))}
    if scalars:
        pieces.append(f"args={scalars}")
    return " ".join(pieces)


def to_spec_time(ts):
    return datetime.fromtimestamp(ts).strftime(SPEC_TIME_FORMAT)

def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def _fmt(x):
    # SPEC-ish numeric formatting
    try:
        return f"{float(x):.12g}"
    except Exception:
        return str(x)
    
import numbers

def _read_scalar(obj):
    # raw numbers (float, int, numpy scalar, etc.)
    if isinstance(obj, numbers.Number):
        return float(obj)

    # callable
    if callable(obj):
        return obj()

    # ophyd-style get()
    get_fn = getattr(obj, "get", None)
    if get_fn is not None:
        try:
            val = get_fn()
            if hasattr(val, "readback"):
                return float(val.readback)
            return float(val)
        except Exception:
            pass

    # ophyd-style read()
    read_fn = getattr(obj, "read", None)
    if read_fn is not None:
        try:
            d = read_fn()
            if isinstance(d, dict) and d:
                v = next(iter(d.values()))["value"]
                return float(v)
        except Exception:
            pass

    return None

def _write_G0_line(fh, g0_items):
    """
    Write one SPEC #G0 line from a fixed ordered list of 11 scalars.
    g0_items: list of ophyd Signals/Positioners or callables returning scalars.
    """
    vals = []
    for it in g0_items:
        v = _read_scalar(it)  # use the same helper you already have (callable/get/read)
        vals.append(v)
    fh.write("#G0 " + " ".join(_fmt(v) for v in vals) + "\n")


def _safe_get_position(obj):
    """
    Try common ophyd conventions to get a scalar position for header #P lines.
    """
    for attr in ("user_readback", "readback", "position"):
        sig = getattr(obj, attr, None)
        if sig is not None:
            try:
                return sig.get()
            except Exception:
                pass
    try:
        return obj.get()
    except Exception:
        return None

class CustomSpecWriter(CallbackBase):
    """
    Minimal, clean SPEC-like writer:
    - File header: #F #E #D #C and curated #O/#o blocks.
    - Per scan: #S #D #C ... #MD ... #P blocks, then #N/#L and data rows.

    Designed for single primary stream and step scans; can be extended.
    """
    def __init__(
        self,
        filepath,
        motor_groups,
        *,
        motors_per_line=8,
        include_md_keys=None,
        x_field_resolver=None,
        data_field_order=None,
        g0_items=None,
        flush=True,
    ):
        """
        filepath: full output path (e.g. /.../spec1.dat)
        motor_groups: list of tuples (group_name, motors)
            motors: list of ophyd objects (Positioners or Signals)
            names used in #O/#o are motor.name unless overridden with ._spec_name attr
        motors_per_line: how many motor names/positions per #O/#P line
        include_md_keys: iterable of RunStart keys to write as #MD key = value
        x_field_resolver: function(start_doc, primary_descriptor_doc, event_doc) -> field_name for x axis
        data_field_order: function(start_doc, primary_descriptor_doc) -> list of data keys (excluding x) in desired order
        """
        self.filepath = filepath
        self.motor_groups = list(motor_groups)
        self.motors_per_line = int(motors_per_line)
        self.include_md_keys = set(include_md_keys or ())
        self.x_field_resolver = x_field_resolver
        self.data_field_order = data_field_order
        self.flush = bool(flush)
        self.g0_items = list(g0_items) if g0_items is not None else []

        self._fh = None
        self._file_header_written = False

        self._start = None
        self._primary_desc = None
        self._primary_desc_uid = None

        # Flatten motor list and names in SPEC order
        self._motors = []
        self._motor_names = []
        for _, motors in self.motor_groups:
            for m in motors:
                self._motors.append(m)
                self._motor_names.append(getattr(m, "_spec_name", m.name))

    def _open(self):
        os.makedirs(os.path.dirname(self.filepath), exist_ok=True)
        self._fh = open(self.filepath, "a", buffering=1)

    def _write_file_header_once(self, start_doc):
        if self._file_header_written:
            return
        if self._fh is None:
            self._open()
        # Only write if file is empty
        if self._fh.tell() != 0:
            self._file_header_written = True
            return

        now = start_doc.get("time", time.time())
        user = getpass.getuser()
        host = socket.gethostname()

        self._fh.write(f"#F {os.path.basename(self.filepath)}\n")
        self._fh.write(f"#E {int(now)}\n")
        self._fh.write(f"#D {to_spec_time(now)}\n")
        self._fh.write(f"#C Bluesky  user = {user}  host = {host}\n")

        # #O/#o motor name blocks, split into multiple lines
        for i, names in enumerate(_chunks(self._motor_names, self.motors_per_line)):
            self._fh.write(f"#O{i} " + "  ".join(names) + "\n")
        for i, names in enumerate(_chunks(self._motor_names, self.motors_per_line)):
            self._fh.write(f"#o{i} " + "  ".join(names) + "\n")

        self._fh.write("\n")
        self._file_header_written = True
        if self.flush:
            self._fh.flush()

    def start(self, doc):
        self._start = doc
        self._primary_desc = None
        self._primary_desc_uid = None
        self._write_file_header_once(doc)

        # Write scan header start immediately (APS-style: #S then #D then #C/#MD/#P/#N/#L later)
        scan_id = doc.get("scan_id", "?")
        # plan_name = doc.get("plan_name", "unknown_plan")
        # plan_args = doc.get("plan_args", {})
        cmd = _format_spec_command(doc)
        self._fh.write(f"\n#S {scan_id} {cmd}\n")
        self._fh.write(f"#D {to_spec_time(doc.get('time', time.time()))}\n")
        ct = doc.get("count_time", None)
        if ct is not None:
            self._fh.write(f"#T {_fmt(ct)}  (Seconds)\n")
        if self.g0_items:
            _write_G0_line(self._fh, self.g0_items)
        self._fh.write(f"#C {to_spec_time(doc.get('time', time.time()))}.  uid = {doc.get('uid', '')}\n")

        # Selected #MD lines (curated, not everything)
        for k in sorted(self.include_md_keys):
            if k in doc:
                self._fh.write(f"#MD {k} = {doc[k]}\n")

        # Write #P blocks using live motor positions aligned to #O order
        positions = [_safe_get_position(m) for m in self._motors]
        for i, vals in enumerate(_chunks(positions, self.motors_per_line)):
            self._fh.write(f"#P{i} " + " ".join(_fmt(v) for v in vals) + "\n")

        if self.flush:
            self._fh.flush()

    def descriptor(self, doc):
        # pick primary stream only
        if doc.get("name") != "primary":
            return
        self._primary_desc = doc
        self._primary_desc_uid = doc["uid"]

        # Decide columns (#L) now that we know data_keys
        # x_field: choose scanned motor field for first column
        if self.x_field_resolver is None:
            # default: use hinted dimension field if present, else first motor in start['motors'], else seq_num
            def _default_x(start, desc, _event):
                hints = start.get("hints", {})
                dims = hints.get("dimensions", [])
                if dims:
                    fields = dims[0][0]
                    for f in fields:
                        if f in desc.get("data_keys", {}):
                            return f
                motors = start.get("motors", ())
                if motors and motors[0] in desc.get("data_keys", {}):
                    return motors[0]
                return "seq_num"
            self.x_field_resolver = _default_x

        if self.data_field_order is None:
            # default: all scalar fields except x, preserve descriptor insertion order
            def _default_order(start, desc):
                x = self.x_field_resolver(start, desc, {})
                out = []
                for k, v in desc["data_keys"].items():
                    if k == x:
                        continue
                    if v.get("shape"):
                        continue
                    out.append(k)
                return out
            self.data_field_order = _default_order

        # Build #N/#L
        x = self.x_field_resolver(self._start, doc, {})
        y_fields = self.data_field_order(self._start, doc)
        # Provide APS-like Epoch_float and Epoch columns
        col_labels = [x, "Epoch_float", "Epoch"] + y_fields
        self._fh.write(f"#N {len(col_labels)}\n")
        self._fh.write("#L " + "  ".join(col_labels) + "\n")
        if self.flush:
            self._fh.flush()

    def event(self, doc):
        if doc.get("descriptor") != self._primary_desc_uid:
            return

        # Determine x and y fields
        x = self.x_field_resolver(self._start, self._primary_desc, doc)
        data = doc.get("data", {})
        ts = doc.get("time", time.time())

        x_val = doc.get("seq_num") if x == "seq_num" else data.get(x, doc.get("seq_num"))
        epoch_float = float(ts)
        epoch_int = int(ts)

        y_fields = self.data_field_order(self._start, self._primary_desc)
        y_vals = [data.get(k, "") for k in y_fields]

        self._fh.write(
            f"{_fmt(x_val)}  {_fmt(epoch_float)}  {epoch_int}  " + "  ".join(_fmt(v) for v in y_vals) + "\n"
        )
        if self.flush:
            self._fh.flush()

    def stop(self, doc):
        # APS-style tail comments
        t = to_spec_time(time.time())
        exit_status = doc.get("exit_status", "")
        reason = doc.get("reason", "")
        self._fh.write(f"#C {t}.  exit_status = {exit_status}\n")
        if exit_status != "success":
            self._fh.write(f"#C {t}.  reason = {reason}\n")
        if self.flush:
            self._fh.flush()
        # keep file open across runs if you want a single spec file; otherwise close here
        # self._fh.close()
