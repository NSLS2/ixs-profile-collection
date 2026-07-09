"""
98-alignment.py
===============
Simplified alignment scans with automatic state recording for the IXS beamline.

Called automatically at IPython startup after all device and plan files.

Usage
-----
    RE(alignment_scan("dcm"))
    RE(alignment_scan("crl", note="after top-up"))
    RE(alignment_scan("ugap"))
    RE(alignment_scan("mcm"))
    RE(alignment_scan("hrm"))
    RE(alignment_scan("ccr"))
    RE(alignment_scan("wcr"))

Beamline state (CRL position, HRM position, ring current) is read
automatically at the end of each scan and stored with the record.

Preconditions
-------------
Each scan keyword has an associated list of precondition checks.  If any
check fails the scan is aborted before any hardware moves.  All failing
conditions are reported together so the operator sees the complete picture.

Design note — state readers vs. precondition checks
----------------------------------------------------
_read_crl_state() and _read_hrm_state() return a state string ("in"/"out"/
"unknown") and are used for recording comparability metadata.

_check_crl_in(), _check_hrm_in(), _check_hrm_out() are precondition
checkers that call the state readers internally — no hardware is read
twice and the position logic lives in exactly one place.

Adding a new scan keyword
-------------------------
1. Add precondition check functions if needed (plain callables,
   return (bool, str)).
2. Add one entry to _SCAN_REGISTRY with the scan parameters and checks list.
3. Call RE(alignment_scan("your_keyword")).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import bluesky.plan_stubs as bps

from utils.alignment_store import AlignmentStore, record_alignment


# ---------------------------------------------------------------------------
# Beamline state readers
# Return a state string; called both for recording metadata and by checks.
# ---------------------------------------------------------------------------

def _read_crl_state() -> str:
    """
    CRL state from crl.y motor position.
    IN  : position < 1 mm
    OUT : position >= 1 mm
    """
    try:
        return "in" if crl.y.position < 1.0 else "out"
    except Exception:
        return "unknown"


def _read_hrm_state() -> str:
    """
    HRM state from hrm2.ux and hrm2.dx positions.
    IN  : both axes within ±2 mm of 0
    OUT : either axis outside ±2 mm of 0
    """
    try:
        if abs(hrm2.ux.position) <= 2.0 and abs(hrm2.dx.position) <= 2.0:
            return "in"
        return "out"
    except Exception:
        return "unknown"


def _read_ring_current() -> float | None:
    """
    Storage-ring current from sr_curr (mA).
    Returns None on failure with a RuntimeWarning.
    """
    try:
        return sr_curr.get()
    except Exception as exc:
        warnings.warn(
            f"Could not read sr_curr: {exc}. "
            "ring_current recorded as None; "
            "comparability keying will use mode 'unknown'.",
            RuntimeWarning,
            stacklevel=2,
        )
        return None


# ---------------------------------------------------------------------------
# Precondition check functions
#
# Each returns (ok: bool, message: str).
# message is shown to the operator only when ok is False.
# Checks that test CRL/HRM state delegate to the state readers above so
# that the position logic lives in exactly one place.
# ---------------------------------------------------------------------------

def _check_crl_in() -> tuple[bool, str]:
    """CRL must be IN the beam (crl.y < 1 mm)."""
    state = _read_crl_state()
    if state == "in":
        return True, ""
    if state == "unknown":
        return False, "Could not read crl.y position."
    pos = crl.y.position
    return False, (
        f"CRL is OUT of the beam (crl.y = {pos:.3f} mm). "
        "This scan requires CRL IN (crl.y < 1 mm)."
    )


def _check_hrm_in() -> tuple[bool, str]:
    """HRM must be IN the beam (both hrm2.ux and hrm2.dx within ±2 mm of 0)."""
    state = _read_hrm_state()
    if state == "in":
        return True, ""
    if state == "unknown":
        return False, "Could not read HRM position."
    ux = hrm2.ux.position
    dx = hrm2.dx.position
    return False, (
        f"HRM is OUT of the beam (hrm2.ux = {ux:.3f}, hrm2.dx = {dx:.3f} mm). "
        "This scan requires HRM IN (both axes within ±2 mm of 0)."
    )


def _check_hrm_out() -> tuple[bool, str]:
    """HRM must be OUT of the beam (either hrm2.ux or hrm2.dx outside ±2 mm of 0)."""
    state = _read_hrm_state()
    if state == "out":
        return True, ""
    if state == "unknown":
        return False, "Could not read HRM position."
    ux = hrm2.ux.position
    dx = hrm2.dx.position
    return False, (
        f"HRM is IN the beam (hrm2.ux = {ux:.3f}, hrm2.dx = {dx:.3f} mm). "
        "This scan requires HRM OUT."
    )


def _check_anc_xtal_y_low() -> tuple[bool, str]:
    """anc_xtal.y must be at around 0 (< 1 mm)."""
    try:
        pos = anc_xtal.y.position
        if pos >= 1.0:
            return False, f"anc_xtal.y = {pos:.3f} mm (must be < 1 mm)."
        return True, ""
    except Exception as exc:
        return False, f"Could not read anc_xtal.y: {exc}."


def _check_anc_xtal_y_high() -> tuple[bool, str]:
    """anc_xtal.y must be > 4.5 mm."""
    try:
        pos = anc_xtal.y.position
        if pos <= 4.5:
            return False, f"anc_xtal.y = {pos:.3f} mm (must be > 4.5 mm)."
        return True, ""
    except Exception as exc:
        return False, f"Could not read anc_xtal.y: {exc}."


def _check_hrm_d5_at_2() -> tuple[bool, str]:
    """hrm2.d5 must be at position 2 (±0.1 mm)."""
    _TARGET = 2.0
    _TOL    = 0.1
    try:
        pos = hrm2.d5.position
        if abs(pos - _TARGET) > _TOL:
            return False, (
                f"hrm2.d5 = {pos:.3f} mm (must be at {_TARGET} ± {_TOL} mm)."
            )
        return True, ""
    except Exception as exc:
        return False, f"Could not read hrm2.d5: {exc}."


def _check_anpd_at(target: float, tol: float = 0.1):
    """Return a check function that verifies anpd is at *target* ± *tol* mm."""
    def _check() -> tuple[bool, str]:
        try:
            pos = anpd.position
            if abs(pos - target) > tol:
                return False, (
                    f"anpd = {pos:.3f} mm (must be at {target:+.0f} ± {tol} mm)."
                )
            return True, ""
        except Exception as exc:
            return False, f"Could not read anpd: {exc}."
    _check.__name__ = f"_check_anpd_at_{target:+.0f}"
    return _check


def _check_whl_at(target: float, tol: float = 0.1):
    """Return a check function that verifies whl is at *target* ± *tol*."""
    def _check() -> tuple[bool, str]:
        try:
            pos = whl.position
            if abs(pos - target) > tol:
                return False, (
                    f"whl = {pos:.3f} (must be at {target:+.0f} ± {tol})."
                )
            return True, ""
        except Exception as exc:
            return False, f"Could not read whl: {exc}."
    _check.__name__ = f"_check_whl_at_{target:+.0f}"
    return _check


def _check_whl_for_mcm() -> tuple[bool, str]:
    """
    whl position for MCM scan depends on HRM state:
      HRM OUT -> whl must be at 2 (±0.1)
      HRM IN  -> whl must be at 0 (±0.1)
    Aborts if HRM state cannot be read.
    """
    _TOL = 0.1
    hrm_state = _read_hrm_state()
    if hrm_state == "unknown":
        return False, (
            "Could not read HRM state — cannot determine required "
            "whl position for MCM scan."
        )
    target = 2 if hrm_state == "out" else 0
    try:
        pos = whl.position
        if abs(pos - target) > _TOL:
            return False, (
                f"whl = {pos:.3f} (HRM is {hrm_state}, "
                f"must be at {target:+.0f} ± {_TOL})."
            )
        return True, ""
    except Exception as exc:
        return False, f"Could not read whl: {exc}."


def _check_det2_em_range_0() -> tuple[bool, str]:
    """det2.em_range must be set to 0."""
    try:
        val = det2.em_range.get()
        if int(val) != 0:
            return False, (
                f"det2.em_range = {val} (must be 0). "
                "Run: det2.em_range.set(0)"
            )
        return True, ""
    except Exception as exc:
        return False, f"Could not read det2.em_range: {exc}."


def _check_bpm1_y_at_zero() -> tuple[bool, str]:
    """bpm1.y must be near 0 (within ±0.1 mm)."""
    _TARGET = 0.0
    _TOL    = 0.1
    try:
        pos = bpm1.y.position
        if abs(pos - _TARGET) > _TOL:
            return False, (
                f"bpm1.y = {pos:.3f} mm (must be at {_TARGET} ± {_TOL} mm)."
            )
        return True, ""
    except Exception as exc:
        return False, f"Could not read bpm1.y: {exc}."


def _make_slit_check(targets: dict[str, float], tol: float = 0.1):
    """
    Return a check function that verifies analyzer_slits are at *targets*
    within *tol* mm.

    Parameters
    ----------
    targets : dict
        e.g. {"top": 2.0, "bottom": -2.0, "outboard": 2.0, "inboard": -2.0}
    tol : float
        Tolerance in mm.
    """
    _BLADES = {
        "top":      lambda: analyzer_slits.top.position,
        "bottom":   lambda: analyzer_slits.bottom.position,
        "outboard": lambda: analyzer_slits.outboard.position,
        "inboard":  lambda: analyzer_slits.inboard.position,
    }
    def _check() -> tuple[bool, str]:
        try:
            bad = []
            for blade, target in targets.items():
                pos = _BLADES[blade]()
                if abs(pos - target) > tol:
                    bad.append(f"{blade} = {pos:.3f} (expected {target:+.3g})")
            if bad:
                expected = ", ".join(
                    f"{b}={v:+.3g}" for b, v in targets.items()
                )
                return False, (
                    "Analyzer slits not at required position: "
                    + ", ".join(bad)
                    + f". Expected {expected} (±{tol} mm)."
                )
            return True, ""
        except Exception as exc:
            return False, f"Could not read analyzer_slits: {exc}."
    _check.__name__ = "_check_analyzer_slits"
    return _check


# Pre-built slit check instances used in the registry
_check_slits_open = _make_slit_check(
    {"top": 2.0, "bottom": -2.0, "outboard": 2.0, "inboard": -2.0},
    tol=0.1,
)
_check_slits_narrow = _make_slit_check(
    {"top": 0.1, "bottom": -0.1, "outboard": 1.0, "inboard": -1.0},
    tol=0.01,
)


# ---------------------------------------------------------------------------
# Alignment store — one shared instance per session
# Log files written to the IXS legacy data directory.
# ---------------------------------------------------------------------------

_LOG_DIR = Path("/nsls2/data/ixs/legacy")

store = AlignmentStore(
    jsonl_path  = _LOG_DIR / "alignment_log.jsonl",
    sqlite_path = _LOG_DIR / "alignment_log.db",
)


# ---------------------------------------------------------------------------
# Scan registry
#
# Maps string keyword -> scan configuration dict.
#
# Keys in each configuration dict:
#   motor      : lambda returning the ophyd positioner to scan
#   start      : relative scan start (passed to dscan)
#   stop       : relative scan stop
#   steps      : number of points
#   detector   : lambda returning the ophyd detector object
#   ct         : count time (seconds)
#   det_ch     : detector channel list for dscan
#   instrument : label stored in AlignmentRecord
#   step       : step label stored in AlignmentRecord
#   det_label  : human-readable detector string stored in AlignmentRecord
#   mot_label  : human-readable motor string stored in AlignmentRecord
#   checks     : list of precondition callables (each returns (bool, str))
#
# motor and detector are lambdas so that device objects are resolved at
# call time (when all startup files have run), not at import time.
# ---------------------------------------------------------------------------

_SCAN_REGISTRY: dict[str, dict] = {
    "dcm": dict(
        motor      = lambda: dcm.p1,
        start      = -80,  stop = 80,  steps = 40,
        detector   = lambda: tm1,
        ct         = 1,    det_ch = [4],
        instrument = "dcm",
        step       = "p1_scan",
        det_label  = "tm1",
        mot_label  = "dcm.p1",
        checks     = [_check_hrm_out, _check_bpm1_y_at_zero],
    ),
    "crl": dict(
        motor      = lambda: crl.y,
        start      = -0.2, stop = 0.2, steps = 20,
        detector   = lambda: tm1,
        ct         = 1,    det_ch = [4],
        instrument = "crl",
        step       = "y_scan",
        det_label  = "tm1",
        mot_label  = "crl.y",
        checks     = [_check_crl_in, _check_bpm1_y_at_zero],
    ),
    "ugap": dict(
        motor      = lambda: ivu22,
        start      = -20,  stop = 20,  steps = 21,
        detector   = lambda: tm1,
        ct         = 1,    det_ch = [4],
        instrument = "ivu22",
        step       = "gap_scan",
        det_label  = "tm1",
        mot_label  = "ivu22",
        checks     = [_check_bpm1_y_at_zero],
    ),
    "mcm": dict(
        motor      = lambda: mcm.y,
        start      = -0.2, stop = 0.2, steps = 40,
        detector   = lambda: det2,
        ct         = 1,    det_ch = [0],
        instrument = "mcm",
        step       = "y_scan",
        det_label  = "det2",
        mot_label  = "mcm.y",
        checks     = [_check_anc_xtal_y_low,
                      _check_slits_open,
                      _check_whl_for_mcm,
                      _check_bpm1_y_at_zero],
    ),
    "hrm": dict(
        motor      = lambda: hrm2.dif,
        start      = -70,  stop = 70,  steps = 51,
        detector   = lambda: det5,
        ct         = 1,    det_ch = [0],
        instrument = "hrm",
        step       = "dif_scan",
        det_label  = "det5",
        mot_label  = "hrm2.dif",
        checks     = [_check_hrm_d5_at_2,
                      _check_hrm_in,
                      _check_bpm1_y_at_zero],
    ),
    "ccr": dict(
        motor      = lambda: analyzer.cfth,
        start      = -150, stop = 150, steps = 31,
        detector   = lambda: det2,
        ct         = 1,    det_ch = [0],
        instrument = "ccr",
        step       = "cfth_scan",
        det_label  = "det2",
        mot_label  = "analyzer.cfth",
        checks     = [_check_anc_xtal_y_high,
                      _check_slits_narrow,
                      _check_anpd_at(40),
                      _check_whl_at(0),
                      _check_hrm_out,
                      _check_det2_em_range_0,
                      _check_bpm1_y_at_zero],
    ),
    "wcr": dict(
        motor      = lambda: analyzer.wfth,
        start      = -20,  stop = 20,  steps = 41,
        detector   = lambda: lambda_det,
        ct         = 1,    det_ch = [0],
        instrument = "wcr",
        step       = "wfth_scan",
        det_label  = "lambda_det",
        mot_label  = "analyzer.wfth",
        checks     = [_check_anpd_at(-90),
                      _check_whl_at(7),
                      _check_slits_narrow,
                      _check_hrm_out,
                      _check_bpm1_y_at_zero],
    ),
}


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def alignment_scan(keyword: str, note: str = ""):
    """
    Run a simplified alignment scan and record the result.

    Preconditions for the requested scan are checked before any hardware
    moves.  If any condition is not met the scan is aborted and all
    failing conditions are reported together.

    Parameters
    ----------
    keyword : str
        Instrument keyword.  One of:
        "dcm", "crl", "ugap", "mcm", "hrm", "ccr", "wcr".
    note : str, optional
        Free-text note stored with the AlignmentRecord.

    Examples
    --------
        RE(alignment_scan("dcm"))
        RE(alignment_scan("crl", note="after top-up"))
        RE(alignment_scan("hrm"))
    """
    # --- validate keyword ---------------------------------------------------
    if keyword not in _SCAN_REGISTRY:
        raise ValueError(
            f"alignment_scan: unknown keyword '{keyword}'. "
            f"Known keywords: {sorted(_SCAN_REGISTRY)}."
        )

    cfg = _SCAN_REGISTRY[keyword]

    # --- precondition checks (all evaluated before any yield) ---------------
    failures = []
    for check in cfg["checks"]:
        ok, msg = check()
        if not ok:
            failures.append(msg)

    if failures:
        print(f"\nalignment_scan('{keyword}') ABORTED — conditions not met:")
        for msg in failures:
            print(f"  ! {msg}")
        return

    # --- run the bare dscan (no motor move afterward) -----------------------
    stats = yield from dscan(
        cfg["motor"](),
        cfg["start"], cfg["stop"], cfg["steps"],
        cfg["detector"](), cfg["ct"],
        det_ch=cfg["det_ch"],
    )

    # --- guard: need a valid peak to record ---------------------------------
    if not stats or stats[0].max is None:
        print(
            f"[alignment_scan] '{keyword}': "
            "scan returned no usable peak — result not recorded."
        )
        return

    s    = stats[0]
    fmax = s.max          # (x_at_max, max_intensity)

    # --- read context -------------------------------------------------------
    scan_uid     = db[-1].start["uid"]
    crl_state    = _read_crl_state()
    hrm_state    = _read_hrm_state()
    ring_current = _read_ring_current()

    # --- record -------------------------------------------------------------
    yield from record_alignment(
        store,
        scan_id       = scan_uid,
        instrument    = cfg["instrument"],
        step          = cfg["step"],
        detector      = cfg["det_label"],
        motor         = cfg["mot_label"],
        max_intensity = fmax[1],
        x_at_max      = fmax[0],
        fwhm          = s.fwhm,
        cen           = s.cen,
        com           = s.com,
        crl_state     = crl_state,
        hrm_state     = hrm_state,
        ring_current  = ring_current,
        note          = note,
    )


# ---------------------------------------------------------------------------
# Diagnostic utilities (call directly in IPython, outside a plan)
# ---------------------------------------------------------------------------

def print_best(instrument: str, step: str,
               crl_state:    str | None   = None,
               hrm_state:    str | None   = None,
               ring_current: float | None = None):
    """
    Print the best record for the exact
    (instrument, step, crl_state, hrm_state, ring_current_mode) key.

    If no exact match exists, lists records available under different
    conditions (FOR DIAGNOSTICS ONLY — not for intensity comparison).

    Parameters
    ----------
    instrument : str
        e.g. "dcm", "hrm"
    step : str
        e.g. "p1_scan", "dif_scan"
    crl_state : str or None
        "in" / "out" / "unknown" / None
    hrm_state : str or None
        "in" / "out" / "unknown" / None
    ring_current : float or None
        Raw ring current in mA.  Bucketed internally to "400mA"/"500mA"/"unknown".
    """
    rec = store.get_best(instrument, step,
                         crl_state=crl_state,
                         hrm_state=hrm_state,
                         ring_current=ring_current)
    if rec:
        print(rec.comparison_summary())
        return

    print(
        f"No records for {instrument}/{step} under conditions "
        f"CRL={crl_state}  HRM={hrm_state}  ring_current={ring_current}."
    )
    others = store.get_best_any_condition(instrument, step)
    if others:
        print(
            "  Records exist under different conditions "
            "(FOR DIAGNOSTICS ONLY — do not use for intensity comparison):"
        )
        for r in others:
            print(
                f"    CRL={r.crl_state}  HRM={r.hrm_state}  "
                f"ring={r.ring_current_mode}  "
                f"best={r.max_intensity:.4g}  @ {r.timestamp}"
            )
    else:
        print(f"  No records at all for {instrument}/{step}.")


def print_history(instrument: str | None = None,
                   step:       str | None = None):
    """Tabular printout of all (or filtered) historical records."""
    print(store.report(instrument=instrument, step=step))


# ---------------------------------------------------------------------------
# Interactive log manager
# ---------------------------------------------------------------------------

_PAGE_SIZE = 20   # rows shown per page


def manage_alignment_log(instrument: str | None = None):
    """
    Interactive browser and editor for the alignment log.

    Displays records in pages of 20.  At each page the operator can move
    to the next or previous page, delete a record by its displayed number,
    or quit.

    Parameters
    ----------
    instrument : str or None
        If given, only records for that instrument are shown
        (e.g. ``"dcm"``, ``"crl"``).  Pass nothing to browse all records.

    Examples
    --------
        manage_alignment_log()
        manage_alignment_log("dcm")
    """
    def _load():
        return store.history(instrument=instrument)

    def _fmt_row(idx: int, r) -> str:
        fwhm_s  = f"{r.fwhm:.4g}"  if r.fwhm  is not None else "—"
        ring_s  = r.ring_current_mode
        return (
            f"{idx:>3}  {r.timestamp:<20}  {r.instrument:<10}  "
            f"{r.step:<15}  {r.max_intensity:>13.4g}  "
            f"{fwhm_s:>8}  {r.crl_state or '—':>5}  "
            f"{r.hrm_state or '—':>5}  {ring_s:>6}"
        )

    _HDR = (
        f"{'#':>3}  {'timestamp':<20}  {'instrument':<10}  "
        f"{'step':<15}  {'max_intensity':>13}  "
        f"{'fwhm':>8}  {'CRL':>5}  {'HRM':>5}  {'ring':>6}"
    )
    _SEP = "-" * len(_HDR)

    records = _load()
    if not records:
        label = f"instrument='{instrument}'" if instrument else "all instruments"
        print(f"No alignment records found for {label}.")
        return

    page = 0

    while True:
        total   = len(records)
        n_pages = max(1, (total + _PAGE_SIZE - 1) // _PAGE_SIZE)
        start   = page * _PAGE_SIZE
        end     = min(start + _PAGE_SIZE, total)
        page_records = records[start:end]

        print()
        label = f"  instrument='{instrument}'" if instrument else "  all instruments"
        print(f"Alignment log —{label}   "
              f"page {page + 1}/{n_pages}   "
              f"records {start + 1}–{end} of {total}")
        print(_SEP)
        print(_HDR)
        print(_SEP)
        for i, r in enumerate(page_records, start=start + 1):
            print(_fmt_row(i, r))
        print(_SEP)

        has_next = page < n_pages - 1
        has_prev = page > 0
        nav = []
        if has_next:
            nav.append("[n]ext")
        if has_prev:
            nav.append("[p]rev")
        nav.append("[d]elete <#>")
        nav.append("[q]uit")
        prompt = "  ".join(nav) + "  > "

        try:
            raw = input(prompt).strip().lower()
        except (EOFError, KeyboardInterrupt):
            print()
            break

        if raw in ("q", "quit", "exit"):
            break

        elif raw in ("n", "next"):
            if has_next:
                page += 1
            else:
                print("  Already on the last page.")

        elif raw in ("p", "prev", "previous"):
            if has_prev:
                page -= 1
            else:
                print("  Already on the first page.")

        elif raw.startswith(("d ", "del ", "delete ")):
            # parse the record number
            parts = raw.split(None, 1)
            try:
                num = int(parts[1])
            except (IndexError, ValueError):
                print("  Usage: d <number>   e.g.  d 3")
                continue

            if num < 1 or num > total:
                print(f"  Number out of range (1–{total}).")
                continue

            target = records[num - 1]
            print()
            print(f"  Record #{num}:")
            print(f"    timestamp   : {target.timestamp}")
            print(f"    scan_id     : {target.scan_id}")
            print(f"    instrument  : {target.instrument}")
            print(f"    step        : {target.step}")
            print(f"    max_intensity: {target.max_intensity:.4g}")
            print(f"    fwhm        : {target.fwhm}")
            print(f"    CRL={target.crl_state}  HRM={target.hrm_state}  "
                  f"ring={target.ring_current_mode}")
            if target.note:
                print(f"    note        : {target.note}")
            print()

            try:
                answer = input(
                    f"  Delete record #{num} "
                    f"({target.instrument}/{target.step} @ {target.timestamp})? "
                    "[y/N]  > "
                ).strip().lower()
            except (EOFError, KeyboardInterrupt):
                print()
                break

            if answer == "y":
                removed = store.delete_record(target)
                if removed:
                    print(f"  Deleted record #{num}.")
                else:
                    print("  Record not found in log — nothing deleted.")
                # reload and clamp page to valid range
                records = _load()
                total   = len(records)
                n_pages = max(1, (total + _PAGE_SIZE - 1) // _PAGE_SIZE)
                if page >= n_pages:
                    page = n_pages - 1
                if not records:
                    print("  Log is now empty.")
                    break
            else:
                print("  Deletion cancelled.")

        else:
            if raw:
                print("  Commands:  n  p  d <#>  q")
