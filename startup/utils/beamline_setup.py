"""Execution framework for the automatic beamline alignment sequence.

Step 1 of the Beamline_Setup work: ordering, interaction modes, failure
handling and reporting only. The module knows nothing about devices, EPICS or
the real alignment plans. The caller injects an ordered step table (see
``Step``); there is no default table, so simulated actions can never be
substituted silently in a production entry point.

``start``/``stop`` only select a range of table steps. No physical-state
checks are done, so ``start=N`` is NOT a validated way to resume a run.
"""

from __future__ import annotations

import datetime
import json
import os
import time
import uuid
from dataclasses import dataclass
from typing import Any, Callable, Optional

from bluesky.utils import RequestAbort, RequestStop

try:  # tabulate is a profile dependency; fall back to plain text without it
    from tabulate import tabulate
except ImportError:  # pragma: no cover
    tabulate = None

KIND_AUTO, KIND_MANUAL, KIND_OPTIONAL = "auto", "manual", "optional"
KINDS = (KIND_AUTO, KIND_MANUAL, KIND_OPTIONAL)

COMPLETED, FAILED, CANCELLED = "completed", "failed", "cancelled"
INTERRUPTED, SKIPPED, NOT_RUN = "interrupted", "skipped", "not_run"

SKIP_MANUAL = "not performed (manual)"
SKIP_OPTIONAL = "skipped (optional)"
OUTSIDE_RANGE = "outside_selected_range"

OUTCOME_COMPLETED = "completed"
OUTCOME_MANUAL_REMAINS = "completed_with_manual_work_incomplete"
OUTCOME_FAILED = "failed"
OUTCOME_CANCELLED = "cancelled"
OUTCOME_INTERRUPTED = "interrupted"
OUTCOME_ERROR = "error"

RANGE_NOTICE = ("start/stop select a range of steps only. No physical-state "
                "validation is done; this is not a validated resume.")
DEFAULT_REPORT_DIR = "beamline_setup_reports"


class AlignmentFailure(RuntimeError):
    """An alignment step failed. The sequence stops."""


class OperatorAbort(RuntimeError):
    """The operator cancelled the sequence."""


@dataclass(frozen=True)
class Step:
    """One entry of the step table.

    id          '1', '2', ... for base steps; '<parent>a' for an optional scan.
    kind        'auto', 'manual' or 'optional'.
    factory     zero-argument callable returning a plan (auto/optional).
    instruction text shown to the operator (manual).
    """

    id: str
    name: str
    kind: str
    factory: Optional[Callable[[], Any]] = None
    instruction: str = ""


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

def _base_number(step_id):
    return int(step_id.rstrip("a"))


def validate_step_table(step_table):
    """Raise ValueError if the table is malformed. Return the base-step count.

    Rules: base steps are numbered 1..N in table order; an optional scan has
    id '<N>a' and directly follows auto base step N; auto/optional steps need
    a callable factory; manual steps need a non-empty instruction.
    """
    if not isinstance(step_table, (list, tuple)) or not step_table:
        raise ValueError("step_table must be a non-empty list of Step")
    n = 0
    prev = None
    for index, step in enumerate(step_table):
        if not isinstance(step, Step):
            raise ValueError(f"step_table[{index}] is not a Step")
        if step.kind not in KINDS:
            raise ValueError(f"step {step.id!r}: invalid kind {step.kind!r}")
        if not isinstance(step.name, str) or not step.name.strip():
            raise ValueError(f"step_table[{index}]: empty name")
        if step.kind == KIND_OPTIONAL:
            if prev is None or prev.kind != KIND_AUTO:
                raise ValueError(
                    f"optional step {step.id!r} must directly follow an auto step")
            expected = f"{n}a"
        else:
            n += 1
            expected = str(n)
        if step.id != expected:
            raise ValueError(
                f"step_table[{index}]: id {step.id!r}, expected {expected!r}")
        if step.kind == KIND_MANUAL:
            if not isinstance(step.instruction, str) or not step.instruction.strip():
                raise ValueError(f"manual step {step.id!r} needs an instruction")
        elif not callable(step.factory):
            raise ValueError(f"step {step.id!r} needs a callable factory")
        prev = step
    return n


def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _validate_arguments(step_table, start, stop, interactive, optional_scans,
                        confirm, report_dir, clock, printer,
                        failure_types=(), abort_types=()):
    count = validate_step_table(step_table)
    if stop is None:
        stop = count
    if not _is_int(start) or not _is_int(stop):
        raise ValueError("start and stop must be integers")
    if not 1 <= start <= stop <= count:
        raise ValueError(
            f"need 1 <= start <= stop <= {count}; got start={start}, stop={stop}")
    if not isinstance(interactive, bool) or not isinstance(optional_scans, bool):
        raise ValueError("interactive and optional_scans must be bool")
    for label, fn in (("confirm", confirm), ("clock", clock), ("printer", printer)):
        if not callable(fn):
            raise ValueError(f"{label} must be callable")
    for label, types in (("failure_types", failure_types),
                         ("abort_types", abort_types)):
        if not isinstance(types, tuple) or not all(
                isinstance(t, type) and issubclass(t, Exception) for t in types):
            raise ValueError(f"{label} must be a tuple of Exception classes")
    report_dir = DEFAULT_REPORT_DIR if report_dir is None else os.fspath(report_dir)
    try:
        os.makedirs(report_dir, exist_ok=True)
    except OSError as exc:
        raise ValueError(f"report_dir {report_dir!r} is not usable: {exc}") from exc
    return stop, report_dir


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def _iso(epoch):
    if epoch is None:
        return None
    return datetime.datetime.fromtimestamp(
        epoch, tz=datetime.timezone.utc).isoformat(timespec="milliseconds")


def _json_safe(value):
    try:
        json.dumps(value)
        return value
    except (TypeError, ValueError):
        return repr(value)


def _new_record(step):
    return {"id": step.id, "name": step.name, "kind": step.kind,
            "status": NOT_RUN, "started_at": None, "ended_at": None,
            "elapsed_s": None, "result": None, "reason": None}


def format_summary(report):
    """Readable text built only from the report dict."""
    opts = report["options"]
    lines = [
        "=" * 72,
        f"Beamline_Setup report: {report['outcome'].upper()}",
        f"  range: steps {opts['start']}..{opts['stop']} of {opts['table_steps']}"
        f"  interactive={opts['interactive']}"
        f"  optional_scans={opts['optional_scans']}",
        f"  {report['range_notice']}",
    ]
    stopped = report["stopped_at"]
    if stopped:
        lines.append(f"  stopped at step {stopped['id']} ({stopped['name']}): "
                     f"{stopped['status']} - {stopped['reason']}")
    manual = report["manual_work"]
    if manual["incomplete"]:
        lines.append("  MANUAL WORK REMAINS INCOMPLETE: steps "
                     + ", ".join(manual["steps"]))
    rows = [[s["id"], s["name"], s["status"],
             "" if s["elapsed_s"] is None else f"{s['elapsed_s']:.2f}",
             s["reason"] or ""] for s in report["steps"]]
    headers = ["id", "step", "status", "elapsed s", "reason"]
    if tabulate is not None:
        lines.append(tabulate(rows, headers=headers))
    else:  # pragma: no cover
        lines.extend(" | ".join(map(str, r)) for r in [headers] + rows)
    for st in report["steps"]:
        res = st["result"]
        if isinstance(res, dict) and "image_path" in res:
            lines.append(f"  step {st['id']} image: "
                         f"{res.get('image_path_relative') or res['image_path']}"
                         + (f" ({res['image_note']})" if res.get("image_note") else ""))
            if res.get("centroid_deviation") is not None:
                lines.append(f"    centroid deviation from targets: "
                             f"{res['centroid_deviation']}")
    if report["report_path"]:
        lines.append(f"  report file: {report['report_path']}")
    if report["report_error"]:
        lines.append(f"  REPORT NOT SAVED: {report['report_error']}")
    lines.append("=" * 72)
    return "\n".join(lines)


def _write_report(report, report_dir, clock):
    """Write the JSON report to a new unique file; return its path."""
    stamp = datetime.datetime.fromtimestamp(
        clock(), tz=datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = os.path.join(
        report_dir, f"beamline_setup_{stamp}_{uuid.uuid4().hex[:8]}.json")
    report["report_path"] = path
    with open(path, "x", encoding="utf-8") as fh:  # "x": never overwrite
        json.dump(report, fh, indent=2)
    return path


def _build_report(run, finished):
    records = run["records"]
    manual = [r["id"] for r in records
              if r["kind"] == KIND_MANUAL and r["reason"] == SKIP_MANUAL]
    stopped = run["stopped"]
    if stopped is None:
        outcome = OUTCOME_MANUAL_REMAINS if manual else OUTCOME_COMPLETED
    else:
        outcome = stopped["outcome"]
    return {
        "schema_version": 1,
        "outcome": outcome,
        "options": run["options"],
        "range_notice": RANGE_NOTICE,
        "started_at": _iso(run["t0"]),
        "ended_at": _iso(finished),
        "elapsed_s": finished - run["t0"],
        "stopped_at": None if stopped is None else
        {k: stopped[k] for k in ("id", "name", "status", "reason")},
        "manual_work": {"incomplete": bool(manual), "steps": manual},
        "steps": records,
        "report_path": None,
        "report_error": None,
    }


def _finalize(run, report_dir, clock, printer):
    """Build, print and save the report. Never raises."""
    try:
        report = _build_report(run, clock())
    except Exception as exc:  # pragma: no cover - defensive
        return {"outcome": OUTCOME_ERROR, "options": run["options"],
                "report_error": f"report build failed: {exc!r}"}
    try:
        _write_report(report, report_dir, clock)
    except Exception as exc:
        report["report_path"] = None
        report["report_error"] = f"{type(exc).__name__}: {exc}"
    try:
        printer(format_summary(report))
    except Exception as exc:
        if report["report_error"] is None:
            report["report_error"] = f"summary not printed: {exc!r}"
    return report


# --------------------------------------------------------------------------
# Execution
# --------------------------------------------------------------------------

def _is_yes(answer):
    if isinstance(answer, bool):
        return answer
    return str(answer).strip().lower() in ("", "y", "yes")


def beamline_setup_plan(step_table, start=1, stop=None, interactive=True,
                        optional_scans=False, confirm=input,
                        report_dir=None, clock=time.time, printer=print,
                        failure_types=(), abort_types=()):
    """Run the injected step table (framework only; no real procedures).

    Arguments are validated, and the table checked, before the plan is
    returned, so nothing runs on bad input. Returns a generator plan whose
    return value is the report dict.

    start, stop    inclusive range of base-step numbers (stop=None: last step).
                   An optional scan follows its parent step's range.
    interactive    confirm before each step. Silent mode never calls confirm,
                   and skips manual steps ("not performed (manual)").
    optional_scans run optional scans when True, else record them skipped.
    confirm        callable(prompt) -> 'yes'/''/True to accept, else decline.
    report_dir     directory for the JSON report (default ./beamline_setup_reports).
    failure_types  extra exception classes handled like AlignmentFailure
                   (e.g. the startup file's own class).
    abort_types    extra exception classes handled like OperatorAbort.
    """
    stop, report_dir = _validate_arguments(
        step_table, start, stop, interactive, optional_scans, confirm,
        report_dir, clock, printer, failure_types, abort_types)
    return _execute(tuple(step_table), start, stop, interactive,
                    optional_scans, confirm, report_dir, clock, printer,
                    (AlignmentFailure,) + failure_types,
                    (OperatorAbort,) + abort_types)


def _execute(table, start, stop, interactive, optional_scans, confirm,
             report_dir, clock, printer, failure_types, abort_types):
    run = {
        "t0": clock(),
        "options": {"start": start, "stop": stop, "interactive": interactive,
                    "optional_scans": optional_scans,
                    "table_steps": _base_number(table[-1].id)},
        "records": [_new_record(s) for s in table],
        "stopped": None,
    }

    def stop_run(rec, status, outcome, reason):
        rec["status"], rec["reason"] = status, reason
        run["stopped"] = {"id": rec["id"], "name": rec["name"], "status": status,
                          "reason": reason, "outcome": outcome}

    for step, rec in zip(table, run["records"], strict=True):
        if not start <= _base_number(step.id) <= stop:
            rec["reason"] = OUTSIDE_RANGE
            continue
        if run["stopped"] is not None:
            s = run["stopped"]
            rec["reason"] = f"prevented_by_{s['status']}_step_{s['id']}"
            continue
        if step.kind == KIND_OPTIONAL and not optional_scans:
            rec["status"], rec["reason"] = SKIPPED, SKIP_OPTIONAL
            continue
        if step.kind == KIND_MANUAL and not interactive:
            rec["status"], rec["reason"] = SKIPPED, SKIP_MANUAL
            continue

        started = None
        try:
            if interactive:
                if step.kind == KIND_MANUAL:
                    prompt = (f"MANUAL step {step.id}: {step.name}\n"
                              f"{step.instruction}\nCompleted? (yes/no): ")
                else:
                    prompt = f"Step {step.id}: {step.name}. Proceed? (yes/no): "
                if not _is_yes(confirm(prompt)):
                    if step.kind == KIND_OPTIONAL:
                        rec["status"] = SKIPPED
                        rec["reason"] = "declined by operator (optional)"
                    else:
                        stop_run(rec, CANCELLED, OUTCOME_CANCELLED,
                                 "declined by operator")
                    continue
            started = clock()
            rec["started_at"] = _iso(started)
            if step.kind == KIND_MANUAL:
                result = None  # operator confirmed completion
            else:
                result = yield from step.factory()
            ended = clock()
            rec.update(status=COMPLETED, ended_at=_iso(ended),
                       elapsed_s=ended - started, result=_json_safe(result))
        except Exception as exc:  # RequestAbort/RequestStop are Exceptions
            ended = clock()
            if started is None:
                rec["started_at"] = _iso(ended)
                started = ended
            rec["ended_at"], rec["elapsed_s"] = _iso(ended), ended - started
            if isinstance(exc, (RequestAbort, RequestStop)):
                stop_run(rec, INTERRUPTED, OUTCOME_INTERRUPTED,
                         f"{type(exc).__name__}")
            elif isinstance(exc, abort_types):
                stop_run(rec, CANCELLED, OUTCOME_CANCELLED,
                         f"operator abort: {exc}")
            elif isinstance(exc, failure_types):
                stop_run(rec, FAILED, OUTCOME_FAILED, f"alignment failure: {exc}")
            else:
                stop_run(rec, FAILED, OUTCOME_ERROR,
                         f"unexpected {type(exc).__name__}: {exc}")
            if isinstance(exc, abort_types + failure_types):
                continue
            report = _finalize(run, report_dir, clock, printer)
            if isinstance(report, dict) and report.get("report_error"):
                try:
                    exc.add_note(f"Beamline_Setup report problem: "
                                 f"{report['report_error']}")
                except Exception:  # pragma: no cover
                    pass
            raise
    return _finalize(run, report_dir, clock, printer)


# --------------------------------------------------------------------------
# Table display and usage
# --------------------------------------------------------------------------

def format_step_table(step_table, start=1, stop=None, optional_scans=False):
    """Text listing of the table; marks what would run for the given options."""
    count = validate_step_table(step_table)
    stop = count if stop is None else stop
    if not (_is_int(start) and _is_int(stop) and 1 <= start <= stop <= count):
        raise ValueError(f"need 1 <= start <= stop <= {count}")
    rows = []
    for step in step_table:
        if not start <= _base_number(step.id) <= stop:
            mark = "outside range"
        elif step.kind == KIND_OPTIONAL and not optional_scans:
            mark = SKIP_OPTIONAL
        else:
            mark = "run"
        note = step.instruction if step.kind == KIND_MANUAL else (
            f"optional scan of step {step.id[:-1]}" if step.kind == KIND_OPTIONAL else "")
        rows.append([step.id, step.name, step.kind, mark, note])
    headers = ["id", "step", "kind", "with these options", "note"]
    if tabulate is not None:
        return tabulate(rows, headers=headers)
    return "\n".join(" | ".join(map(str, r)) for r in [headers] + rows)  # pragma: no cover


def print_step_table(step_table, start=1, stop=None, optional_scans=False,
                     printer=print):
    """Print the step table. Runs nothing. Returns the text."""
    text = format_step_table(step_table, start, stop, optional_scans)
    printer(text)
    return text


USAGE = """\
Beamline_Setup / beamline_setup_plan usage
------------------------------------------
  RE(Beamline_Setup())                       interactive, confirm every step
  RE(Beamline_Setup(optional_scans=True))    also run the optional scans
  RE(Beamline_Setup(start=4, stop=8))        run a range of steps
  print_step_table(...) / Beamline_Setup_table()   show the table, run nothing

Arguments
  start, stop     inclusive range of base-step numbers (stop=None: last step).
                  Range selection ONLY: no physical-state checks are made and
                  this is not a validated resume.
  interactive     True: confirm before every step. False (silent): no prompts at all;
                  manual steps are skipped ("not performed (manual)").
  optional_scans  True: run optional scans (ids like 6a); False: skip them.
  confirm         callable(prompt) -> yes/no (default: input).
  report_dir      folder for the JSON report (unique file per run).

Behaviour
  declined required step -> sequence cancelled, later steps not run
  AlignmentFailure       -> step failed, sequence stops
  OperatorAbort          -> step cancelled, sequence stops
  other exception / RequestAbort / RequestStop -> partial report, re-raised
  No retries and no restoration moves are made by the framework.
Statuses: completed failed cancelled interrupted skipped not_run
"""


def beamline_setup_usage(printer=print):
    """Print the usage text. Returns it."""
    printer(USAGE)
    return USAGE
