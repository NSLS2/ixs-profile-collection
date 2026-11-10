# MCM Alignment Plan: Project Scope

## 1. Purpose and status

Develop a Bluesky plan to align the MCM using horizontal and vertical intensity scans, building on the first-attempt code:

- `mcm_setup`
- `mcm_setup_prep`
- `mcm_setup_post`

These exist in `startup/97-setup_plans.py`.

**Status: updated requirements for review. No authorization to change implementation or operate beamline hardware.**

## 2. Working rules

- Proceed step by step; propose and obtain approval before implementation changes.
- Inspect the existing code before relying on reported device mappings or behavior.
- Do not invent device names, settings, correction formulas, or undocumented units.
- Reuse existing fitting and HRM-state detection procedures where suitable.
- Verify using simulated devices before live operation.
- Keep numerical settings configurable, using the confirmed values below as initial defaults.
- Execute device changes through Bluesky and await completion.

## 3. Devices and reference procedures

| Function | Device or reference |
|---|---|
| Analyzer height | `anc_xtal.y`, mm |
| Analyzer slit blades | `analyzer_slits.top/bottom/outboard/inboard` |
| Photodiode position | `anpd` |
| Absorption wheel | `whl` |
| Air pad | `airpad` |
| Detector measurement range (sensitivity: 0 less sensitive, 1 more sensitive) | `det2.em_range` |
| Intensity | `det2.current1.mean_value` |
| MCM horizontal position | `mcm.x`, mm |
| MCM vertical position | `mcm.y`, mm |
| Sample-tower horizontal correction | `sample_stage.tx`, mm |
| HRM angular correction | `hrm2.dch`, degrees |
| HRM vertical correction | `hrm2.dy`, mm |
| Spectrometer angle | `spec.tth` |
| HRM-state detection | `hrm2.ux` and `hrm2.dx` readbacks; see below |
| Peak fitting | `calculate_max_value` (`startup/97-setup_plans.py`), used in `mcm_setup()` |

Units of `hrm2.dch` (degrees), `sample_stage.tx` (mm) and `hrm2.dy` (mm) are confirmed.

### HRM-state detection

The state is determined automatically from the `hrm2.ux` and `hrm2.dx` readbacks:

| `hrm2.ux` and `hrm2.dx` | State |
|---|---|
| Both within 0 ± 2 mm | HRM in |
| Both within −20 ± 2 mm | HRM out |
| Anything else | Inconsistent |

"Anything else" includes one axis in and one out, either axis at an intermediate position, and a failed readback. An inconsistent state stops the plan before any preparation move and produces diagnostics, including both readbacks, the thresholds and the state.

The state is read once, at the start of the plan. The HRM is not expected to move during alignment.

The existing `_read_hrm_state()` (`startup/98-alignment.py:73`) treats every non-"in" case as "out" and has no inconsistent state, so the plan needs its own three-state check. `98-alignment.py` is not modified.

## 4. Options and initial defaults

| Parameter | Default |
|---|---:|
| `restore_positions` | `False` |
| `interactive` | `True` |
| `align_horizontal` | `True` |
| `align_vertical` | `True` |
| `MCM_XPOS` | −0.941 mm |
| Analyzer lowering displacement | −4 mm |
| Allowed lowered analyzer position | Strictly between 0 and 1 mm |
| Scan half-range, either axis | 0.2 mm |
| Points per scan, either axis | 41 |
| Acquisition time per point | 1 s |
| Horizontal alignment tolerance | Absolute error < 0.002 mm |
| HRM-in vertical alignment tolerance | Absolute error < 0.002 mm |
| Maximum scans per iterative axis | 4 |
| Low-intensity threshold | `1e5` |
| Maximum sample-tower correction magnitude | 0.02 mm |
| Maximum HRM vertical correction magnitude | 0.2 mm |
| Maximum HRM angular correction term magnitude | 0.008° |

The correction factors, two-step angular offsets, tolerances, scan settings, and limits must be configurable.

No additional 0.04 mm correction limit applies to `mcm.x` or `mcm.y`; that proposal was withdrawn.

## 5. Preparation

Save the entry values required for restoration before changing them:

- Analyzer height.
- Four analyzer slit positions.
- Photodiode position.
- Absorption-wheel setting.
- Detector measurement range.

For HRM-in alignment, capture the current `mcm.y` as `MCM_YPOS` at the start and hold it fixed throughout the procedure. This position is expected to have been established previously under HRM-out conditions.

Perform preparation:

0. Determine the HRM state (section 3). An inconsistent state stops the plan before any move.
1. Handle analyzer height:
   - If `0 < anc_xtal.y < 1 mm`, leave it unchanged.
   - Otherwise calculate `y_target = y_initial - 4 mm`.
   - Require `0 < y_target < 1 mm` before moving.
   - After moving, verify that the readback is also strictly within this interval.
   - Approximately 0.5 mm is an expected result, not an absolute target.
2. Set analyzer slit positions:
   - Top: +2 mm.
   - Bottom: −2 mm.
   - Outboard: +2 mm.
   - Inboard: −2 mm.
3. Move `anpd` to 0.
4. Apply the HRM-dependent measurement settings:

| Setting | HRM in | HRM out |
|---|---:|---:|
| `whl` | 0 | 2 |
| `det2.em_range` | 1 | 0 |

5. Set `airpad` to 1 and await completion. Then move `spec.tth` to 0. Exclude `spec.tth` and `airpad` from restoration.
6. Await completion of preparation settings before acquiring intensity.
7. Check `det2.current1.mean_value` once after preparation and before scanning.
   - If below `1e5`, stop with diagnostics and failure restoration.
   - Do not apply this threshold to individual scan points or fitted peak intensities.

If any preparation move fails, abort the alignment. This is an alignment failure and triggers failure restoration (section 10). If the analyzer-height check fails before any move, nothing has changed and nothing needs restoring.

## 6. Horizontal alignment

If enabled, scan `mcm.x` around `MCM_XPOS`, using the configured range, point count, and acquisition time.

`MCM_XPOS` is fixed, and `mcm.x` must be at `MCM_XPOS` at all times outside a scan. Before every horizontal scan, move `mcm.x` to `MCM_XPOS` using the existing motion interface (`bps.mv`); no additional done-signal or position verification is added. The existing `dscan` is relative, so it is centered on this position.

Use `calculate_max_value` with the exact scan UID to determine `xmax`. The scan UID is captured from the scan's Start document and used for the fit (not `uid=-1`). After the scan and before retrieving that UID, flush the Tiled writer (`flush_tiled()`, `startup/02-tiled_writer.py`). Its completion behavior (whether the run is readable after the flush) must be verified in the live environment.

After each horizontal scan, return `mcm.x` to `MCM_XPOS`, default −0.941 mm unless explicitly overridden.

Calculate:

```python
dxmax = xmax - MCM_XPOS
```

Convergence requires:

```python
abs(dxmax) < horizontal_tolerance
```

### HRM out

Apply a relative sample-tower correction:

```python
sample_stage_tx_correction = dxmax  # mm
```

Require:

```python
abs(sample_stage_tx_correction) <= 0.02  # mm
```

### HRM in

Define the angular correction term:

```python
angular_term = dxmax / 10  # degrees, with dxmax in mm
```

Require:

```python
abs(angular_term) <= 0.008  # degrees
```

If a correction is needed:

| Condition | Relative `hrm2.dch` moves |
|---|---|
| `dxmax < 0` | `-dxmax / 10` degrees |
| `dxmax > 0` | First −0.01°, await completion, then `0.009 - dxmax / 10` degrees |

The 0.008° limit applies to `abs(dxmax / 10)`. It does not apply separately to the fixed two-step offsets or their net displacement. Configured motor travel limits still apply to every target.

### Iteration

- Allow at most four scans.
- After scans 1–3, apply a valid correction if outside tolerance, then rescan.
- The fourth scan must confirm convergence.
- If the fourth scan is outside tolerance, apply no further correction; report non-convergence and perform failure restoration.
- Do not apply corrections after convergence.

## 7. Vertical alignment

After horizontal alignment succeeds or is skipped, offer vertical alignment in interactive mode.

Use 41 points, ±0.2 mm, and 1 s per point by default. Determine `ymax` using `calculate_max_value`.

### HRM out

- Save `mcm.y` immediately before the vertical scan.
- Center the scan on that position.
- After a successful fit, move `mcm.y` to `ymax`.
- Complete this axis without repeating the scan.
- On failure, return `mcm.y` to its saved pre-scan position.

### HRM in

- Center every vertical scan on the fixed `MCM_YPOS` captured at the start.
- Return `mcm.y` to `MCM_YPOS` after each scan.
- Calculate:

```python
dymax = ymax - MCM_YPOS
hrm_dy_correction = -5 * dymax  # mm
```

- Require:

```python
abs(hrm_dy_correction) <= 0.2  # mm
```

- Convergence requires:

```python
abs(dymax) < vertical_tolerance
```

- Use the same maximum of four scans as horizontal alignment.
- The fourth scan must confirm convergence; do not apply another correction if it fails.

After vertical completion, finish without another horizontal scan.

## 8. Interaction and silent operation

### Interactive mode

- Ask whether to perform horizontal alignment before starting that axis.
- Ask whether to perform vertical alignment before starting that axis.
- Declining an axis means "skipped by operator" and counts as successful completion of that stage.
- Both axes may be skipped.
- Request confirmation before every alignment correction, including the HRM-out move of `mcm.y` to its fitted maximum.
- Show the motor, current position, proposed relative displacement, and target.
- One confirmation authorizes both moves of a two-step `hrm2.dch` correction.
- Declining a correction is an operator abort, not a successful skip.

### Silent mode

With `interactive=False`:

- Issue no prompts.
- Use `align_horizontal` and `align_vertical` to select axes.
- Apply valid corrections automatically.
- Continue producing results and failure diagnostics.

## 9. Fit validation and correction limits

- A failed fit stops alignment immediately.
- A fit fails when `calculate_max_value` raises `ValueError` because the sampled maximum is at the first or last scan point.
- `calculate_max_value` is used unchanged. There is no check for inverted curvature or a flat curve. For a non-concave three-point window (for example a tie at the maximum) it can return a window-edge value without error. This is a known and accepted limitation.
- A flat-curve (minimum contrast) check is a possible new requirement. It is pending definition and is not part of the current scope.
- Use a large `sampling` value (for example 200) so the resampling step is well below the tolerance. Resampling improves numerical resolution only, not measurement accuracy. `mcm_setup` currently uses `sampling=5`, which quantizes the result to 0.005 mm.
- Comparing the fitted maximum with `PeakStats.max` may be used as a consistency check. It does not establish run identity.
- Do not silently substitute a raw maximum for a failed fit.
- Check correction magnitude and motor travel limits before applying a correction.
- Validate both targets before starting a two-step angular correction.
- If a correction exceeds a limit, fail without applying it.
- Do not clip corrections to their limits.

The absolute correction limits replace the earlier proposed 5% rule.

## 10. Completion, failure, and interruption

| Outcome | Behavior |
|---|---|
| Success, `restore_positions=False` | Retain the measurement configuration, including the lowered analyzer |
| Success, `restore_positions=True` | Restore saved analyzer height, slit positions, `anpd`, `whl`, and detector range |
| Alignment failure | Restore saved analyzer height, slit positions, `anpd`, `whl`, and detector range regardless of the restoration option |
| Operator interruption, abort, or declined correction | Leave positions and settings as they are; perform no restoration |

Additional rules:

- Preserve sample-tower and HRM corrections already applied, including after failure.
- On failure, return `mcm.x` to `MCM_XPOS`.
- On HRM-in failure, return `mcm.y` to the fixed `MCM_YPOS`.
- On HRM-out vertical failure, return `mcm.y` to its pre-scan position.
- Keep `spec.tth` at 0 once set; never restore its entry value. `airpad` is also not restored.
- Successful axis skips use normal success restoration behavior.
- Distinguish operator interruption from failure so interruption cannot accidentally trigger restoration.
- Report restoration failures without hiding the original failure.

## 11. Diagnostics and results

Every failure must produce diagnostics, including in silent mode:

- Failed stage and reason.
- HRM state.
- Axis and scan iteration, where applicable.
- Available scan identifier.
- Available fitted maximum, reference position, and alignment error.
- Proposed or last applied correction.
- Relevant limit or threshold.
- Restoration outcome and any restoration errors.

On completion, report:

- Overall status.
- Per-axis status: aligned or skipped.
- Scan maxima and correction history.
- Final positions.
- Whether preparation settings were restored or retained.

## 12. Existing-code review

Findings from the repository, checked against the previously reported issues:

- **Missing branches:** `mcm_setup` implements only the horizontal scan. `s2` (vertical) is accepted but does nothing, and there is no HRM-in branch.
- **Correction loop:** The horizontal correction loop is commented out. `mcm_setup` computes `dxmax = MCM_XPOS - xmax`, prints it, and moves nothing. The scope's `xmax - MCM_XPOS` supersedes this sign.
- **Bare return:** The low-intensity abort in `mcm_setup_prep` uses a bare `return` (`None`), so `mcm_setup` fails with a `TypeError` on `aret[0]`.
- **Restoration:** `mcm_setup_post` is never called. It restores `anc_xtal.y`, `whl` (to 0) and `anpd` (to -90), and sets the slits to fixed values (1, -1, 1.5, -1.5), not to the saved entry values.
- **Device changes outside the RunEngine:** `airpad.set(1)` and `det2.em_range.set(0)` in `mcm_setup_prep` are not awaited. `98-alignment.py` already sets `det2.em_range` through `bps.mv` in `_ccr_temporary_settings_wrapper`, which also restores through `bpp.finalize_wrapper`.
- **MCM move status:** `MCMBase.set` (`startup/10-optics.py:112-148`) reports success after a 3 s `MOVE_TIMEOUT` even if `InPos` was not reached. Moves of `mcm.x` and `mcm.y` need a readback check.
- **HRM-state thresholds:** `mcm_setup_prep` treats `hrm2.ux > -5` or `hrm2.dx > -5` as in. `_read_hrm_state()` uses ±2 mm of 0 for in, and everything else is out. Neither has an inconsistent state. The three-state rule in section 3 supersedes both.
- **Wheel setting:** `_check_whl_for_mcm` (`98-alignment.py:227`) requires `whl = 2` with HRM out and `0` with HRM in. The scope now uses the same values (HRM out: `whl` 2, range 0; HRM in: `whl` 0, range 1), so `alignment_scan('mcm')` passes its `whl` precondition after this preparation. The check and `docs/alignment_tracking_status.md` are not changed.
- **Slit settings:** No conflict. `_check_slits_open` targets +2, -2, +2, -2 with a ±0.1 mm tolerance, which agrees with the preparation values. The ±0.1 is a tolerance, not an opening.
- **Fitting:** `calculate_max_value` fits a quadratic through three points around the sampled maximum (`delta=1`) and resamples at 200 points. It raises `ValueError` only when the maximum is at either end. It does not detect a flat curve, an inverted curve, or noise. It is used unchanged, with no added checks (section 9). It reads the run from `trc[uid]`, and `mcm_setup` currently calls it with `uid=-1` and `sampling=5`.
- **Detector channel mapping:** `dscan` defaults to `det_ch=[0]`, while the fit reads `det2_current1_mean_value`. Confirmed: channel 0 corresponds to `det2_current1_mean_value`.

The confirmed requirements in this scope resolve behavioral conflicts.

## 13. Deliverables and milestones

Each implementation stage requires approval before proceeding.

1. Update the repository scope document with these confirmed requirements.
2. Review the existing source and report any remaining ambiguities.
3. Propose the implementation structure and configuration interface.
4. Implement and test fitting validation, corrections, limits, and convergence without hardware.
5. Implement the Bluesky plan and restoration behavior.
6. Verify with simulated devices:
   - HRM in and out, and an inconsistent HRM state.
   - Positive and negative corrections.
   - Convergence and fourth-scan non-convergence.
   - Failed fit and out-of-range maximum.
   - Low intensity.
   - Correction limits and motor failures.
   - Success restoration enabled and disabled.
   - Failure restoration and restoration errors.
   - Operator interruption and declined corrections.
   - Axis skipping and silent operation.
7. Produce operating instructions and a proposed live commissioning procedure.

## 14. Open items

- **Flat-curve (minimum contrast) check:** pending definition. Not part of the current scope.
- **Writer flush semantics:** `flush_tiled()` raises `RuntimeError` on timeout or a dead worker thread. Whether the run is readable after it returns could not be checked in the source (`BufferingWrapper` was not inspectable). Verify with a live read-back test.
- **Settling time** after the `em_range` and `whl` changes, before the one-shot intensity check: not defined.
- **Interactive `input()`** inside a plan blocks the RunEngine loop and is not queue-server compatible: decision pending.
- **Low-intensity threshold** `1e5` is kept for both HRM states for now.

Resolved:

- `det_ch=[0]` in `dscan` corresponds to `det2_current1_mean_value`.
- `det2.em_range` is detector sensitivity: 0 is less sensitive (HRM out), 1 is more sensitive (HRM in, where intensity is low).
- `MCMBase` and its motion interface are unchanged. No extra done-signal checks or independent position verification are added.
- The analyzer is left unchanged if `0 < y_initial < 1 mm`.

## 15. Implementation status

**Milestone 1 (done):** `startup/utils/mcm_alignment.py` (hardware-free logic) and `tests/test_mcm_alignment_logic.py`.

**Milestone 2 (done, simulation only, not run on hardware):** plan layer in `startup/97-setup_plans.py` and `tests/test_mcm_setup_plan.py`. Only the old `mcm_setup_prep`, `mcm_setup_post` and `mcm_setup` were replaced. `calculate_max_value`, `98-alignment.py`, `10-optics.py` (`MCMBase`) and the alignment-tracking docs are unchanged.

Usage:

```python
RE(mcm_setup())                                              # interactive
RE(mcm_setup(interactive=False, restore_positions=True))     # silent
```

Keyword overrides set `McmAlignConfig` fields (unknown names raise `TypeError`). The previous `mcm_setup(s1, s2)` signature is replaced by `align_horizontal` and `align_vertical`. The plan returns a result dict (status `success`, `failed` or `aborted`, HRM state, axes, scan history with UIDs, corrections, final positions, restoration outcome and errors) and prints the same report.

Decisions made during implementation, for review:

- **Prompts:** `input()` is the default `ask` function, as in the existing macros, and it can be replaced (`ask=`). Prompts require an explicit `yes`; an empty answer declines. Declining an axis skips it. Declining a correction is an operator abort (no restoration).
- **Interruption:** `RequestAbort` and `RequestStop` (RunEngine abort or stop) pass through without restoration. Every other exception, including a failed device status, is an alignment failure and restores. A failed-status failure reports the underlying cause.
- **Restoration scope:** only items preparation changed are restored (analyzer height only if it was moved). Every item is attempted and errors are collected, so one failed restore does not stop the others or hide the original failure. `airpad` and `spec.tth` are never restored.
- **Scan flow:** `mcm.x` (or `mcm.y` for HRM-in) is moved to the fixed reference before every scan and after it. The scan UID is taken from the Start document, `flush_tiled()` is called after the scan, and `calculate_max_value(uid=...)` runs with `delta=1` and `sampling=200` (configurable `fit_delta`, `fit_sampling`). A failed flush or fit is an alignment failure.
- **Limits:** motor travel limits are taken from `motor.limits` when available (a (0, 0) pair means unset) and are checked together with the correction limits. Corrections are never clipped.
- **Not implemented, as agreed:** settling time, flat-curve check, inverted-curvature check, `MCMBase` done-signal or independent position verification.
- **Wheel consistency:** the `whl` values agree with `_check_whl_for_mcm` in `98-alignment.py` (HRM out 2, HRM in 0).

Simulation coverage (`tests/test_mcm_setup_plan.py`): HRM in and out, positive and negative corrections including the two-step `dch` order, convergence and fourth-scan non-convergence, fit and flush failures, low intensity, correction limits, prep move failure, restoration error not hiding the original failure, exact UID with flush before fit, interactive prompts and one confirmation for a two-step move, declined correction, RunEngine abort without restoration, axis skipping, silent mode, `restore_positions` on and off, the already-lowered analyzer, inconsistent HRM state and analyzer-interval boundaries.

Still required before any live use: the writer flush read-back check (section 14), a live commissioning procedure with operator confirmation at each hardware-moving stage, and your explicit authorization to operate hardware.
