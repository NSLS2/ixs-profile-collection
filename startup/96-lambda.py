import time

from ophyd.areadetector.cam import Lambda750kCam
from ophyd import Component as Cpt
from ophyd import ROIPlugin, TransformPlugin, EpicsSignal, EpicsSignalRO
from ophyd.areadetector.base import EpicsSignalWithRBV as SignalWithRBV
from ophyd.areadetector.cam import CamBase
from ophyd.areadetector import ADComponent as ADCpt, DetectorBase
from nslsii.ad33 import StatsPluginV33, SingleTriggerV33
from ophyd.areadetector.plugins import PluginBase


class PluginCV(PluginBase):
    comp_vision_function1 = ADCpt(EpicsSignal, 'CompVisionFunction1')
    input1 = ADCpt(EpicsSignal, 'Input1')

class LambdaDetector(DetectorBase):
    _html_docs = ['lambda.html']
    cam = Cpt(Lambda750kCam, 'cam1:')

class Lambda(SingleTriggerV33, LambdaDetector):
    # MR20200122: created all dirs recursively in /nsls2/jpls/data/lambda/
    # from 2020 to 2030 with 777 permissions, owned by xf12id1 user.
    # tiff = Cpt(TIFFPluginWithFileStore,
    #            suffix="TIFF1:",
    #            write_path_template="/nsls2/xf12id1g/data/lambda/%Y/%m/%d/",
    #            read_path_template="/nsls2/xf12id1g/data/lambda/%Y/%m/%d/",
    #            root='/nsls2/xf12id1g/data') 
    cv1 = Cpt(PluginCV, 'CV1:')
    
    roi1 = Cpt(ROIPlugin, 'ROI1:')
    roi2 = Cpt(ROIPlugin, 'ROI2:')
    roi3 = Cpt(ROIPlugin, 'ROI3:')
    roi4 = Cpt(ROIPlugin, 'ROI4:')
    roi5 = Cpt(ROIPlugin, 'ROI5:')
    roi6 = Cpt(ROIPlugin, 'ROI6:')
    roi7 = Cpt(ROIPlugin, 'ROI7:')


    stats1 = Cpt(StatsPluginV33, 'Stats1:', read_attrs=['total'])
    stats2 = Cpt(StatsPluginV33, 'Stats2:', read_attrs=['total'])
    stats3 = Cpt(StatsPluginV33, 'Stats3:', read_attrs=['total'])
    stats4 = Cpt(StatsPluginV33, 'Stats4:', read_attrs=['total'])
    stats5 = Cpt(StatsPluginV33, 'Stats5:', read_attrs=['total'])
    stats6 = Cpt(StatsPluginV33, 'Stats6:', read_attrs=['total'])
    stats7 = Cpt(StatsPluginV33, 'Stats7:', read_attrs=['total'])


    trans1 = Cpt(TransformPlugin, 'Trans1:')

    low_thr = Cpt(EpicsSignal, 'cam1:LowEnergyThreshold')
    hig_thr = Cpt(EpicsSignal, 'cam1:HighEnergyThreshold')
    oper_mode = Cpt(EpicsSignal, 'cam1:OperatingMode')
    sysreset = Cpt(EpicsSignal, ':SYSRESET')

lambda_det = Lambda('XF:10IDC-BI{Lambda-Cam:1}', name='lambda_det')
for j in range(1, 8):
    getattr(lambda_det, f'stats{j}').kind = 'normal'
lambda_det.stats7.total.kind = 'hinted'


# Impose Stats4 to be ROI4 if in the future we need to exclude bad pixels
def set_defaut_stat_roi():
    yield from bps.mv(lambda_det.stats1.nd_array_port, 'ROI1')
    yield from bps.mv(lambda_det.stats2.nd_array_port, 'ROI2')
    yield from bps.mv(lambda_det.stats3.nd_array_port, 'ROI3')
    yield from bps.mv(lambda_det.stats4.nd_array_port, 'ROI4')


def set_lambda_exposure(exposure):
    # Sets the Lambda detector exposure time (exposure)
    det = lambda_det
    yield from bps.mv(det.cam.acquire_time, exposure, det.cam.acquire_period, exposure)


# Timing parameters for the acquisition-verification steps in
# setup_lambda_detector(). These exist because IOC/plugin-derived readbacks
# such as cam1:NumImagesCounter_RBV and cam1:ArrayCounter can lag the
# completion of an acquisition by up to ~1 s. A single read after a fixed
# delay is therefore prone to false failures; the verification steps instead
# poll for a stable/increased value within a bounded timeout.
#
# Do not reduce these values merely to make tests run faster -- tests should
# monkeypatch these module-level names instead.
_LAMBDA_VERIFY_SETTLE = 0.5
_LAMBDA_VERIFY_TIMEOUT = 5.0
_LAMBDA_VERIFY_POLL = 0.2
_LAMBDA_VERIFY_CONSECUTIVE = 2
_LAMBDA_DUAL_THRESHOLD_ATTEMPTS = 3
_LAMBDA_REBOOT_TIMEOUT = 30.0
_LAMBDA_REBOOT_POLL = 0.2
_LAMBDA_REBOOT_CONSECUTIVE = 2


def _lambda_reboot_signals(det):
    """Return the signals that must be usable before setup can start."""
    return (
        det.sysreset,
        det.cam.acquire,
        det.cam.acquire_time,
        det.cam.acquire_period,
        det.cam.operating_mode,
        det.cam.detector_state,
        det.cam.num_images_counter,
        det.cam.array_counter,
        det.low_thr,
        det.hig_thr,
        det.cv1.comp_vision_function1,
        det.cv1.input1,
    )


def _reboot_lambda_ioc(
    det,
    *,
    timeout=_LAMBDA_REBOOT_TIMEOUT,
    poll_interval=_LAMBDA_REBOOT_POLL,
    consecutive=_LAMBDA_REBOOT_CONSECUTIVE,
):
    """Reboot the Lambda IOC and wait for all setup-critical PVs.

    A real disconnect must be observed before reconnection is accepted. All
    waiting is performed with Bluesky sleeps so the RunEngine remains
    responsive while Channel Access reconnects in the background.
    """
    signals = _lambda_reboot_signals(det)
    initially_unavailable = [signal.name for signal in signals if not signal.connected]
    if initially_unavailable:
        raise RuntimeError(
            "Lambda IOC reboot failed: signals unavailable before reset: "
            + ", ".join(initially_unavailable)
        )

    print(f"Rebooting Lambda IOC via {det.sysreset.pvname}")
    det.sysreset.put(1)

    deadline = time.monotonic() + timeout
    while det.cam.detector_state.connected:
        if time.monotonic() >= deadline:
            raise RuntimeError(
                f"Lambda IOC reboot failed: no disconnect observed within "
                f"{timeout} s"
            )
        yield from bps.sleep(poll_interval)
    print("Lambda IOC disconnected; waiting for critical PVs to reconnect")

    deadline = time.monotonic() + timeout
    unavailable = []
    while True:
        unavailable = [signal.name for signal in signals if not signal.connected]
        if not unavailable:
            unreadable = []
            for signal in signals:
                try:
                    signal.get()
                except Exception:
                    unreadable.append(signal.name)
            unavailable.extend(unreadable)
        if not unavailable:
            break
        if time.monotonic() >= deadline:
            raise RuntimeError(
                f"Lambda IOC reboot failed: IOC did not reconnect within "
                f"{timeout} s; unavailable signals: {', '.join(unavailable)}"
            )
        yield from bps.sleep(poll_interval)

    deadline = time.monotonic() + timeout
    idle_matches = 0
    last_state = None
    while True:
        last_state = det.cam.detector_state.get(as_string=True)
        if last_state == "Idle":
            idle_matches += 1
            if idle_matches >= consecutive:
                print("Lambda IOC is online and detector state is Idle")
                return
        else:
            idle_matches = 0
        if time.monotonic() >= deadline:
            raise RuntimeError(
                f"Lambda IOC reboot failed: detector did not become Idle "
                f"within {timeout} s; last state: {last_state}"
            )
        yield from bps.sleep(poll_interval)


def _wait_for_stable_value(
    signal,
    expected,
    *,
    timeout=_LAMBDA_VERIFY_TIMEOUT,
    poll_interval=_LAMBDA_VERIFY_POLL,
    consecutive=_LAMBDA_VERIFY_CONSECUTIVE,
    label=None,
):
    """Bluesky plan: poll `signal` until it reads `expected` for `consecutive`
    consecutive reads, or raise RuntimeError after `timeout` seconds.

    Any read that differs from `expected` resets the consecutive-match
    counter. This tolerates asynchronous EPICS/plugin propagation delay
    (e.g. NumImagesCounter_RBV lagging acquisition completion) while still
    requiring the value to have genuinely settled rather than transiently
    matching once.

    Never blocks the RunEngine with time.sleep(); all waiting is done via
    `yield from bps.sleep(...)`. `time.monotonic()` is used only to bound
    the total elapsed wall-clock time.
    """
    label = label or signal.name
    deadline = time.monotonic() + timeout
    matches = 0
    last_value = None

    while True:
        last_value = signal.get()
        if last_value == expected:
            matches += 1
            if matches >= consecutive:
                return last_value
        else:
            matches = 0

        if time.monotonic() >= deadline:
            raise RuntimeError(
                f"{label}: expected stable value {expected}; "
                f"last observed value {last_value} after {timeout} s"
            )

        yield from bps.sleep(poll_interval)


def _wait_for_increase(
    signal,
    baseline,
    *,
    timeout=_LAMBDA_VERIFY_TIMEOUT,
    poll_interval=_LAMBDA_VERIFY_POLL,
    consecutive=_LAMBDA_VERIFY_CONSECUTIVE,
    label=None,
):
    """Bluesky plan: poll `signal` until it reads a value greater than
    `baseline` for `consecutive` consecutive reads, or raise RuntimeError
    after `timeout` seconds.

    The condition checked on each read is simply `value > baseline` (not
    `value > previous_reading`), so e.g. baseline=100 with observed reads
    100, 101, 101 satisfies consecutive=2 -- the signal does not need to
    keep incrementing on every poll, only to have moved past baseline and
    stayed there.

    Never blocks the RunEngine with time.sleep(); all waiting is done via
    `yield from bps.sleep(...)`.
    """
    label = label or signal.name
    deadline = time.monotonic() + timeout
    matches = 0
    last_value = None

    while True:
        last_value = signal.get()
        if last_value > baseline:
            matches += 1
            if matches >= consecutive:
                return last_value
        else:
            matches = 0

        if time.monotonic() >= deadline:
            raise RuntimeError(
                f"{label}: expected value greater than {baseline}; "
                f"last observed value {last_value} after {timeout} s"
            )

        yield from bps.sleep(poll_interval)


def _start_lambda_acquisition(det):
    """Start one acquisition and allow its readbacks to begin updating.

    The Lambda Acquire readback does not provide reliable completion
    semantics for EpicsSignal.set(1), so initiate the acquisition with a
    direct EPICS write. The bounded NumImagesCounter and ArrayCounter
    polling in setup_lambda_detector() determines whether a new acquisition
    actually completed.
    """
    det.cam.acquire.put(1)
    yield from bps.sleep(_LAMBDA_VERIFY_SETTLE)


def _wait_for_stable_acquisition_result(
    det,
    baseline,
    *,
    timeout=_LAMBDA_VERIFY_TIMEOUT,
    poll_interval=_LAMBDA_VERIFY_POLL,
    consecutive=_LAMBDA_VERIFY_CONSECUTIVE,
    label=None,
):
    """Wait for new array data and a stable pair of acquisition counters.

    A result is complete when ArrayCounter is greater than `baseline` and
    the (NumImagesCounter, ArrayCounter) pair is unchanged for `consecutive`
    reads. This distinguishes a completed one-image acquisition, which can
    be retried promptly, from an acquisition that has not produced data yet.
    """
    label = label or "Lambda acquisition"
    deadline = time.monotonic() + timeout
    previous = None
    matches = 0
    last_result = (det.cam.num_images_counter.get(), det.cam.array_counter.get())

    while True:
        last_result = (
            det.cam.num_images_counter.get(),
            det.cam.array_counter.get(),
        )
        if last_result[1] > baseline:
            if last_result == previous:
                matches += 1
            else:
                previous = last_result
                matches = 1
            if matches >= consecutive:
                return last_result
        else:
            previous = None
            matches = 0

        if time.monotonic() >= deadline:
            raise RuntimeError(
                f"{label}: no stable acquisition result after {timeout} s; "
                f"last NumImagesCounter={last_result[0]}, "
                f"ArrayCounter={last_result[1]}, baseline={baseline}"
            )

        yield from bps.sleep(poll_interval)


def _verify_dual_threshold_acquisition(
    det,
    *,
    image_step,
    array_step,
    attempts=_LAMBDA_DUAL_THRESHOLD_ATTEMPTS,
):
    """Acquire until one attempt produces two images and new array data.

    The first acquisition after a Lambda mode/plugin change can retain the
    previous one-image behavior. Each attempt is fully verified before a
    retry. A stable one-image result is rejected promptly rather than waiting
    for the full timeout; no ArrayCounter movement still waits until timeout.
    """
    last_error = None

    for attempt in range(1, attempts + 1):
        counter_before = det.cam.array_counter.get()
        step = image_step.removesuffix(" failed")
        yield from _start_lambda_acquisition(det)

        try:
            images, array_counter = yield from _wait_for_stable_acquisition_result(
                det,
                counter_before,
                timeout=_LAMBDA_VERIFY_TIMEOUT,
                poll_interval=_LAMBDA_VERIFY_POLL,
                consecutive=_LAMBDA_VERIFY_CONSECUTIVE,
                label=f"{step} attempt {attempt}/{attempts}",
            )
        except RuntimeError as exc:
            if det.cam.array_counter.get() <= counter_before:
                last_error = RuntimeError(
                    f"{array_step}: expected ArrayCounter greater than "
                    f"{counter_before}; last observed value "
                    f"{det.cam.array_counter.get()} after "
                    f"{_LAMBDA_VERIFY_TIMEOUT} s"
                )
            else:
                last_error = exc
            print(
                f"{step}: acquisition attempt {attempt}/{attempts} did not "
                f"produce new stable counter data; retrying. "
                f"Details: {last_error}"
            )
            continue

        if images == 2 and array_counter > counter_before:
            print(
                f"{step}: acquisition attempt {attempt}/{attempts} accepted "
                f"with {images} images and ArrayCounter increase "
                f"{array_counter - counter_before}"
            )
            return

        last_error = RuntimeError(
            f"{image_step}: expected 2 images with new array data; observed "
            f"NumImagesCounter={images}, ArrayCounter={array_counter}, "
            f"baseline={counter_before}"
        )
        print(
            f"{step}: acquisition attempt {attempt}/{attempts} produced "
            f"{images} image(s), ArrayCounter increase "
            f"{array_counter - counter_before}; retrying"
        )

    raise RuntimeError(
        f"{image_step}: no valid two-image acquisition after {attempts} "
        f"attempts; last failure: {last_error}"
    ) from last_error


def setup_lambda_detector():
    """Configure the Lambda detector for standard use.

    Applies settings and verifications in the order of the standard setup procedure:

    Phase 1 - ContinuousReadWrite mode (steps 1-6):
        AcquireTime=1s, AcquirePeriod=1s, LowEnergyThreshold=4.5 keV,
        OperatingMode=ContinuousReadWrite. Verifies 1 image per acquisition
        and that ArrayCounter increments.

    Phase 2 - DualThreshold mode (steps 7-10):
        OperatingMode=DualThreshold, LowEnergyThreshold=4.5 keV,
        HighEnergyThreshold=11.0 keV. Verifies 2 images per acquisition.

    Phase 3 - ADCompVision + final verification (steps 12-16):
        CV1:CompVisionFunction1=Subtract, CV1:Input1=1.
        Verifies 2 images per acquisition and that ArrayCounter increments.

    Verification steps that depend on IOC/plugin-derived readbacks
    (NumImagesCounter_RBV, ArrayCounter) poll for a stable or increased
    value within a bounded timeout (see _wait_for_stable_value and
    _wait_for_increase) rather than reading once after a fixed delay, since
    those readbacks can lag acquisition completion by up to ~1 s.

    Raises RuntimeError if any verification step fails to reach/stabilize at
    the expected value within the timeout.
    """
    det = lambda_det

    yield from _reboot_lambda_ioc(det)

    # Step 1
    yield from bps.mv(det.cam.acquire_time, 1)
    # Step 2
    yield from bps.mv(det.cam.acquire_period, 1)
    # Step 3
    yield from bps.mv(det.low_thr, 4.5)
    # Step 5
    yield from bps.mv(det.cam.operating_mode, 'ContinuousReadWrite')

    # Steps 4 and 6 verify the same ContinuousReadWrite acquisition.
    counter_before = det.cam.array_counter.get()
    yield from _start_lambda_acquisition(det)

    # Step 4: verify 1 image per acquisition in ContinuousReadWrite mode
    yield from _wait_for_stable_value(
        det.cam.num_images_counter,
        1,
        timeout=_LAMBDA_VERIFY_TIMEOUT,
        poll_interval=_LAMBDA_VERIFY_POLL,
        consecutive=_LAMBDA_VERIFY_CONSECUTIVE,
        label=f"Step 4 failed: {det.cam.num_images_counter.name}",
    )

    # Step 6: verify ArrayCounter activity from the Step 4 acquisition
    yield from _wait_for_increase(
        det.cam.array_counter,
        counter_before,
        timeout=_LAMBDA_VERIFY_TIMEOUT,
        poll_interval=_LAMBDA_VERIFY_POLL,
        consecutive=_LAMBDA_VERIFY_CONSECUTIVE,
        label=f"Step 6 failed: {det.cam.array_counter.name}",
    )
    print(
        "Phase 1 verified: images=1, "
        f"ArrayCounter +{det.cam.array_counter.get() - counter_before}"
    )

    # Step 7
    yield from bps.mv(det.cam.operating_mode, 'DualThreshold')
    # Step 8
    yield from bps.mv(det.low_thr, 4.5)
    # Step 9
    yield from bps.mv(det.hig_thr, 11.0)

    # Step 10: verify one new 2-image acquisition in DualThreshold mode.
    yield from _verify_dual_threshold_acquisition(
        det,
        image_step="Step 10 failed",
        array_step="Step 10 ArrayCounter failed",
    )

    # Step 12
    yield from bps.mv(det.cv1.comp_vision_function1, 'Subtract')
    # Step 13
    yield from bps.mv(det.cv1.input1, 1)

    # Steps 15 and 16 verify the same post-CV acquisition.
    yield from _verify_dual_threshold_acquisition(
        det,
        image_step="Step 15 failed",
        array_step="Step 16 failed",
    )
    print("Lambda setup completed successfully")
