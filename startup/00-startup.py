###############################################################################
# TODO: remove this block once https://github.com/bluesky/ophyd/pull/959 is
# merged/released.
from datetime import datetime
from ophyd.signal import EpicsSignalBase, EpicsSignal, DEFAULT_CONNECTION_TIMEOUT

def print_now():
    return datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S.%f')

def wait_for_connection_base(self, timeout=DEFAULT_CONNECTION_TIMEOUT):
    '''Wait for the underlying signals to initialize or connect'''
    if timeout is DEFAULT_CONNECTION_TIMEOUT:
        timeout = self.connection_timeout
    # print(f'{print_now()}: waiting for {self.name} to connect within {timeout:.4f} s...')
    start = time.time()
    try:
        self._ensure_connected(self._read_pv, timeout=timeout)
        # print(f'{print_now()}: waited for {self.name} to connect for {time.time() - start:.4f} s.')
    except TimeoutError:
        if self._destroyed:
            raise DestroyedError('Signal has been destroyed')
        raise

def wait_for_connection(self, timeout=DEFAULT_CONNECTION_TIMEOUT):
    '''Wait for the underlying signals to initialize or connect'''
    if timeout is DEFAULT_CONNECTION_TIMEOUT:
        timeout = self.connection_timeout
    # print(f'{print_now()}: waiting for {self.name} to connect within {timeout:.4f} s...')
    start = time.time()
    self._ensure_connected(self._read_pv, self._write_pv, timeout=timeout)
    # print(f'{print_now()}: waited for {self.name} to connect for {time.time() - start:.4f} s.')

EpicsSignalBase.wait_for_connection = wait_for_connection_base
EpicsSignal.wait_for_connection = wait_for_connection
###############################################################################

from ophyd.signal import EpicsSignalBase
# EpicsSignalBase.set_default_timeout(timeout=10, connection_timeout=10)  # old style
EpicsSignalBase.set_defaults(timeout=10, connection_timeout=10)  # new style

import nslsii
nslsii.configure_base(get_ipython().user_ns, 'ixs')

# After the above call, you will now have the following in your namespace:
#
#	RE : RunEngine
#	db : databroker
#	sd : SupplementalData
#	pbar_manager : ProgressBarManager
#	bec : BestEffortCallback
#	peaks : bec.peaks
#	plt : matplotlib.pyplot
#	np : numpy
#	bc : bluesky.callbacks
#	bp : bluesky.plans
#	bps : bluesky.plan_stubs
#	mv : bluesky.plan_stubs.mv
#	mvr : bluesky.plan_stubs.mvr
#	mov : bluesky.plan_stubs.mov
#	movr : bluesky.plan_stubs.movr
#	bpp : bluesky.preprocessors


# At the end of every run, verify that files were saved and
# print a confirmation message.
from bluesky.callbacks.broker import verify_files_saved
# RE.subscribe(post_run(verify_files_saved), 'stop')


# Optional: set any metadata that rarely changes.
# RE.md['beamline_id'] = 'YOUR_BEAMLINE_HERE'

# Uncomment the following lines to turn on verbose messages for
# debugging.
# import logging
# ophyd.logger.setLevel(logging.DEBUG)
# logging.basicConfig(level=logging.DEBUG)

# New figure title so no overplot.
def relabel_fig(fig, new_label):
    fig.set_label(new_label)
    fig.canvas.manager.set_window_title(fig.get_label())

# Implementing grid on plots:
import matplotlib as mpl
mpl.rcParams['axes.grid'] = True

RE.md['beamline_id'] = 'IXS'
RE.md['owner'] = 'xf10id'
RE.md['group'] = 'ixs'

from pyOlog.ophyd_tools import get_all_positioners

def relabel_motors():
    mtrs = get_all_positioners()
    for mtr in mtrs:
        attr = getattr(mtr, 'user_readback')
        attr.name = mtr.name
