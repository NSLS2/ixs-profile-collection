from ophyd.signal import EpicsSignalBase

EpicsSignalBase.set_defaults(timeout=10, connection_timeout=10)  # new style

import nslsii

nslsii.configure_base(get_ipython().user_ns, "ixs")

# After the above call, you will now have the following in your namespace:
#
# 	RE : RunEngine
# 	db : databroker
# 	sd : SupplementalData
# 	pbar_manager : ProgressBarManager
# 	bec : BestEffortCallback
# 	peaks : bec.peaks
# 	plt : matplotlib.pyplot
# 	np : numpy
# 	bc : bluesky.callbacks
# 	bp : bluesky.plans
# 	bps : bluesky.plan_stubs
# 	mv : bluesky.plan_stubs.mv
# 	mvr : bluesky.plan_stubs.mvr
# 	mov : bluesky.plan_stubs.mov
# 	movr : bluesky.plan_stubs.movr
# 	bpp : bluesky.preprocessors


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

mpl.rcParams["axes.grid"] = True

RE.md["beamline_id"] = "IXS"
RE.md["owner"] = "xf10id"
RE.md["group"] = "ixs"

from pyOlog.ophyd_tools import get_all_positioners


def relabel_motors():
    mtrs = get_all_positioners()
    for mtr in mtrs:
        attr = getattr(mtr, "user_readback")
        attr.name = mtr.name

# ## Live specfile exporting
import time
from event_model import RunRouter
from suitcase.specfile import Serializer


def spec_factory(name, doc):
    if not spec_factory.enabled:
        return [], []
    directory = "/nsls2/data/ixs/legacy/specfiles/"
    file_prefix = "spec_" + time.strftime("%Y-%m-%d")
    spec_cb = Serializer(directory, file_prefix=file_prefix, flush=True)
    return [spec_cb], []


spec_factory.enabled = True

spec_router = RunRouter([spec_factory])
RE.subscribe(spec_router)
