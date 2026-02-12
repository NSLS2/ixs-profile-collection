from ophyd.signal import EpicsSignalBase

EpicsSignalBase.set_defaults(timeout=10, connection_timeout=10)  # new style

import nslsii
import os

nslsii.configure_base(
    get_ipython().user_ns,
    'ixs',
    publish_documents_with_kafka=False,
    call_returns_result=True,
    redis_url="info.ixs.nsls2.bnl.gov"
)

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
from bluesky.utils import PersistentDict

# RE.subscribe(post_run(verify_files_saved), 'stop')


# Optional: set any metadata that rarely changes.
# RE.md['beamline_id'] = 'YOUR_BEAMLINE_HERE'

# Uncomment the following lines to turn on verbose messages for
# debugging.
# import logging
# ophyd.logger.setLevel(logging.DEBUG)
# logging.basicConfig(level=logging.DEBUG)

bec.disable_plots()
# bec.enable_plots()

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
RE.md = PersistentDict('/nsls2/data/ixs/shared/config/bluesky/profile_collection/md')

from pyOlog.ophyd_tools import get_all_positioners


def relabel_motors():
    mtrs = get_all_positioners()
    for mtr in mtrs:
        attr = getattr(mtr, "user_readback")
        attr.name = mtr.name

# ## Live specfile exporting
import time
from event_model import RunRouter
# from suitcase.specfile import Serializer
from utils.CustomSpecWriter import CustomSpecWriter

# def spec_factory(name, doc):
#     if not spec_factory.enabled:
#         return [], []
#     directory = "/nsls2/data/ixs/legacy/specfiles/"

#     spec_cb = Serializer(directory, file_prefix=spec_factory.prefix, flush=True)
#     return [spec_cb], []

def my_spec_factory(name, doc):
    if name != "start" or not my_spec_factory.enabled:
        return [], []

    directory = getattr(my_spec_factory, "directory", "/nsls2/data/ixs/legacy/specfiles/")
    prefix = getattr(my_spec_factory, "prefix", "spec_test")
    filepath = os.path.join(directory, f"{prefix}.dat")

    cb = CustomSpecWriter(
        filepath=filepath,
        motor_groups=motor_groups,
        motors_per_line=8,  # or 4 if you want #O0/#O1 splitting
        include_md_keys={"uid", "detectors", "motors", "num_points", "num_intervals", "plan_pattern", "plan_pattern_args"},
        g0_items=g0_items,
        # Optional: force x field to match your scanned axis naming conventions
        # x_field_resolver=...,
        # data_field_order=...,
        flush=True,
    )
    return [cb], []


# spec_factory.enabled = True
my_spec_factory.enabled = True # Set to False to disable spec file writing without removing the callback

# Check if the 'spec_file' key exists and is not empty
if RE.md.get('spec_file'):
    config_file = RE.md['spec_file']
    directory, prefix = os.path.split(config_file)
    prefix = os.path.splitext(prefix)[0]  # remove .spec if present
    my_spec_factory.directory = directory
    my_spec_factory.prefix = prefix
    # spec_factory.directory = directory
    # spec_factory.prefix = prefix
else:
    # spec_factory.prefix = "spec_test"
    my_spec_factory.prefix = "spec_test"

# spec_router = RunRouter([spec_factory])
# RE.subscribe(spec_router)
spec_router = RunRouter([my_spec_factory])
RE.subscribe(spec_router)
