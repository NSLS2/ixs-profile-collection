import logging
from bluesky.standard_config import *  # gs, etc.
import matplotlib.pyplot as plt
plt.ion()
from bluesky import qt_kicker
qt_kicker.install_qt_kicker()
from databroker import DataBroker as db, get_events, get_images, get_table

from epics import caput, caget

# connect olog

gs.RE.logbook = olog_wrapper(olog_client, ['Data Acquisition'])
RE=gs.RE
from bluesky.scientific_callbacks import plot_peak_stats
# from chxtools.xfuncs import *
# from chxtools.plot import plot1
from bluesky.scans import Count,AbsScan, DeltaScan
