from ophyd_async.core import init_devices, print_connection_status
from ophyd_async.epics.adcore import ADWriterFactory

from pdftools.detectors import XSPDetector
from nslsii.ophyd_async.providers import NSLS2PathProvider

path_provider = NSLS2PathProvider(RE.md)

with init_devices(timeout=2.0, on_connect_hook=print_connection_status):
    lambda_flex = XSPDetector("XF:10ID1-ES{XSPD-Det:1}", ADWriterFactory.hdf(path_provider), name="lambdaflex-1")
