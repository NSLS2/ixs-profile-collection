#gs.DETS = [s1_1, s1_2, s1_3, s1_4, s2_1, s2_2, s2_3, s2_4]
all_objs = globals()
counters = [counter for counter in all_objs.values() if isinstance(counter, EpicsSignal)]

gs.DETS = counters


import logging

from ophyd.session import get_session_manager

sessionmgr = get_session_manager()
sessionmgr['olog_client'] = olog_client
print('These positioners are disconnected:')
print([k for k, v in sessionmgr.get_positioners().items() if not v.connected])

# metadata set at startup
gs.RE.md['owner'] = 'xf10id'
gs.RE.md['group'] = 'ixs'
gs.RE.md['beamline_id'] = 'IXS'
#gs.RE.md['custom'] = {}



def print_scanid(name, doc):
    if name == 'start':
        print('Scan ID:', doc['scan_id'])
        print('Unique ID:', doc['uid'])

def print_md(name, doc):
    if name == 'start':
        print('Metadata:\n', repr(doc))

gs.RE.subscribe('start', print_scanid)
from ophyd.commands import wh_pos, log_pos, mov, movr



