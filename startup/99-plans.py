def check_zero():
    ssxop = 0
    E0
    print('scanning zero')
    yield from bps.mov(tth, 0, phi, 0, ssx, -1)
    sample_pos = yield from bps.read(sample_stage)
    print(sample_pos)
    dets = [det1]
    yield from bps.mov(whl, 5)
    yield from bps.mov(hrmE, E0 - 20)
    for d in dets:
        # set the exposure times
        pass

    local_peaks = bluesky.callbacks.fitting.PeakStats(hrmE.energy.name, 
                                                      dets[0].hints['fields'][0])
    plan = bpp.subs_wrapper(rel_scan(dets, hrmE, 0, 40, 200, md={'reason': 'alignment'}), 
                            [local_peaks])
                                                      
    yield from plan
    cen = local_peaks.cen
    target = 0.2 * (cen // .2)
    yield from bps.mov(hrmE, target)

