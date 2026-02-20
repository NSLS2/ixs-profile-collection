# sd.baseline.append(spec)
sd.baseline.append(sample_stage)
sd.baseline.append(hklps)
sd.baseline.append(hkl_params)

motor_groups = [
        ("spec", [spec.tth, spec.th, spec.chi, spec.phi]),
        ("sample_stage", [sample_stage.sx, sample_stage.sy, sample_stage.sz]),
        ("sample_prime", [sam_prime.xp, sam_prime.zp]),
        ("undulator", [ivu22]),
        ("energy", [blE.energy, dcmE.energy, hrmE.energy]),
    ]

g0_items = [
        hklps.sc.LAMBDA,       # 1) lambda
        hkl_params.alpha,      # 2) alpha
        hkl_params.beta,       # 3) beta
        hkl_params.omega,      # 4) omega
        hkl_params.mu,         # 5) mu
        hkl_params.gam,        # 6) gam
        hkl_params.sa,         # 7) sa
        hkl_params.azimuth,    # 8) azimuth
]

g1_items = [
    hklps.sc.g_aa, hklps.sc.g_bb, hklps.sc.g_cc,
    hklps.sc.g_al, hklps.sc.g_be, hklps.sc.g_ga,

    hklps.sc.g_aa_s, hklps.sc.g_bb_s, hklps.sc.g_cc_s,
    hklps.sc.g_al_s, hklps.sc.g_be_s, hklps.sc.g_ga_s,

    hklps.sc.g_h0, hklps.sc.g_k0, hklps.sc.g_l0,
    hklps.sc.g_u00, hklps.sc.g_u01, hklps.sc.g_u02,
    hklps.sc.g_u03, hklps.sc.g_u04, hklps.sc.g_u05,

    hklps.sc.g_h1, hklps.sc.g_k1, hklps.sc.g_l1,
    hklps.sc.g_u10, hklps.sc.g_u11, hklps.sc.g_u12,
    hklps.sc.g_u13, hklps.sc.g_u14, hklps.sc.g_u15,
]

q_items = [hklps.H, hklps.K, hklps.L]