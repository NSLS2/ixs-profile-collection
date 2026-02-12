# sd.baseline.append(spec)
sd.baseline.append(sample_stage)
sd.baseline.append(hklps)
sd.baseline.append(hkl_params)

motor_groups = [
        ("spec", [spec.tth, spec.th, spec.chi, spec.phi]),
        ("sample_stage", [sample_stage.sx, sample_stage.sy, sample_stage.sz]),
    ]

g0_items = [
        hklps.H,                 # 1) H
        hklps.K,                 # 2) K
        hklps.L,                 # 3) L
        hklps.sc.LAMBDA,       # 4) lambda
        hkl_params.alpha,      # 5) alpha
        hkl_params.beta,       # 6) beta
        hkl_params.omega,      # 7) omega
        hkl_params.mu,         # 8) mu
        hkl_params.gam,        # 9) gam
        hkl_params.sa,         # 10) sa
        hkl_params.azimuth,    # 11) azimuth
]