"""Tests for the pure analysis functions in startup/utils/XRayEye.py."""
import importlib.util
import pathlib

import numpy as np
import pytest
from PIL import Image

SRC = pathlib.Path(__file__).resolve().parent.parent / "startup" / "utils" / "XRayEye.py"


@pytest.fixture(scope="module")
def xe():
    spec = importlib.util.spec_from_file_location("XRayEye", SRC)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def gauss_img(sx=4.0, sy=4.0, cx=300, cy=200, amp=1000.0, bg=10.0, shape=(484, 644)):
    yy, xx = np.indices(shape)
    return (bg + amp * np.exp(-((xx - cx) ** 2 / (2 * sx**2) + (yy - cy) ** 2 / (2 * sy**2)))
            ).astype(np.float32)


def test_fwhm_gaussian(xe):
    x = np.linspace(-20, 20, 801)
    s = 3.0
    assert xe.fwhm(x, np.exp(-x**2 / (2 * s**2))) == pytest.approx(2.3548 * s, rel=1e-2)


def test_fwhm_flat_is_nan(xe):
    x = np.arange(50.0)
    assert np.isnan(xe.fwhm(x, np.ones(50)))


def test_find_max(xe):
    assert xe.find_max(gauss_img()) == (300.0, 200.0)


def test_centroid_with_background(xe):
    c = xe.centroid(gauss_img(cx=300.3, cy=199.6, bg=50.0), (295, 205))
    assert c[0] == pytest.approx(300.3, abs=0.1)
    assert c[1] == pytest.approx(199.6, abs=0.1)


def test_has_peak_rejects_noise(xe):
    rng = np.random.default_rng(0)
    noise = rng.normal(100, 2, (484, 644)).astype(np.float32)
    assert not xe.has_peak(noise, (300, 200))
    assert xe.has_peak(gauss_img(), (300, 200))


def test_build_filename(xe):
    assert xe.build_filename("/d", "XRY", 7) == "/d/XRY_007.tif"


def test_load_tiff_roundtrip(xe, tmp_path):
    arr = (np.arange(12, dtype=np.uint16).reshape(3, 4) * 100)
    p = tmp_path / "a.tif"
    Image.fromarray(arr).save(p)
    out = xe.load_tiff(str(p))
    assert out.shape == (3, 4)
    assert np.array_equal(out, arr.astype(np.float32))


def test_measure_image_fwhm(xe):
    m = xe.measure_image(gauss_img(sx=5.0, sy=3.0), px_um=0.38, mode="auto")
    assert m["fwhm_x"] == pytest.approx(2.3548 * 5.0 * 0.38, rel=0.03)
    assert m["fwhm_y"] == pytest.approx(2.3548 * 3.0 * 0.38, rel=0.03)


def test_focus_scan(xe, tmp_path):
    widths = {1: 8.0, 2: 5.0, 3: 2.0, 4: 6.0}  # file 3 sharpest; file 5 missing
    for n, s in widths.items():
        amp = 1000.0 + (500.0 if n == 2 else 0.0)
        Image.fromarray(gauss_img(sx=s, amp=amp).astype(np.uint16)).save(
            tmp_path / f"XRY_{n:03d}.tif")
    res = xe.focus_scan(str(tmp_path), "XRY", 1, 5, "x", 0.38)
    assert list(res["n"]) == [1, 2, 3, 4]
    assert res["n_wmin"] == 3
    assert res["n_imax"] == 2
    assert [n for n, _ in res["skipped"]] == [5]
