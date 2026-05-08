import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from epics import caget

DEXELA = "XF:10IDD-ES{Dexela:1}"
IMAGE = f"{DEXELA}image1:"

ARRAY_DATA_PV = f"{IMAGE}ArrayData"
SIZE_X_PV = f"{IMAGE}ArraySize0_RBV"
SIZE_Y_PV = f"{IMAGE}ArraySize1_RBV"


def read_dexela_image():
    sx = int(caget(SIZE_X_PV))
    sy = int(caget(SIZE_Y_PV))
    data = caget(ARRAY_DATA_PV, as_numpy=True)

    if data is None:
        raise RuntimeError("Could not read ArrayData PV")

    flat = np.asarray(data).ravel()
    expected = sx * sy
    if flat.size < expected:
        raise RuntimeError(f"ArrayData has {flat.size} values, expected {expected}")

    img = flat[:expected].reshape((sy, sx)).astype(np.float32)
    return img


def raw_to_css(x_raw, y_raw, ny):
    x_css = x_raw
    y_css = (ny - 1) - y_raw
    return x_css, y_css


def css_to_raw(x_css, y_css, ny):
    x_raw = x_css
    y_raw = (ny - 1) - y_css
    return x_raw, y_raw


def display_dexela_image():
    img = read_dexela_image()
    ny, nx = img.shape

    print("shape:", img.shape)
    print("dtype:", img.dtype)
    print("min:", img.min())
    print("max:", img.max())

    vmin = np.percentile(img, 1)
    vmax = np.percentile(img, 99.5)

    # flip only for display so the shown image matches CSS-style coordinates
    img_css = np.flipud(img)

    plt.figure(figsize=(10, 8))
    plt.imshow(
        img_css,
        cmap="gray",
        origin="lower",
        extent=[0, nx - 1, 0, ny - 1],
        vmin=vmin,
        vmax=vmax,
    )
    plt.xlabel("X (CSS pixels)")
    plt.ylabel("Y (CSS pixels)")
    plt.title("Dexela image displayed in CSS coordinate convention")
    plt.colorbar()
    plt.tight_layout()
    plt.show()

    return img

if __name__ == "__main__":
    img = display_dexela_image()