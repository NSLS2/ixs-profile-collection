"""X-Ray Eye: beam spot viewer / FWHM / focus-scan analysis for the XrayEye camera.

Modelled on DexelaCalc.py. The pure analysis functions have no Tk or EPICS
imports at module level so they can be unit tested.
"""
import os

import numpy as np

CAMERA = "XF:10IDD-BI{XrayEye-Cam:3}"
IMAGE = f"{CAMERA}image1:"
ARRAY_DATA_PV = f"{IMAGE}ArrayData"
SIZE_X_PV = f"{IMAGE}ArraySize0_RBV"
SIZE_Y_PV = f"{IMAGE}ArraySize1_RBV"

EXPECTED_SHAPE = (484, 644)  # (rows, cols), Mako G-030B
DEFAULT_PX_UM = 0.38
DEFAULT_FOLDER = "/IXS/CAMERA/KB"
DEFAULT_PREFIX = "XRY"
DEFAULT_DIGITS = 3
DEFAULT_EXT = "tif"


# ----------------------------------------------------------------------
# Pure analysis functions
# ----------------------------------------------------------------------
def fwhm(x, y, level=0.5):
    """Full width at `level` of the maximum, with linear sub-pixel crossings.

    Returns nan if the profile does not fall below the level on both sides.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3 or y.size != x.size:
        return float("nan")
    ipk = int(np.argmax(y))
    thr = level * y[ipk]
    if y[ipk] <= 0:
        return float("nan")

    # left crossing
    left = None
    for i in range(ipk, 0, -1):
        if y[i - 1] < thr <= y[i]:
            left = x[i - 1] + (thr - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
            break
    right = None
    for i in range(ipk, y.size - 1):
        if y[i + 1] < thr <= y[i]:
            right = x[i] + (y[i] - thr) * (x[i + 1] - x[i]) / (y[i] - y[i + 1])
            break
    if left is None or right is None:
        return float("nan")
    return float(right - left)


def _roi_bounds(shape, point, half_size):
    ny, nx = shape
    xc, yc = int(round(point[0])), int(round(point[1]))
    return (max(0, yc - half_size), min(ny, yc + half_size + 1),
            max(0, xc - half_size), min(nx, xc + half_size + 1))


def _border_median(roi):
    if roi.size == 0:
        return None
    parts = [roi[0, :], roi[-1, :]]
    if roi.shape[0] > 2:
        parts += [roi[1:-1, 0], roi[1:-1, -1]]
    return float(np.median(np.concatenate([np.ravel(p) for p in parts])))


def has_peak(img, point, half_size=20, ratio=2.0):
    """True if the ROI maximum is at least `ratio` times the border median."""
    y0, y1, x0, x1 = _roi_bounds(img.shape, point, half_size)
    roi = img[y0:y1, x0:x1]
    bg = _border_median(roi)
    if bg is None:
        return False
    return float(roi.max()) >= ratio * bg


def centroid(img, point, half_size=20):
    """Background-subtracted centre of mass in a ROI. Returns (x, y) or None."""
    y0, y1, x0, x1 = _roi_bounds(img.shape, point, half_size)
    roi = np.asarray(img[y0:y1, x0:x1], dtype=np.float64)
    bg = _border_median(roi)
    if bg is None:
        return None
    w = np.clip(roi - bg, 0.0, None)
    total = w.sum()
    if total <= 0:
        return None
    yy, xx = np.indices(roi.shape)
    return (x0 + float((xx * w).sum() / total), y0 + float((yy * w).sum() / total))


def find_max(img):
    """(x, y) of the brightest pixel."""
    r, c = np.unravel_index(int(np.argmax(img)), img.shape)
    return (float(c), float(r))


def profiles(img, point, px_um=1.0):
    """Horizontal and vertical cuts through `point`.

    Returns (x_um, ix, y_um, iy).
    """
    ny, nx = img.shape
    c = min(max(int(round(point[0])), 0), nx - 1)
    r = min(max(int(round(point[1])), 0), ny - 1)
    return (px_um * np.arange(nx), np.asarray(img[r, :], dtype=float),
            px_um * np.arange(ny), np.asarray(img[:, c], dtype=float))


def build_filename(folder, prefix, n, digits=DEFAULT_DIGITS, ext=DEFAULT_EXT):
    return os.path.join(folder, f"{prefix}_{n:0{digits}d}.{ext.lstrip('.')}")


def load_tiff(path):
    from PIL import Image
    with Image.open(path) as im:
        return np.asarray(im, dtype=np.float32)


def measure_image(img, px_um, mode="click", point=None, half_size=20):
    """Locate the beam and return dict(point, fwhm_x, fwhm_y, imax) or None."""
    if mode == "auto" or point is None:
        pt = find_max(img)
    else:
        pt = centroid(img, point, half_size)
        if pt is None:
            return None
    xx, ix, yy, iy = profiles(img, pt, px_um)
    # background-subtract each profile by its minimum, as the old scripts did
    return {
        "point": pt,
        "fwhm_x": fwhm(xx, ix - ix.min()),
        "fwhm_y": fwhm(yy, iy - iy.min()),
        "imax": float(img.max()),
    }


def focus_scan(folder, prefix, n1, n2, direction, px_um,
               digits=DEFAULT_DIGITS, ext=DEFAULT_EXT):
    """Width and intensity vs. file number using the brightest pixel per image.

    direction: "x" or "y". Returns dict with arrays n, width, imax, the best
    file numbers, and a list of (n, reason) for skipped files.
    """
    ns, widths, imaxs, skipped = [], [], [], []
    for n in range(n1, n2 + 1):
        path = build_filename(folder, prefix, n, digits, ext)
        try:
            img = load_tiff(path)
        except Exception as exc:  # missing/unreadable file
            skipped.append((n, str(exc)))
            continue
        m = measure_image(img, px_um, mode="auto")
        w = m["fwhm_x"] if direction.lower() == "x" else m["fwhm_y"]
        ns.append(n)
        widths.append(w)
        imaxs.append(m["imax"])
    ns = np.array(ns)
    widths = np.array(widths, dtype=float)
    imaxs = np.array(imaxs, dtype=float)
    n_wmin = n_imax = None
    if ns.size:
        if np.any(np.isfinite(widths)):
            n_wmin = int(ns[np.nanargmin(widths)])
        n_imax = int(ns[np.argmax(imaxs)])
    return {"n": ns, "width": widths, "imax": imaxs,
            "n_wmin": n_wmin, "n_imax": n_imax, "skipped": skipped}


# ----------------------------------------------------------------------
# EPICS
# ----------------------------------------------------------------------
def read_camera_image():
    from epics import caget
    sx = caget(SIZE_X_PV)
    sy = caget(SIZE_Y_PV)
    data = caget(ARRAY_DATA_PV, as_numpy=True)
    if sx is None or sy is None:
        raise RuntimeError("Could not read image size PVs")
    if data is None:
        raise RuntimeError("Could not read ArrayData PV")
    sx, sy = int(sx), int(sy)
    flat = np.asarray(data).ravel()
    if flat.size < sx * sy:
        raise RuntimeError(f"ArrayData has {flat.size} values, expected {sx * sy}")
    return flat[:sx * sy].reshape((sy, sx)).astype(np.float32)


# ----------------------------------------------------------------------
# GUI
# ----------------------------------------------------------------------
class XRayEyeGUI:
    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        self.tk, self.ttk = tk, ttk
        self.root = root
        root.title("X-Ray Eye")
        root.geometry("1500x950")

        self.img = None
        self.point = None
        self.pick_mode = tk.StringVar(value="click")
        self.px_um = tk.StringVar(value=str(DEFAULT_PX_UM))
        self.half = tk.StringVar(value="20")
        self.folder = tk.StringVar(value=DEFAULT_FOLDER)
        self.prefix = tk.StringVar(value=DEFAULT_PREFIX)
        self.digits = tk.StringVar(value=str(DEFAULT_DIGITS))
        self.ext = tk.StringVar(value=DEFAULT_EXT)
        self.n1 = tk.StringVar(value="1")
        self.n2 = tk.StringVar(value="10")
        self.direction = tk.StringVar(value="x")
        self.center_var = tk.StringVar(value="Center: --")
        self.fwx_var = tk.StringVar(value="FWHM X: --")
        self.fwy_var = tk.StringVar(value="FWHM Y: --")
        self.imax_var = tk.StringVar(value="Imax: --")
        self.status = tk.StringVar(value="Load an image")

        main = ttk.Frame(root, padding=8)
        main.pack(fill="both", expand=True)
        ctl = ttk.Frame(main, width=330)
        ctl.pack(side="left", fill="y", padx=(0, 12))
        ctl.pack_propagate(False)
        view = ttk.Frame(main)
        view.pack(side="left", fill="both", expand=True)

        bold = ("Arial", 11, "bold")
        ttk.Label(ctl, text="X-Ray Eye", font=("Arial", 12, "bold")).pack(anchor="w", pady=(0, 6))
        ttk.Button(ctl, text="Load from camera", command=self.load_camera).pack(fill="x", pady=2)
        ttk.Button(ctl, text="Load file...", command=self.load_file).pack(fill="x", pady=2)
        ttk.Separator(ctl).pack(fill="x", pady=6)

        ttk.Label(ctl, text="Pick mode", font=bold).pack(anchor="w")
        ttk.Radiobutton(ctl, text="Click + ROI centroid", variable=self.pick_mode,
                        value="click").pack(anchor="w")
        ttk.Radiobutton(ctl, text="Auto (brightest pixel)", variable=self.pick_mode,
                        value="auto", command=self.analyze).pack(anchor="w")
        self._field(ctl, "Pixel size (um)", self.px_um)
        self._field(ctl, "ROI half size (px)", self.half)
        ttk.Separator(ctl).pack(fill="x", pady=6)

        for v in (self.center_var, self.fwx_var, self.fwy_var, self.imax_var):
            ttk.Label(ctl, textvariable=v, font=bold).pack(anchor="w", pady=1)
        ttk.Separator(ctl).pack(fill="x", pady=6)

        ttk.Label(ctl, text="Focus scan (file range)", font=bold).pack(anchor="w")
        self._field(ctl, "Folder", self.folder)
        self._field(ctl, "Prefix", self.prefix)
        self._field(ctl, "Digits", self.digits)
        self._field(ctl, "Extension", self.ext)
        self._field(ctl, "First #", self.n1)
        self._field(ctl, "Last #", self.n2)
        row = ttk.Frame(ctl)
        row.pack(anchor="w")
        ttk.Radiobutton(row, text="X width", variable=self.direction, value="x").pack(side="left")
        ttk.Radiobutton(row, text="Y width", variable=self.direction, value="y").pack(side="left")
        ttk.Button(ctl, text="Run scan analysis", command=self.run_scan).pack(fill="x", pady=4)
        ttk.Label(ctl, textvariable=self.status, foreground="blue",
                  wraplength=320, justify="left").pack(anchor="w", pady=6)

        self.fig = Figure(figsize=(10.5, 8.5), dpi=100)
        gs = self.fig.add_gridspec(2, 2, height_ratios=[3, 1])
        self.ax = self.fig.add_subplot(gs[0, :])
        self.axx = self.fig.add_subplot(gs[1, 0])
        self.axy = self.fig.add_subplot(gs[1, 1])
        self.canvas = FigureCanvasTkAgg(self.fig, master=view)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.mpl_connect("button_press_event", self.on_click)

    def _field(self, parent, label, var):
        row = self.ttk.Frame(parent)
        row.pack(fill="x", pady=1)
        self.ttk.Label(row, text=label, width=18).pack(side="left")
        self.ttk.Entry(row, textvariable=var, width=18).pack(side="left", fill="x", expand=True)

    def _err(self, title, msg):
        from tkinter import messagebox
        messagebox.showerror(title, msg)

    def _float(self, var, name):
        try:
            return float(var.get())
        except ValueError:
            raise ValueError(f"{name} must be a number")

    # ---- loading
    def _set_image(self, img, source):
        if img.shape != EXPECTED_SHAPE:
            self.status.set(f"Warning: shape {img.shape}, expected {EXPECTED_SHAPE}. {source}")
        else:
            self.status.set(f"Loaded {source}")
        self.img = img
        self.point = None
        self.analyze()

    def load_camera(self):
        try:
            img = read_camera_image()
        except Exception as exc:
            self._err("EPICS read error", str(exc))
            return
        self._set_image(img, "camera image")

    def load_file(self):
        from tkinter import filedialog
        path = filedialog.askopenfilename(
            initialdir=self.folder.get() if os.path.isdir(self.folder.get()) else None,
            filetypes=[("TIFF", "*.tif *.tiff"), ("All", "*.*")])
        if not path:
            return
        try:
            img = load_tiff(path)
        except Exception as exc:
            self._err("File error", str(exc))
            return
        self._set_image(img, os.path.basename(path))

    # ---- analysis
    def on_click(self, event):
        if self.img is None or event.inaxes != self.ax or event.xdata is None:
            return
        if self.pick_mode.get() != "click":
            return
        try:
            half = int(self._float(self.half, "ROI half size"))
        except ValueError as exc:
            self._err("Input error", str(exc))
            return
        pt = (event.xdata, event.ydata)
        if not has_peak(self.img, pt, half):
            self.status.set("Selection canceled: no peak found.")
            return
        self.point = pt
        self.analyze()

    def analyze(self):
        if self.img is None:
            return
        try:
            px = self._float(self.px_um, "Pixel size")
            half = int(self._float(self.half, "ROI half size"))
        except ValueError as exc:
            self._err("Input error", str(exc))
            return
        mode = self.pick_mode.get()
        if mode == "click" and self.point is None:
            self.redraw(None)
            self.status.set("Click on the beam to measure it")
            return
        res = measure_image(self.img, px, mode, self.point, half)
        if res is None:
            self.status.set("Could not compute centroid.")
            self.redraw(None)
            return
        self.center_var.set(f"Center: x={res['point'][0]:.1f}, y={res['point'][1]:.1f}")
        self.fwx_var.set(f"FWHM X: {res['fwhm_x']:.2f} um")
        self.fwy_var.set(f"FWHM Y: {res['fwhm_y']:.2f} um")
        self.imax_var.set(f"Imax: {res['imax']:.0f}")
        self.redraw(res, px)

    def redraw(self, res, px=1.0):
        self.ax.clear()
        self.axx.clear()
        self.axy.clear()
        ny, nx = self.img.shape
        self.ax.imshow(self.img, cmap="gray", origin="upper",
                       vmin=np.percentile(self.img, 1), vmax=np.percentile(self.img, 99.5),
                       interpolation="nearest")
        self.ax.set_title("XrayEye image (as stored, origin top-left)")
        if res is not None:
            x0, y0 = res["point"]
            self.ax.scatter([x0], [y0], s=500, facecolors="none", edgecolors="lime")
            xx, ix, yy, iy = profiles(self.img, res["point"], px)
            self.axx.plot(xx, ix)
            self.axx.set_xlabel("X (um)")
            self.axx.set_title(f"FWHM = {res['fwhm_x']:.2f} um")
            self.axy.plot(yy, iy)
            self.axy.set_xlabel("Y (um)")
            self.axy.set_title(f"FWHM = {res['fwhm_y']:.2f} um")
        self.fig.tight_layout()
        self.canvas.draw()

    # ---- focus scan
    def run_scan(self):
        try:
            px = self._float(self.px_um, "Pixel size")
            n1 = int(self._float(self.n1, "First #"))
            n2 = int(self._float(self.n2, "Last #"))
            digits = int(self._float(self.digits, "Digits"))
        except ValueError as exc:
            self._err("Input error", str(exc))
            return
        if n2 < n1:
            self._err("Input error", "Last # must be >= First #")
            return
        res = focus_scan(self.folder.get(), self.prefix.get(), n1, n2,
                         self.direction.get(), px, digits, self.ext.get())
        if res["n"].size == 0:
            self._err("Focus scan", "No readable files in range.")
            return
        d = self.direction.get().upper()
        msg = (f"Min {d} width at file {res['n_wmin']}; "
               f"max intensity at file {res['n_imax']}")
        if res["skipped"]:
            msg += f"; skipped {len(res['skipped'])} file(s)"
        self.status.set(msg)

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        win = self.tk.Toplevel(self.root)
        win.title("Focus scan")
        fig = Figure(figsize=(8, 7), dpi=100)
        a1 = fig.add_subplot(211)
        a2 = fig.add_subplot(212)
        a1.plot(res["n"], res["width"], "o-")
        a1.set_ylabel(f"{d} width (um)")
        a2.plot(res["n"], res["imax"], "o-")
        a2.set_ylabel("Intensity")
        a2.set_xlabel("File #")
        fig.tight_layout()
        c = FigureCanvasTkAgg(fig, master=win)
        c.get_tk_widget().pack(fill="both", expand=True)
        c.draw()


def main():
    import tkinter as tk
    root = tk.Tk()
    XRayEyeGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
