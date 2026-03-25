import tkinter as tk
from tkinter import ttk, messagebox

import numpy as np
from epics import caget

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


DEXELA = "XF:10IDD-ES{Dexela:1}"
IMAGE = f"{DEXELA}image1:"

ARRAY_DATA_PV = f"{IMAGE}ArrayData"
SIZE_X_PV = f"{IMAGE}ArraySize0_RBV"
SIZE_Y_PV = f"{IMAGE}ArraySize1_RBV"


class DexelaMatplotlibGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Dexela EPICS Viewer")
        self.root.geometry("1500x950")

        self.raw_image = None
        self.display_image = None
        self.marker_screen = None       # last selected COM, still useful as a generic “current” selection
        self.direct_point = None        # stored direct beam position in CSS coordinates
        self.reflected_point = None     # stored reflected beam position in CSS coordinates
        self.measure_mode = tk.StringVar(value="direct")
        self.image_artist = None
        self.px_size = 0.0748   # mm/pixel
        self.L = 140.0          # mm

        self.build_ui()

    def build_ui(self):
        main = ttk.Frame(self.root, padding=8)
        main.pack(fill="both", expand=True)

        controls = ttk.Frame(main, width=320)
        controls.pack(side="left", fill="y", padx=(0, 12))
        controls.pack_propagate(False)

        viewer = ttk.Frame(main)
        viewer.pack(side="left", fill="both", expand=True)

        ttk.Label(controls, text="Dexela image viewer", font=("Arial", 12, "bold")).pack(anchor="w", pady=(0, 8))
        ttk.Button(controls, text="Load image from EPICS", command=self.load_epics_image).pack(fill="x", pady=2)
        ttk.Button(controls, text="Clear marker", command=self.clear_marker).pack(fill="x", pady=2)

        ttk.Separator(controls, orient="horizontal").pack(fill="x", pady=8)

        ttk.Label(controls, text="Measurement mode", font=("Arial", 10, "bold")).pack(anchor="w")
        ttk.Radiobutton(controls, text="Measure direct beam", variable=self.measure_mode, value="direct").pack(anchor="w")
        ttk.Radiobutton(controls, text="Measure reflected beam", variable=self.measure_mode, value="reflected").pack(anchor="w")

        ttk.Separator(controls, orient="horizontal").pack(fill="x", pady=8)

        # self.pv_label = tk.StringVar(value=f"Image PV: {ARRAY_DATA_PV}")
        # self.shape_label = tk.StringVar(value="Image shape: --")
        # self.coord_label = tk.StringVar(value="Selected point: --")
        self.direct_label = tk.StringVar(value="Direct beam: --")
        self.reflected_label = tk.StringVar(value="Reflected beam: --")
        self.twotheta_label = tk.StringVar(value="2Theta: --")
        self.tilt_label = tk.StringVar(value="Tilt: --")
        self.hover_label = tk.StringVar(value="Cursor: --")
        self.status_label = tk.StringVar(value="Load image from EPICS.")

        # ttk.Label(controls, textvariable=self.pv_label, wraplength=320, justify="left").pack(anchor="w", pady=2)
        # ttk.Label(controls, textvariable=self.shape_label, wraplength=320, justify="left").pack(anchor="w", pady=2)
        # ttk.Label(controls, textvariable=self.coord_label, font=("Arial", 10), wraplength=320, justify="left").pack(anchor="w", pady=(10, 2))
        ttk.Label(controls, textvariable=self.direct_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))
        ttk.Label(controls, textvariable=self.reflected_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))
        ttk.Label(controls, textvariable=self.hover_label, wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))
        ttk.Label(controls, textvariable=self.status_label, foreground="blue", wraplength=320, justify="left").pack(anchor="w", pady=(10, 2))
        ttk.Label(controls, textvariable=self.twotheta_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(6, 2))
        ttk.Label(controls, textvariable=self.tilt_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))

        # help_text = "Click on image to define beam position"
        # ttk.Label(controls, text=help_text, wraplength=320, justify="left").pack(anchor="w", pady=(10, 0))

        self.fig = Figure(figsize=(10.5, 8.5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor("black")
        self.ax.set_xlabel("X (pixels)")
        self.ax.set_ylabel("Y (pixels)")
        self.ax.set_title("Dexela image displayed in CSS coordinate convention")

        self.canvas = FigureCanvasTkAgg(self.fig, master=viewer)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.mpl_connect("button_press_event", self.on_plot_click)
        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)

    def update_calculations(self):
        if self.direct_point is None or self.reflected_point is None:
            self.twotheta_label.set("2Theta: --")
            self.tilt_label.set("Tilt: --")
            return

        x0, y0 = self.direct_point
        x1, y1 = self.reflected_point

        dx_px = x1 - x0
        dy_px = y1 - y0
        dr_px = np.hypot(dx_px, dy_px)
        dr_mm = dr_px * self.px_size

        two_theta_rad = np.arctan2(dr_mm, self.L)
        two_theta_deg = np.degrees(two_theta_rad)

        tilt_rad = np.arctan2((y1 - y0), (x0 - x1))
        tilt_deg = np.degrees(tilt_rad)

        self.twotheta_label.set(f"2Theta: {two_theta_deg:.5f} deg")
        self.tilt_label.set(f"Tilt: {tilt_deg:.5f} deg")

    def read_dexela_image(self) -> np.ndarray:
        sx = caget(SIZE_X_PV)
        sy = caget(SIZE_Y_PV)
        data = caget(ARRAY_DATA_PV, as_numpy=True)

        if sx is None or sy is None:
            raise RuntimeError("Could not read image size PVs")
        if data is None:
            raise RuntimeError("Could not read ArrayData PV")

        sx = int(sx)
        sy = int(sy)
        flat = np.asarray(data).ravel()
        expected = sx * sy

        if flat.size < expected:
            raise RuntimeError(f"ArrayData has {flat.size} values, expected {expected}")

        # Keep this identical to the earlier working reader.
        img = flat[:expected].reshape((sy, sx)).astype(np.float32)
        return img

    def load_epics_image(self):
        try:
            self.raw_image = self.read_dexela_image()
        except Exception as exc:
            messagebox.showerror("EPICS read error", str(exc))
            return

        self.display_image = np.flipud(self.raw_image)
        self.status_label.set("Click on image to define beam position")
        self.redraw_image()

    def redraw_image(self):
        if self.display_image is None:
            return

        ny, nx = self.display_image.shape
        vmin = np.percentile(self.raw_image, 1)
        vmax = np.percentile(self.raw_image, 99.5)

        self.ax.clear()
        self.ax.set_facecolor("black")

        self.image_artist = self.ax.imshow(
            self.display_image,
            cmap="gray",
            origin="lower",
            extent=[0, nx - 1, 0, ny - 1],
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )

        self.ax.set_xlabel("X (pixels)")
        self.ax.set_ylabel("Y (pixels)")
        self.ax.set_title("Dexela Image")
        self.ax.set_xlim(0, nx - 1)
        self.ax.set_ylim(0, ny - 1)

        if self.direct_point is not None:
            self.draw_marker(self.direct_point, color="lime")

        if self.reflected_point is not None:
            self.draw_marker(self.reflected_point, color="red")

        self.canvas.draw()

    def is_valid_click_region(self, point_css, half_size=20):
        if self.display_image is None:
            return False

        x_css, y_css = point_css
        ny, nx = self.display_image.shape

        xc = int(round(x_css))
        yc = int(round(y_css))

        x0 = max(0, xc - half_size)
        x1 = min(nx, xc + half_size + 1)
        y0 = max(0, yc - half_size)
        y1 = min(ny, yc + half_size + 1)

        roi = self.display_image[y0:y1, x0:x1]
        if roi.size == 0:
            return False

        max_intensity = float(np.max(roi))

        border_values = []
        if roi.shape[0] >= 1:
            border_values.append(roi[0, :])
            if roi.shape[0] > 1:
                border_values.append(roi[-1, :])
        if roi.shape[1] >= 1:
            if roi.shape[0] > 2:
                border_values.append(roi[1:-1, 0])
                if roi.shape[1] > 1:
                    border_values.append(roi[1:-1, -1])
            elif roi.shape[1] > 1:
                border_values.append(roi[:, 0])
                border_values.append(roi[:, -1])

        if not border_values:
            return False

        border = np.concatenate([np.ravel(v) for v in border_values if np.size(v) > 0])
        if border.size == 0:
            return False

        background = float(np.median(border))
        return max_intensity >= 1.0 * background

    def on_plot_click(self, event):
        if self.display_image is None:
            return
        if event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        x_css = float(event.xdata)
        y_css = float(event.ydata)

        ny, nx = self.display_image.shape
        x_css = min(max(x_css, 0.0), nx - 1.0)
        y_css = min(max(y_css, 0.0), ny - 1.0)
        point = (x_css, y_css)

        if not self.is_valid_click_region(point):
            self.marker_screen = None
            # self.coord_label.set("Selected point: --")
            self.status_label.set("Selection canceled: no peak is found.")
            self.redraw_image()
            return

        com_point = self.compute_com(point)
        if com_point is None:
            self.marker_screen = None
            # self.coord_label.set("Selected point: --")
            self.status_label.set("Selection canceled: could not compute center of mass.")
            self.redraw_image()
            return

        self.marker_screen = com_point

        mode = self.measure_mode.get()
        if mode == "direct":
            self.direct_point = com_point
            self.direct_label.set(f"Direct beam: x={com_point[0]:.2f}, y={com_point[1]:.2f}")
        elif mode == "reflected":
            self.reflected_point = com_point
            self.reflected_label.set(f"Reflected beam: x={com_point[0]:.2f}, y={com_point[1]:.2f}")

        self.status_label.set("Click on image to define beam position")
        self.update_calculations()
        self.redraw_image()

    def on_mouse_move(self, event):
        if self.display_image is None:
            self.hover_label.set("Cursor: --")
            return

        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            self.hover_label.set("Cursor: --")
            return

        ny, nx = self.display_image.shape

        x_css = min(max(float(event.xdata), 0.0), nx - 1.0)
        y_css = min(max(float(event.ydata), 0.0), ny - 1.0)

        ix = int(round(x_css))
        iy = int(round(y_css))

        ix = min(max(ix, 0), nx - 1)
        iy = min(max(iy, 0), ny - 1)

        intensity = float(self.display_image[iy, ix])

        self.hover_label.set(
            f"Cursor: x={x_css:.2f}, y={y_css:.2f}, I={intensity:.1f}"
        )

    def draw_marker(self, point_css, color="red"):
        x_css, y_css = point_css

        self.ax.scatter(
            [x_css], [y_css],
            s=700,
            facecolors="none",
            edgecolors=color,
            linewidths=1.5,
            zorder=10,
        )

    def compute_com(self, point_css, half_size=20):
        if self.display_image is None:
            return None

        x_css, y_css = point_css
        ny, nx = self.display_image.shape

        xc = int(round(x_css))
        yc = int(round(y_css))

        x0 = max(0, xc - half_size)
        x1 = min(nx, xc + half_size + 1)
        y0 = max(0, yc - half_size)
        y1 = min(ny, yc + half_size + 1)

        roi = self.display_image[y0:y1, x0:x1].astype(np.float64)
        if roi.size == 0:
            return None

        border_values = []
        if roi.shape[0] >= 1:
            border_values.append(roi[0, :])
            if roi.shape[0] > 1:
                border_values.append(roi[-1, :])
        if roi.shape[1] >= 1:
            if roi.shape[0] > 2:
                border_values.append(roi[1:-1, 0])
                if roi.shape[1] > 1:
                    border_values.append(roi[1:-1, -1])
            elif roi.shape[1] > 1:
                border_values.append(roi[:, 0])
                border_values.append(roi[:, -1])

        if not border_values:
            return None

        border = np.concatenate([np.ravel(v) for v in border_values if np.size(v) > 0])
        if border.size == 0:
            return None

        background = float(np.median(border))
        weights = roi - background
        weights[weights < 0] = 0.0
        total = weights.sum()
        if total <= 0:
            return None

        yy, xx = np.indices(roi.shape)
        x_com = x0 + float((xx * weights).sum() / total)
        y_com = y0 + float((yy * weights).sum() / total)
        return (x_com, y_com)

    def clear_marker(self):
        self.marker_screen = None
        self.direct_point = None
        self.reflected_point = None
        self.direct_label.set("Direct beam: --")
        self.reflected_label.set("Reflected beam: --")
        self.status_label.set("Click on image to define beam position")
        self.update_calculations()
        self.redraw_image()


def main():
    root = tk.Tk()
    DexelaMatplotlibGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
