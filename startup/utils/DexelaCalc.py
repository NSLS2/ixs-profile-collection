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
        self.marker_screen = None
        self.show_markers = True

        # basic mode points
        self.direct_point = None
        self.reflected_point = None

        # extension mode points / derived values
        self.extension_reflected_reference = None   # used only for shift calibration
        self.extension_reflected_current = None     # used for later extension measurements
        self.extension_direct_inferred = None
        self.extension_shift = None   # (dx, dy) in CSS pixels

        # four-state operation mode
        self.measure_mode = tk.StringVar(value="direct_basic")

        self.image_artist = None

        self.px_size = 0.0748   # mm/pixel
        self.L = 140.0          # mm

        self.build_ui()

    def build_ui(self):
        main = ttk.Frame(self.root, padding=8)
        main.pack(fill="both", expand=True)

        style = ttk.Style()

        style.configure("Big.TButton", font=("Arial", 12))
        style.configure("Mode.TLabel", font=("Arial", 11, "bold"))
        style.configure("Mode.TRadiobutton", font=("Arial", 11))

        controls = ttk.Frame(main, width=320)
        controls.pack(side="left", fill="y", padx=(0, 12))
        controls.pack_propagate(False)

        viewer = ttk.Frame(main)
        viewer.pack(side="left", fill="both", expand=True)

        ttk.Label(controls, text="Dexela image viewer", font=("Arial", 12, "bold")).pack(anchor="w", pady=(0, 8))
        ttk.Button(controls, text="Load image from Dexela", command=self.load_epics_image, style="Big.TButton").pack(fill="x", pady=2)
        ttk.Button(controls, text="Clear settings", command=self.clear_marker, style="Big.TButton").pack(fill="x", pady=2)

        ttk.Separator(controls, orient="horizontal").pack(fill="x", pady=8)

        ttk.Label(controls, text="Measurement mode", style="Mode.TLabel").pack(anchor="w")
        ttk.Radiobutton(controls, text="Basic: set direct beam", variable=self.measure_mode, value="direct_basic", style="Mode.TRadiobutton", command=self.update_active_geometry_label).pack(anchor="w")
        ttk.Radiobutton(controls, text="Basic: measure reflected", variable=self.measure_mode, value="reflected_basic", style="Mode.TRadiobutton", command=self.update_active_geometry_label).pack(anchor="w")
        ttk.Radiobutton(controls, text="Extension: calibrate shift", variable=self.measure_mode, value="extension_calibrate", style="Mode.TRadiobutton", command=self.update_active_geometry_label).pack(anchor="w")
        ttk.Radiobutton(controls, text="Extension: measure reflected", variable=self.measure_mode, value="extension_measure", style="Mode.TRadiobutton", command=self.update_active_geometry_label).pack(anchor="w")

        # ttk.Separator(controls, orient="horizontal").pack(fill="x", pady=8)

        self.active_geometry_label = tk.StringVar(value="Active geometry: BASIC")
        self.direct_label = tk.StringVar(value="Direct beam: --")
        self.reflected_label = tk.StringVar(value="Basic reflected beam: --")
        self.extension_ref_label = tk.StringVar(value="Extension reflected ref: --")
        self.extension_current_label = tk.StringVar(value="Extension reflected current: --")
        self.extension_shift_label = tk.StringVar(value="Extension shift: --")
        self.extension_direct_label = tk.StringVar(value="Inferred direct beam: --")
        self.twotheta_label = tk.StringVar(value="2Theta: --")
        self.tilt_label = tk.StringVar(value="Tilt: --")
        self.hover_label = tk.StringVar(value="Cursor: --")
        self.status_label = tk.StringVar(value="Click on image to define beam position")

        ttk.Label(controls, textvariable=self.direct_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(8, 2))
        ttk.Label(controls, textvariable=self.reflected_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))

        ttk.Separator(controls, orient="horizontal").pack(fill="x", pady=6)

        ttk.Label(controls, textvariable=self.extension_ref_label, font=("Arial", 10), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))
        ttk.Label(controls, textvariable=self.extension_current_label, font=("Arial", 10), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))
        ttk.Label(controls, textvariable=self.extension_shift_label, font=("Arial", 10), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))
        ttk.Label(controls, textvariable=self.extension_direct_label, font=("Arial", 10), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))

        ttk.Separator(controls, orient="horizontal").pack(fill="x", pady=6)

        ttk.Label(controls, textvariable=self.active_geometry_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(8, 4))
        ttk.Label(controls, textvariable=self.twotheta_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(6, 2))
        ttk.Label(controls, textvariable=self.tilt_label, font=("Arial", 11, "bold"), wraplength=320, justify="left").pack(anchor="w", pady=(2, 2))
        ttk.Label(controls, textvariable=self.hover_label, font=("Arial", 11), wraplength=320, justify="left").pack(anchor="w", pady=(8, 2))
        ttk.Label(controls, textvariable=self.status_label, font=("Arial", 11), foreground="blue", wraplength=320, justify="left").pack(anchor="w", pady=(10, 2))

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
        self.update_active_geometry_label()

    def update_active_geometry_label(self):
        mode = self.measure_mode.get()
        if mode in ("direct_basic", "reflected_basic"):
            self.active_geometry_label.set("Active geometry: BASIC")
        elif mode in ("extension_calibrate", "extension_measure"):
            self.active_geometry_label.set("Active geometry: EXTENSION")
        else:
            self.active_geometry_label.set("Active geometry: --")

    def update_calculations(self):
        mode = self.measure_mode.get()

        if mode in ("direct_basic", "reflected_basic"):
            direct = self.direct_point
            reflected = self.reflected_point
        elif mode in ("extension_calibrate", "extension_measure"):
            direct = self.extension_direct_inferred
            reflected = self.extension_reflected_current
        else:
            direct = None
            reflected = None

        if direct is None or reflected is None:
            self.twotheta_label.set("2Theta: --")
            self.tilt_label.set("Tilt: --")
            return

        x0, y0 = direct
        x1, y1 = reflected

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
        self.marker_screen = None
        self.show_markers = False
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

        if self.show_markers:
            mode = self.measure_mode.get()

            if mode == "direct_basic":
                if self.direct_point is not None:
                    self.draw_marker(self.direct_point, color="lime")

            elif mode == "reflected_basic":
                if self.reflected_point is not None:
                    self.draw_marker(self.reflected_point, color="red")

            elif mode == "extension_calibrate":
                if self.extension_reflected_reference is not None:
                    self.draw_marker(self.extension_reflected_reference, color="orange")

            elif mode == "extension_measure":
                if self.extension_reflected_current is not None:
                    self.draw_marker(self.extension_reflected_current, color="lime")

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

        mode = self.measure_mode.get()

        if mode == "reflected_basic" and self.direct_point is None:
            self.status_label.set("Basic reflected measurement requires direct beam position first.")
            return

        if mode == "extension_calibrate" and (self.direct_point is None or self.reflected_point is None):
            self.status_label.set("Extension calibration requires basic direct and reflected beam positions first.")
            return

        if mode == "extension_measure" and self.extension_direct_inferred is None:
            self.status_label.set("Extension measurement requires extension calibration first.")
            return

        if not self.is_valid_click_region(point):
            self.marker_screen = None
            self.status_label.set("Selection canceled: no peak is found.")
            self.redraw_image()
            return

        com_point = self.compute_com(point)
        if com_point is None:
            self.marker_screen = None
            self.status_label.set("Selection canceled: could not compute center of mass.")
            self.redraw_image()
            return
        
        com_point = (round(com_point[0], 1), round(com_point[1], 1))

        self.show_markers = True
        self.marker_screen = com_point
        mode = self.measure_mode.get()

        if mode == "direct_basic":
            self.direct_point = com_point
            self.direct_label.set(f"Direct beam: x={com_point[0]:.2f}, y={com_point[1]:.2f}")
            self.status_label.set("Direct beam updated")

        elif mode == "reflected_basic":
            self.reflected_point = com_point
            self.reflected_label.set(f"Basic reflected beam: x={com_point[0]:.2f}, y={com_point[1]:.2f}")
            self.status_label.set("Basic reflected beam updated")

        elif mode == "extension_calibrate":
            self.extension_reflected_reference = com_point
            self.extension_ref_label.set(
                f"Extension reflected ref: x={com_point[0]:.2f}, y={com_point[1]:.2f}"
            )

            dx = com_point[0] - self.reflected_point[0]
            dy = com_point[1] - self.reflected_point[1]
            self.extension_shift = (dx, dy)
            self.extension_shift_label.set(f"Extension shift: dx={dx:.2f}, dy={dy:.2f}")

            inferred_x = self.direct_point[0] + dx
            inferred_y = self.direct_point[1] + dy
            self.extension_direct_inferred = (inferred_x, inferred_y)
            self.extension_direct_label.set(
                f"Inferred direct beam: x={inferred_x:.2f}, y={inferred_y:.2f}"
            )

            self.extension_reflected_current = com_point
            self.extension_current_label.set(
                f"Extension reflected current: x={com_point[0]:.2f}, y={com_point[1]:.2f}"
            )

            self.status_label.set("Extension shift calibrated")

        elif mode == "extension_measure":
            self.extension_reflected_current = com_point
            self.extension_current_label.set(
                f"Extension reflected current: x={com_point[0]:.2f}, y={com_point[1]:.2f}"
            )
            self.status_label.set("Extension reflected beam updated")

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

        self.extension_reflected_reference = None
        self.extension_reflected_current = None
        self.extension_direct_inferred = None
        self.extension_shift = None

        self.direct_label.set("Direct beam: --")
        self.reflected_label.set("Basic reflected beam: --")
        self.extension_ref_label.set("Extension reflected ref: --")
        self.extension_current_label.set("Extension reflected current: --")
        self.extension_shift_label.set("Extension shift: --")
        self.extension_direct_label.set("Inferred direct beam: --")

        self.status_label.set("Click on image to define beam position")
        self.update_calculations()
        self.redraw_image()


def main():
    root = tk.Tk()
    DexelaMatplotlibGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
