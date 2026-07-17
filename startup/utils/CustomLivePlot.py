import numpy as np
import matplotlib.pyplot as plt

from bluesky.callbacks import CallbackBase
from bluesky.callbacks.mpl_plotting import QtAwareCallback


def _calc_com(x, y):
    """
    Center of mass using only non-negative finite signal values.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        return np.nan

    x = x[m]
    y = y[m]

    y = np.clip(y, 0, None)
    s = y.sum()
    if s <= 0:
        return np.nan

    return np.sum(x * y) / s


def _calc_fwhm(x, y):
    """
    Approximate FWHM from sampled data using linear interpolation
    around the half-maximum crossings.

    Returns
    -------
    (fwhm, x_center)
        fwhm     : width at half maximum
        x_center : midpoint between left/right half-max crossings

    Returns (np.nan, np.nan) if the curve does not cross half-maximum twice.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(m) < 3:
        return np.nan, np.nan

    x = x[m]
    y = y[m]

    ymax = np.nanmax(y)
    if not np.isfinite(ymax) or ymax <= 0:
        return np.nan, np.nan

    half = ymax / 2.0
    above = y >= half
    idx = np.where(above)[0]

    if len(idx) < 2:
        return np.nan, np.nan

    left_i = idx[0]
    right_i = idx[-1]

    # left crossing
    if left_i == 0:
        x_left = x[0]
    else:
        x0, x1 = x[left_i - 1], x[left_i]
        y0, y1 = y[left_i - 1], y[left_i]
        if y1 == y0:
            x_left = x1
        else:
            x_left = x0 + (half - y0) * (x1 - x0) / (y1 - y0)

    # right crossing
    if right_i == len(x) - 1:
        x_right = x[-1]
    else:
        x0, x1 = x[right_i], x[right_i + 1]
        y0, y1 = y[right_i], y[right_i + 1]
        if y1 == y0:
            x_right = x0
        else:
            x_right = x0 + (half - y0) * (x1 - x0) / (y1 - y0)

    fwhm = x_right - x_left
    x_center = 0.5 * (x_left + x_right)

    return fwhm, x_center

def _calc_max(x, y):
    """
    Return (x_at_max, y_max).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        return np.nan, np.nan

    x = x[m]
    y = y[m]

    if len(y) == 0:
        return np.nan, np.nan

    i = np.argmax(y)
    return x[i], y[i]

class CustomLivePlot(QtAwareCallback):
    """
    Live 1-D plot for selected detector channels on a shared persistent figure.

    The axes are rebuilt at the start of every scan (``fig.clf()`` +
    ``fig.add_subplot``), so the figure window can be shared with
    ``CustomLiveMesh`` without leaving stale axes references if the window
    is closed and reopened between scans.

    Parameters
    ----------
    y_fields : list[str]
        Event-data field names to plot.
    x : str
        Event-data field name for x-axis.
    fig : matplotlib.figure.Figure
        Persistent figure shared with mesh scans (e.g. ``myfig``). Required.
    legend_keys : dict[str, str], optional
        Map field name -> legend label.
    stream_name : str, optional
        Stream name to listen to, default 'primary'.
    clear_on_start : bool, optional
        If True, clear the figure at the start of each run.
    show_stats : bool, optional
        If True, show MAX and FWHM in a text box.
    update_every : int, optional
        Redraw every N accepted events. Default 1.
    xlabel : str, optional
        Override x-axis label.
    ylabel : str, optional
        Override y-axis label.
    title : str, optional
        Plot title.
    """
    def __init__(
        self,
        y_fields,
        x,
        *,
        fig,
        legend_keys=None,
        stream_name="primary",
        clear_on_start=True,
        show_stats=True,
        update_every=1,
        xlabel=None,
        ylabel="signal",
        title="Live detector channels",
        use_teleporter=None,
    ):
        super().__init__(use_teleporter=use_teleporter)

        if fig is None:
            raise ValueError("fig must be a persistent matplotlib Figure")

        self.y_fields = list(y_fields)
        self.x_field = x
        self.fig = fig
        self._fig_num = fig.get_label()
        self._fig_size = tuple(fig.get_size_inches())
        self.ax = None
        self.legend_keys = legend_keys or {}
        self.stream_name = stream_name
        self.clear_on_start = clear_on_start
        self.show_stats = show_stats
        self.update_every = max(1, int(update_every))
        self.xlabel = xlabel or x
        self.ylabel = ylabel
        self.title = title

        self._descriptor_uids = set()
        self._xdata = []
        self._ydata = {field: [] for field in self.y_fields}
        self._lines = {}
        self._stats_text = None
        self._event_count = 0
        self._run_uid = None

    def start(self, doc):
        self._run_uid = doc.get("uid")
        self._descriptor_uids.clear()
        self._xdata = []
        self._ydata = {field: [] for field in self.y_fields}
        self._lines = {}
        self._stats_text = None
        self._event_count = 0

        # Recreate the figure window if it was closed since the last scan.
        self.fig = plt.figure(num=self._fig_num, figsize=self._fig_size, clear=False)

        if self.clear_on_start:
            self.fig.clf()
        self.ax = self.fig.add_subplot(1, 1, 1)

        for field in self.y_fields:
            label = self.legend_keys.get(field, field)
            line, = self.ax.plot([], [], marker="o", linestyle="-", label=label)
            self._lines[field] = line

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)

        self.ax.set_autoscale_on(True)
        self.ax.margins(y=0.15)

        if self.y_fields:
            self.ax.legend(loc="best")

        if self.show_stats:
            self._stats_text = self.ax.text(
                0.02,
                0.98,
                "",
                transform=self.ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
            )

        self._draw()

    def descriptor(self, doc):
        if doc.get("name") == self.stream_name:
            self._descriptor_uids.add(doc["uid"])

    def event(self, doc):
        if doc.get("descriptor") not in self._descriptor_uids:
            return

        data = doc.get("data", {})
        if self.x_field not in data:
            return

        x = data[self.x_field]
        self._xdata.append(x)

        for field in self.y_fields:
            self._ydata[field].append(data.get(field, np.nan))

        self._event_count += 1

        xarr = np.asarray(self._xdata, dtype=float)
        stats_lines = []

        for field in self.y_fields:
            yarr = np.asarray(self._ydata[field], dtype=float)
            self._lines[field].set_data(xarr, yarr)

            if self.show_stats:
                # com = _calc_com(xarr, yarr)
                fwhm, xfwhm = _calc_fwhm(xarr, yarr)
                xmax, ymax = _calc_max(xarr, yarr)

                label = self.legend_keys.get(field, field)
                stats_lines.append(f"{label}: MAX=({xmax:.5g}, {ymax:.5g}), FWHM={fwhm:.5g} at {xfwhm:.5g}")

        if self.show_stats and self._stats_text is not None:
            self._stats_text.set_text("\n".join(stats_lines))

        if self._event_count % self.update_every == 0:
            self.ax.relim()
            self.ax.autoscale_view()
            self._draw()

    def stop(self, doc):
        self.ax.relim()
        self.ax.autoscale_view()
        self._draw()

    def _draw(self):
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        # self.fig.show()

def select_detector_fields(det, channels=[0]):
    hinted = list(det.hints["fields"])
    return [hinted[i] for i in channels]


class CustomLiveMesh(QtAwareCallback):
    """
    Live 2-D image plot for selected detector channels during a mesh scan.

    Uses ``ax.imshow`` so each scan cell is a solid-filled rectangle with no
    gaps between points — identical in appearance to Bluesky's ``LiveGrid``.
    One subplot is created per field in *y_fields*.  The subplot layout is
    rebuilt each scan (``fig.clf()`` is called in ``start()``), so the number
    of subplots adapts dynamically.

    Parameters
    ----------
    y_fields : list[str]
        Event-data field names to plot (one subplot each).
    outer_field : str
        Event-data field name for the slow (outer) motor — plotted on the
        y-axis.
    inner_field : str
        Event-data field name for the fast (inner) motor — plotted on the
        x-axis.
    slow_steps : int
        Number of points on the slow (outer) axis — sets the number of image
        rows.
    fast_steps : int
        Number of points on the fast (inner) axis — sets the number of image
        columns.
    fig : matplotlib.figure.Figure
        Persistent figure shared with 1-D scans (e.g. ``myfig``).  Required.
    legend_keys : dict[str, str], optional
        Map field name -> subplot title / label.
    stream_name : str, optional
        Stream to listen to.  Default ``'primary'``.
    clear_on_start : bool, optional
        If True (default), call ``fig.clf()`` at the start of each scan.
    update_every : int, optional
        Redraw every N accepted events.  Default 1.
    cmap : str, optional
        Matplotlib colormap name.  Default ``'viridis'``.
    title : str, optional
        Figure suptitle.
    """

    def __init__(
        self,
        y_fields,
        outer_field,
        inner_field,
        slow_steps,
        fast_steps,
        *,
        fig,
        legend_keys=None,
        stream_name="primary",
        clear_on_start=True,
        update_every=1,
        cmap="viridis",
        title="Live Mesh Scan",
        use_teleporter=None,
    ):
        super().__init__(use_teleporter=use_teleporter)

        if fig is None:
            raise ValueError("fig must be a persistent matplotlib Figure")

        self.y_fields = list(y_fields)
        self.outer_field = outer_field
        self.inner_field = inner_field
        self.slow_steps = max(1, int(slow_steps))
        self.fast_steps = max(1, int(fast_steps))
        self.fig = fig
        self._fig_num = fig.get_label()
        self._fig_size = tuple(fig.get_size_inches())
        self.legend_keys = legend_keys or {}
        self.stream_name = stream_name
        self.clear_on_start = clear_on_start
        self.update_every = max(1, int(update_every))
        self.cmap = cmap
        self.title = title

        self._descriptor_uids = set()
        self._Z = {}
        self._images = {}
        self._axes = {}
        self._colorbars = {}
        self._inner_vals = []
        self._outer_vals = []
        self._event_count = 0

    # ------------------------------------------------------------------
    def start(self, doc):
        self._descriptor_uids.clear()
        self._Z = {}
        self._images = {}
        self._axes = {}
        self._colorbars = {}
        self._inner_vals = []
        self._outer_vals = []
        self._event_count = 0

        # Recreate the figure window if it was closed since the last scan.
        self.fig = plt.figure(num=self._fig_num, figsize=self._fig_size, clear=False)

        if self.clear_on_start:
            self.fig.clf()

        # Colormap with a distinct colour for not-yet-acquired cells.
        cmap_obj = plt.get_cmap(self.cmap).copy()
        cmap_obj.set_bad("lightgray")

        n = max(len(self.y_fields), 1)
        for i, field in enumerate(self.y_fields):
            ax = self.fig.add_subplot(1, n, i + 1)
            self._axes[field] = ax

            Z_init = np.full((self.slow_steps, self.fast_steps), np.nan)
            self._Z[field] = Z_init

            label = self.legend_keys.get(field, field)
            im = ax.imshow(
                Z_init,
                origin="lower",           # row 0 at bottom (Cartesian convention)
                aspect="auto",            # stretch to fill axes
                interpolation="nearest",  # sharp cell boundaries, no blending
                cmap=cmap_obj,
                vmin=0,
                vmax=1,                   # placeholder; replaced on first real update
            )
            self._images[field] = im

            ax.set_xlabel(self.inner_field)
            ax.set_ylabel(self.outer_field)
            ax.set_title(label)

        self.fig.suptitle(self.title)
        self._draw()

    # ------------------------------------------------------------------
    def descriptor(self, doc):
        if doc.get("name") == self.stream_name:
            self._descriptor_uids.add(doc["uid"])

    # ------------------------------------------------------------------
    def event(self, doc):
        if doc.get("descriptor") not in self._descriptor_uids:
            return

        # Use the 1-indexed seq_num to determine the grid cell.
        # seq_num counts primary-stream events from 1 regardless of other streams.
        idx = doc.get("seq_num", 1) - 1
        row = idx // self.fast_steps
        col = idx % self.fast_steps

        if row >= self.slow_steps or col >= self.fast_steps:
            return   # safety guard against out-of-range events

        data = doc.get("data", {})

        for field in self.y_fields:
            self._Z[field][row, col] = data.get(field, np.nan)

        # Collect actual motor positions for axis extent computation.
        inner_val = data.get(self.inner_field)
        outer_val = data.get(self.outer_field)
        if inner_val is not None:
            self._inner_vals.append(float(inner_val))
        if outer_val is not None:
            self._outer_vals.append(float(outer_val))

        self._event_count += 1

        if self._event_count % self.update_every == 0:
            self._update_plots()
            self._draw()

    # ------------------------------------------------------------------
    def stop(self, doc):
        self._update_plots()

        # Add colorbars on scan completion (once per field).
        for field in self.y_fields:
            if field not in self._colorbars and field in self._images:
                im = self._images[field]
                ax = self._axes[field]
                try:
                    cb = self.fig.colorbar(im, ax=ax)
                    self._colorbars[field] = cb
                except Exception:
                    pass

        self._draw()

    # ------------------------------------------------------------------
    def _compute_extent(self):
        """Return [x0, x1, y0, y1] from collected motor positions, or None."""
        iv = [v for v in self._inner_vals if np.isfinite(v)]
        ov = [v for v in self._outer_vals if np.isfinite(v)]
        if not iv or not ov:
            return None

        x0, x1 = min(iv), max(iv)
        y0, y1 = min(ov), max(ov)
        dx = (x1 - x0) / (self.fast_steps - 1) if (self.fast_steps > 1 and x1 > x0) else 1.0
        dy = (y1 - y0) / (self.slow_steps - 1) if (self.slow_steps > 1 and y1 > y0) else 1.0

        # Extent uses cell edges, not centres (half-step padding on each side).
        return [x0 - dx / 2, x1 + dx / 2, y0 - dy / 2, y1 + dy / 2]

    # ------------------------------------------------------------------
    def _update_plots(self):
        extent = self._compute_extent()

        for field in self.y_fields:
            if field not in self._images:
                continue

            im = self._images[field]
            ax = self._axes[field]
            Z  = self._Z[field]

            im.set_data(Z)

            valid = Z[np.isfinite(Z)]
            if len(valid) >= 2 and valid.min() < valid.max():
                im.set_clim(float(valid.min()), float(valid.max()))

            if extent is not None:
                im.set_extent(extent)
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])

    # ------------------------------------------------------------------
    def _draw(self):
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()