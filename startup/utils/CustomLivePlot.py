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
    Live 2D plot for selected detector channels on one persistent axes.

    Parameters
    ----------
    y_fields : list[str]
        Event-data field names to plot.
    x : str
        Event-data field name for x-axis.
    ax : matplotlib.axes.Axes
        Existing persistent axes. Required.
    legend_keys : dict[str, str], optional
        Map field name -> legend label.
    stream_name : str, optional
        Stream name to listen to, default 'primary'.
    clear_on_start : bool, optional
        If True, clear axes at the start of each run.
    show_stats : bool, optional
        If True, show COM and FWHM in a text box.
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
        ax,
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

        if ax is None:
            raise ValueError("ax must be a persistent matplotlib Axes")

        self.y_fields = list(y_fields)
        self.x_field = x
        self.ax = ax
        self.legend_keys = legend_keys or {}
        self.stream_name = stream_name
        self.clear_on_start = clear_on_start
        self.show_stats = show_stats
        self.update_every = max(1, int(update_every))
        self.xlabel = xlabel or x
        self.ylabel = ylabel
        self.title = title

        self._fig = self.ax.figure
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

        if self.clear_on_start:
            self.ax.cla()

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
        self._fig.canvas.draw_idle()
        self._fig.canvas.flush_events()
        # self._fig.show()

def select_detector_fields(det, channels=[0]):
    hinted = list(det.hints["fields"])
    return [hinted[i] for i in channels]


class CustomLiveMesh(QtAwareCallback):
    """
    Live 2-D scatter plot for selected detector channels during a mesh scan.

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
    fig : matplotlib.figure.Figure
        Persistent figure (e.g. ``mymeshfig``).  Required.
    legend_keys : dict[str, str], optional
        Map field name -> subplot title / label.
    stream_name : str, optional
        Stream to listen to.  Default ``'primary'``.
    clear_on_start : bool, optional
        If True (default), call ``fig.clf()`` at the start of each scan.
    update_every : int, optional
        Redraw every N accepted events.  Default 10.
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
        *,
        fig,
        legend_keys=None,
        stream_name="primary",
        clear_on_start=True,
        update_every=10,
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
        self.fig = fig
        self.legend_keys = legend_keys or {}
        self.stream_name = stream_name
        self.clear_on_start = clear_on_start
        self.update_every = max(1, int(update_every))
        self.cmap = cmap
        self.title = title

        self._descriptor_uids = set()
        self._inner_data = {}
        self._outer_data = {}
        self._z_data = {}
        self._scatters = {}
        self._axes = {}
        self._colorbars = {}
        self._event_count = 0

    # ------------------------------------------------------------------
    def start(self, doc):
        self._descriptor_uids.clear()
        self._inner_data = {f: [] for f in self.y_fields}
        self._outer_data = {f: [] for f in self.y_fields}
        self._z_data = {f: [] for f in self.y_fields}
        self._scatters = {}
        self._axes = {}
        self._colorbars = {}
        self._event_count = 0

        if self.clear_on_start:
            self.fig.clf()

        n = max(len(self.y_fields), 1)
        for i, field in enumerate(self.y_fields):
            ax = self.fig.add_subplot(1, n, i + 1)
            self._axes[field] = ax

            label = self.legend_keys.get(field, field)
            sc = ax.scatter([], [], c=np.array([]), cmap=self.cmap, s=4)
            sc.set_clim(0, 1)          # prevent warnings on empty scatter
            self._scatters[field] = sc

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

        data = doc.get("data", {})
        if self.inner_field not in data or self.outer_field not in data:
            return

        inner_val = data[self.inner_field]
        outer_val = data[self.outer_field]

        for field in self.y_fields:
            self._inner_data[field].append(inner_val)
            self._outer_data[field].append(outer_val)
            self._z_data[field].append(data.get(field, np.nan))

        self._event_count += 1

        if self._event_count % self.update_every == 0:
            self._update_plots()
            self._draw()

    # ------------------------------------------------------------------
    def stop(self, doc):
        self._update_plots()

        # Add colorbars on scan completion (once per field)
        for field in self.y_fields:
            if field not in self._colorbars and field in self._scatters:
                sc = self._scatters[field]
                ax = self._axes[field]
                try:
                    cb = self.fig.colorbar(sc, ax=ax)
                    self._colorbars[field] = cb
                except Exception:
                    pass

        self._draw()

    # ------------------------------------------------------------------
    def _update_plots(self):
        for field in self.y_fields:
            inner = np.asarray(self._inner_data[field], dtype=float)
            outer = np.asarray(self._outer_data[field], dtype=float)
            z = np.asarray(self._z_data[field], dtype=float)

            if len(inner) == 0:
                continue

            sc = self._scatters[field]
            ax = self._axes[field]

            sc.set_offsets(np.c_[inner, outer])
            sc.set_array(z)

            valid_z = z[np.isfinite(z)]
            if len(valid_z) >= 1:
                zmin, zmax = float(valid_z.min()), float(valid_z.max())
                if zmin < zmax:
                    sc.set_clim(zmin, zmax)

            def _lim(arr):
                lo, hi = float(arr.min()), float(arr.max())
                span = hi - lo
                margin = 0.05 * span if span > 0 else 0.5
                return lo - margin, hi + margin

            ax.set_xlim(*_lim(inner))
            ax.set_ylim(*_lim(outer))

    # ------------------------------------------------------------------
    def _draw(self):
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()