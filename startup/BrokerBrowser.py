from tkinter import ttk, messagebox
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class BrokerBrowser(tk.Tk):
    def __init__(self, db):
        super().__init__()
        self.title("Databroker Browser")
        self.geometry("1400x750")
        self.db = db
        self.results = []
        self.current_hdr = None
        self.current_table = None

        # ========= Main Paned Window (3 columns) =========
        paned = ttk.Panedwindow(self, orient="horizontal")
        paned.pack(fill="both", expand=True)

        # Left pane (1/4 width)
        left = ttk.Frame(paned, padding=5)
        paned.add(left, weight=1)

        # Middle pane (1/4 width)
        middle = ttk.Frame(paned, padding=10)
        paned.add(middle, weight=1)

        # Right pane (1/2 width) – plot window
        right = ttk.Frame(paned, padding=10)
        paned.add(right, weight=2)

        # ========= LEFT PANE =========
        ttk.Label(left, text="Search type:").pack(anchor="w")
        self.query_type = tk.StringVar(value="metadata")
        ttk.Combobox(left, textvariable=self.query_type,
                     values=["metadata", "scan_id", "uid"]).pack(fill="x")

        ttk.Label(left, text="Search term:").pack(anchor="w")
        self.query_entry = ttk.Entry(left)
        self.query_entry.pack(fill="x")

        ttk.Button(left, text="Search", command=self.search_runs).pack(pady=4)

        # Treeview with only 2 columns
        self.tree = ttk.Treeview(
            left,
            columns=("uid", "scan_id"),
            show="headings",
            height=20
        )
        for col in ("uid", "scan_id"):
            self.tree.heading(col, text=col)
            self.tree.column(col, width=160 if col == "uid" else 80)
        self.tree.pack(fill="both", expand=True)
        self.tree.bind("<<TreeviewSelect>>", self.show_details)

        # ========= MIDDLE PANE (motors, detectors) =========
        self.plan_label = ttk.Label(middle, text="Plan: ---")
        self.plan_label.pack(anchor="w", pady=5)

        ttk.Label(middle, text="Motors:").pack(anchor="w")
        self.motor_list = tk.Listbox(middle, height=6)
        self.motor_list.pack(fill="x", pady=4)
        self.motor_list.bind("<<ListboxSelect>>", self.update_plot)

        ttk.Label(middle, text="Detectors:").pack(anchor="w")
        self.det_list = tk.Listbox(middle, selectmode="extended", height=10)
        self.det_list.pack(fill="x", pady=4)
        self.det_list.bind("<<ListboxSelect>>", self.update_plot)

        self.points_label = ttk.Label(middle, text="Number of points: ---")
        self.points_label.pack(anchor="w", pady=5)

        # ========= RIGHT PANE (matplotlib plot) =========
        self.fig, self.ax = plt.subplots(figsize=(5, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    # ===========================================
    # SEARCH
    # ===========================================
    def search_runs(self):
        qtype = self.query_type.get()
        query = self.query_entry.get().strip()

        self.tree.delete(*self.tree.get_children())
        self.results = []

        try:
            if qtype == "metadata":
                key, value = query.split("=")
                headers = self.db(**{key.strip(): value.strip()})
            elif qtype == "scan_id":
                headers = self.db(scan_id=int(query))
            elif qtype == "uid":
                headers = [self.db[query]]
            else:
                headers = []

            for hdr in headers:
                uid = hdr.start.get("uid", "?")
                scan_id = hdr.start.get("scan_id", "?")
                self.results.append(hdr)
                self.tree.insert("", "end", values=(uid, scan_id))

        except Exception as e:
            messagebox.showerror("Error", str(e))

    # ===========================================
    # DETAILS PANEL
    # ===========================================
    def show_details(self, event):
        sel = self.tree.selection()
        if not sel:
            return

        uid = self.tree.item(sel[0], "values")[0]
        try:
            hdr = self.db[uid]
            self.current_hdr = hdr
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load: {e}")
            return

        # Begin data access outside Tk calls
        plan = hdr.start.get("plan_name", "unknown")
        motors = hdr.start.get("motors", [])
        detectors = hdr.start.get("detectors", [])
        num_points = hdr.start.get("num_points", "?")

        # Load the table, which may be heavy
        try:
            self.current_table = hdr.table()
        except Exception as e:
            messagebox.showerror("Error", f"Table error: {e}")
            return

        # Tk updates
        self.plan_label.config(text=f"Plan: {plan}")
        self.points_label.config(text=f"Number of points: {num_points}")

        # Motors
        self.motor_list.delete(0, tk.END)
        for m in motors:
            self.motor_list.insert(tk.END, m)
        if motors:
            self.motor_list.selection_set(0)
            self.motor_list.activate(0)

        # Detectors
        self.det_list.delete(0, tk.END)
        if detectors:
            # Detector columns are expanded names inside table
            root = detectors[0]
            tbl_cols = self.current_table.columns
            det_cols = [c for c in tbl_cols if root in c]

            for d in det_cols:
                self.det_list.insert(tk.END, d)

            # Default: last detector selected
            if det_cols:
                last = len(det_cols) - 1
                self.det_list.selection_set(last)
                self.det_list.activate(last)

        self.update_plot()

    # ===========================================
    # PLOTTING
    # ===========================================
    def update_plot(self, event=None):
        if self.current_table is None:
            return

        breakpoint()
        # Selected motor
        motor_sel = self.motor_list.curselection()
        if not motor_sel:
            return
        motor = self.motor_list.get(motor_sel[0])

        # Selected detectors (multiple)
        det_indices = self.det_list.curselection()
        if not det_indices:
            return
        detectors = [self.det_list.get(i) for i in det_indices]

        tab = self.current_table
        if motor not in tab.columns:
            return

        x = tab[motor].values

        # Clear plot
        self.ax.clear()

        # Cycle through detectors with distinct colors
        colors = plt.cm.tab10.colors
        for i, det in enumerate(detectors):
            if det not in tab.columns:
                continue
            y = tab[det].values
            self.ax.plot(x, y, label=det, color=colors[i % len(colors)])

        self.ax.set_xlabel(motor)
        self.ax.set_ylabel("Detector signal")
        self.ax.legend()
        self.ax.grid(True)

        self.canvas.draw()


def OpenBrokerBrowser():
# Opens the Databroker Browser application.
    while True:
        app = BrokerBrowser(db)  # pass your existing databroker here
        app.mainloop()  # runs until window is closed
        # ask user if they want to reopen
        answer = input("Open the app again? [y/n]: ").strip().lower()
        if answer != 'y':
            break
