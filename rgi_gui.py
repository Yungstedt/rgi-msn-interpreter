"""
RG-I MS^n Interpreter -- desktop front end for rgi_msn.py.

Standard library only (Tkinter + csv), so it packages into a small single-file
executable and runs on any machine with no scientific Python stack installed.

Run from source:   python rgi_gui.py
"""

from __future__ import annotations

import csv
import os
import queue
import sys
import threading
import tkinter as tk
import traceback
from tkinter import filedialog, messagebox, ttk

import rgi_msn as R

APP = "RG-I MS\u207f Interpreter"
PRESETS = {
    "DP4  (celery lyase digest)": "dGalA,Rha,GalA,Rha",
    "DP6  (A. mucilage lyase digest)": "dGalA,Rha,GalA,Rha,GalA,Rha",
    "DP8": "dGalA,Rha,GalA,Rha,GalA,Rha,GalA,Rha",
    "DP10": "dGalA,Rha,GalA,Rha,GalA,Rha,GalA,Rha,GalA,Rha",
}
DERIVS = {"trideuteroacetyl  (+45.0294)": "dAc",
          "propionyl  (+56.0262)": "Prop",
          "acetyl  (+42.0106)": "Ac",
          "methyl  (+14.0157)": "Me",
          "none (underivatized)": "OH"}
ADDUCTS = {"[M+Na]+": "Na", "[M+H]+": "H", "[M+K]+": "K", "[M-H]-": "-H"}


# --------------------------------------------------------------- spectrum ---

class SpectrumCanvas(tk.Canvas):
    """Stick spectrum: predicted fragments up, observed peaks down."""

    def __init__(self, master, **kw):
        super().__init__(master, bg="white", highlightthickness=1,
                         highlightbackground="#c8c8c8", **kw)
        self.pred: list = []
        self.obs: list = []
        self.bind("<Configure>", lambda e: self.redraw())

    def set_data(self, pred, obs):
        self.pred, self.obs = list(pred), list(obs)
        self.redraw()

    def redraw(self):
        self.delete("all")
        W, H = self.winfo_width(), self.winfo_height()
        if W < 50 or H < 50:
            return
        L, Rm, T, B = 46, 14, 14, 26
        mid = (T + H - B) / 2
        xs = [m for m, _, _ in self.pred] + [m for m, _ in self.obs]
        if not xs:
            self.create_text(W / 2, H / 2, text="no fragments computed yet",
                             fill="#999999", font=("Segoe UI", 9))
            return
        lo, hi = min(xs), max(xs)
        if hi - lo < 1:
            lo, hi = lo - 1, hi + 1
        pad = (hi - lo) * 0.04
        lo, hi = lo - pad, hi + pad

        def X(m):
            return L + (m - lo) / (hi - lo) * (W - L - Rm)

        self.create_line(L, mid, W - Rm, mid, fill="#000000")
        for i in range(6):
            m = lo + (hi - lo) * i / 5
            x = X(m)
            self.create_line(x, H - B, x, H - B + 4, fill="#666666")
            self.create_text(x, H - B + 13, text=f"{m:.0f}",
                             fill="#444444", font=("Segoe UI", 7))
        self.create_text(L - 6, T + 6, text="predicted", anchor="w",
                         fill="#1f6fb4", font=("Segoe UI", 7))
        self.create_text(L - 6, H - B - 6, text="observed", anchor="w",
                         fill="#c4531a", font=("Segoe UI", 7))

        half = mid - T - 10
        for m, lab, hot in self.pred:
            x = X(m)
            self.create_line(x, mid, x, mid - (half if hot else half * 0.55),
                             fill="#1f6fb4" if hot else "#b9c7d4",
                             width=2 if hot else 1)
            if hot:
                self.create_text(x, mid - half - 4, text=lab, angle=90,
                                 anchor="w", fill="#1f6fb4",
                                 font=("Segoe UI", 7))
        if self.obs:
            top = max(i for _, i in self.obs) or 1.0
            for m, inten in self.obs:
                x = X(m)
                h = (mid - T - 10) * (inten / top if top else 1.0)
                self.create_line(x, mid, x, mid + max(h, 3), fill="#c4531a",
                                 width=2)


# -------------------------------------------------------------------- app ---

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP)
        self.geometry("1140x760")
        self.minsize(960, 640)
        self.q: queue.Queue = queue.Queue()
        self.rows: list = []
        self.headers: list = []
        self._build()
        self.after(100, self._drain)

    # ---- layout --------------------------------------------------------

    def _build(self):
        try:
            ttk.Style().theme_use("vista")
        except tk.TclError:
            pass
        root = ttk.Frame(self, padding=8)
        root.pack(fill="both", expand=True)
        left = ttk.Frame(root)
        left.pack(side="left", fill="y", padx=(0, 8))
        right = ttk.Frame(root)
        right.pack(side="left", fill="both", expand=True)

        # --- structure ---
        g = ttk.LabelFrame(left, text="Oligomer", padding=8)
        g.pack(fill="x")
        ttk.Label(g, text="Preset").grid(row=0, column=0, sticky="w")
        self.preset = ttk.Combobox(g, values=list(PRESETS), width=28,
                                   state="readonly")
        self.preset.current(1)
        self.preset.grid(row=0, column=1, sticky="ew", pady=2)
        self.preset.bind("<<ComboboxSelected>>", self._apply_preset)
        ttk.Label(g, text="Backbone").grid(row=1, column=0, sticky="w")
        self.backbone = ttk.Entry(g, width=30)
        self.backbone.insert(0, PRESETS["DP6  (A. mucilage lyase digest)"])
        self.backbone.grid(row=1, column=1, sticky="ew", pady=2)
        ttk.Label(g, text="non-reducing \u2192 reducing; dGalA = \u03944,5-unsaturated",
                  foreground="#777777", font=("Segoe UI", 7)
                  ).grid(row=2, column=0, columnspan=2, sticky="w")
        ttk.Label(g, text="Blocking group").grid(row=3, column=0, sticky="w")
        self.deriv = ttk.Combobox(g, values=list(DERIVS), width=28,
                                  state="readonly")
        self.deriv.current(0)
        self.deriv.grid(row=3, column=1, sticky="ew", pady=2)
        ttk.Label(g, text="Native O-acetyls").grid(row=4, column=0, sticky="w")
        self.noac = tk.Spinbox(g, from_=0, to=4, width=5)
        self.noac.delete(0, "end")
        self.noac.insert(0, "1")
        self.noac.grid(row=4, column=1, sticky="w", pady=2)
        ttk.Label(g, text="Gal branch").grid(row=5, column=0, sticky="w")
        self.branch = ttk.Entry(g, width=30)
        self.branch.grid(row=5, column=1, sticky="ew", pady=2)
        ttk.Label(g, text="blank, or residue index on Rha O-4, e.g.  1",
                  foreground="#777777", font=("Segoe UI", 7)
                  ).grid(row=6, column=0, columnspan=2, sticky="w")
        ttk.Label(g, text="Adduct").grid(row=7, column=0, sticky="w")
        self.adduct = ttk.Combobox(g, values=list(ADDUCTS), width=12,
                                   state="readonly")
        self.adduct.current(0)
        self.adduct.grid(row=7, column=1, sticky="w", pady=2)
        g.columnconfigure(1, weight=1)

        # --- acquisition ---
        a = ttk.LabelFrame(left, text="Acquisition", padding=8)
        a.pack(fill="x", pady=(8, 0))
        ttk.Label(a, text="MS\u207f precursor path").grid(row=0, column=0, sticky="w")
        self.path = ttk.Entry(a, width=18)
        self.path.grid(row=0, column=1, sticky="ew", pady=2)
        ttk.Label(a, text="blank = MS\u00b2; '1087.43' = MS\u00b3 of that ion",
                  foreground="#777777", font=("Segoe UI", 7)
                  ).grid(row=1, column=0, columnspan=2, sticky="w")
        ttk.Label(a, text="Cleavage depth").grid(row=2, column=0, sticky="w")
        self.depth = tk.Spinbox(a, from_=1, to=3, width=5)
        self.depth.delete(0, "end")
        self.depth.insert(0, "1")
        self.depth.grid(row=2, column=1, sticky="w", pady=2)
        self.xring = tk.BooleanVar(value=True)
        ttk.Checkbutton(a, text="include cross-ring (A/X) cleavages",
                        variable=self.xring).grid(row=3, column=0, columnspan=2,
                                                  sticky="w")
        ttk.Label(a, text="Tolerance (mmu)").grid(row=4, column=0, sticky="w")
        self.tol = ttk.Entry(a, width=8)
        self.tol.insert(0, "5")
        self.tol.grid(row=4, column=1, sticky="w", pady=2)
        a.columnconfigure(1, weight=1)

        # --- peaks ---
        p = ttk.LabelFrame(left, text="Observed peak list", padding=8)
        p.pack(fill="both", expand=True, pady=(8, 0))
        self.peaks = tk.Text(p, width=34, height=11, font=("Consolas", 9))
        self.peaks.pack(fill="both", expand=True)
        ttk.Label(p, text="one peak per line:  m/z  [intensity]",
                  foreground="#777777", font=("Segoe UI", 7)).pack(anchor="w")
        pb = ttk.Frame(p)
        pb.pack(fill="x", pady=(4, 0))
        ttk.Button(pb, text="Load\u2026", command=self.load_peaks).pack(side="left")
        ttk.Button(pb, text="Clear",
                   command=lambda: self.peaks.delete("1.0", "end")).pack(side="left", padx=4)

        # --- actions ---
        b = ttk.Frame(left)
        b.pack(fill="x", pady=8)
        ttk.Button(b, text="Predict fragments", command=self.do_predict).pack(fill="x", pady=1)
        ttk.Button(b, text="Interpret peak list", command=self.do_interpret).pack(fill="x", pady=1)
        ttk.Button(b, text="Minimal mixture", command=self.do_mixture).pack(fill="x", pady=1)
        ttk.Button(b, text="Export table to CSV\u2026", command=self.export).pack(fill="x", pady=1)

        # --- results ---
        self.status = ttk.Label(right, text="ready", foreground="#555555")
        self.status.pack(anchor="w")
        self.canvas = SpectrumCanvas(right, height=190)
        self.canvas.pack(fill="x", pady=(4, 6))
        tf = ttk.Frame(right)
        tf.pack(fill="both", expand=True)
        self.tree = ttk.Treeview(tf, show="headings", height=18)
        vs = ttk.Scrollbar(tf, orient="vertical", command=self.tree.yview)
        hs = ttk.Scrollbar(tf, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vs.set, xscrollcommand=hs.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vs.grid(row=0, column=1, sticky="ns")
        hs.grid(row=1, column=0, sticky="ew")
        tf.rowconfigure(0, weight=1)
        tf.columnconfigure(0, weight=1)

    def _apply_preset(self, _evt=None):
        self.backbone.delete(0, "end")
        self.backbone.insert(0, PRESETS[self.preset.get()])

    # ---- form reading --------------------------------------------------

    def cfg(self):
        bb = [t.strip() for t in self.backbone.get().split(",") if t.strip()]
        bad = [t for t in bb if t not in R.RESIDUES]
        if bad:
            raise ValueError("unknown residue(s): %s\nallowed: %s"
                             % (", ".join(bad), ", ".join(R.RESIDUES)))
        if not bb:
            raise ValueError("backbone is empty")
        branches = {}
        txt = self.branch.get().strip()
        if txt:
            i = int(txt)
            if not 0 <= i < len(bb) or bb[i] != "Rha":
                raise ValueError("Gal branch index must point at a Rha residue "
                                 "(0 = non-reducing end)")
            branches[(i, 4)] = "Gal"
        path = [float(t) for t in self.path.get().replace(",", " ").split()]
        return dict(bb=bb, acyl=DERIVS[self.deriv.get()],
                    noac=int(self.noac.get()), branches=branches,
                    adduct=ADDUCTS[self.adduct.get()], path=path,
                    depth=int(self.depth.get()), tol=float(self.tol.get()),
                    xring=bool(self.xring.get()))

    def peak_list(self):
        out = []
        for line in self.peaks.get("1.0", "end").splitlines():
            line = line.replace(",", " ").strip()
            if not line or line[0].isalpha():
                continue
            parts = line.split()
            out.append((float(parts[0]),
                        float(parts[1]) if len(parts) > 1 else 1.0))
        return out

    def candidates(self, c):
        return R.enumerate_isomers(c["bb"], c["noac"], acyl=c["acyl"],
                                   branches=c["branches"])

    # ---- background plumbing -------------------------------------------

    def run(self, fn, note):
        self.status.config(text=note + " \u2026", foreground="#1f6fb4")
        threading.Thread(target=self._worker, args=(fn,), daemon=True).start()

    def _worker(self, fn):
        try:
            self.q.put(("ok", fn()))
        except Exception as exc:                      # surfaced in the UI
            self.q.put(("err", "".join(traceback.format_exception_only(
                type(exc), exc)).strip()))

    def _drain(self):
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind == "err":
                    self.status.config(text="error", foreground="#b00020")
                    messagebox.showerror(APP, payload)
                else:
                    headers, rows, pred, note = payload
                    self.show(headers, rows)
                    self.canvas.set_data(pred, self.peak_list())
                    self.status.config(text=note, foreground="#2b7a2b")
        except queue.Empty:
            pass
        self.after(100, self._drain)

    def show(self, headers, rows):
        self.headers, self.rows = headers, rows
        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = headers
        for h in headers:
            self.tree.heading(h, text=h)
            w = 300 if h in ("label", "assignments", "isomer") else 95
            self.tree.column(h, width=w, anchor="w", stretch=False)
        for r in rows:
            self.tree.insert("", "end", values=r)

    # ---- actions -------------------------------------------------------

    def do_predict(self):
        c = self.cfg()

        def job():
            cands = self.candidates(c)
            hdr = ["isomer", "m/z", "ion", "formula", "residues"]
            rows, pred = [], []
            for s in cands:
                _, frs = R.fragments_for(s, c["path"], c["depth"], c["adduct"],
                                         cross_ring=c["xring"])
                for f in frs:
                    mz = f.mz(c["adduct"])
                    rows.append([s.name, f"{mz:.4f}", R.label(f, s),
                                 R.formula(f), R.residue_span(f)])
                    if s is cands[0]:
                        pred.append((mz, R.label(f, s), len(R.label(f, s).split()) <= 1))
            rows.sort(key=lambda r: (r[0], float(r[1])))
            return hdr, rows, pred, ("%d fragments for %d candidate isomer(s)"
                                     % (len(rows), len(cands)))
        self.run(job, "predicting fragments")

    def do_interpret(self):
        c = self.cfg()
        pk = self.peak_list()
        if not pk:
            messagebox.showwarning(APP, "Enter or load an observed peak list first.")
            return

        def job():
            cands = self.candidates(c)
            res = R.interpret(pk, cands, c["path"], c["depth"], c["adduct"],
                              c["tol"])
            hdr = ["m/z", "intensity", "explained by", "n isomers",
                   "n assignments", "unique?", "best assignment"]
            rows = []
            for row in res.peak_table:
                first = next((v for k, v in row.items()
                              if k.startswith("assign[") and v), "")
                rows.append([f"{row['mz']:.4f}", f"{row['intensity']:g}",
                             row["explained_by"][:60],
                             row["n_candidates_explaining"],
                             row["n_assignments"],
                             "yes" if row["discriminating"] else "",
                             first[:90]])
            rows.append(["", "", "", "", "", "", ""])
            for n, sc in sorted(res.scores.items(),
                                key=lambda kv: -kv[1]["discriminating_fraction"]):
                rows.append(["RANK", n, "", sc["peaks_explained"],
                             sc["n_theoretical"], "",
                             "score %.3f" % sc["discriminating_fraction"]])
            top = res.ranking[0][0] if res.ranking else ""
            _, frs = R.fragments_for(next(s for s in cands if s.name == top),
                                     c["path"], c["depth"], c["adduct"])
            pred = [(f.mz(c["adduct"]), R.label(f, cands[0]), False) for f in frs]
            hot = {p[0] for p in pk}
            pred = [(m, l, any(abs(m - h) * 1000 <= c["tol"] for h in hot))
                    for m, l, _ in pred]
            return hdr, rows, pred, "top candidate: %s" % top
        self.run(job, "interpreting peaks")

    def do_mixture(self):
        c = self.cfg()
        pk = self.peak_list()
        if not pk:
            messagebox.showwarning(APP, "Enter or load an observed peak list first.")
            return

        def job():
            cands = self.candidates(c)
            mm = R.minimal_mixture(pk, cands, c["path"], c["depth"],
                                   c["adduct"], c["tol"], max_components=3)
            hdr = ["solution", "components", "peaks covered"]
            rows = []
            for i, sol in enumerate(mm["solutions"], 1):
                cov = set().union(*(mm["coverage"][n] for n in sol))
                rows.append([i, " + ".join(sol), "%d/%d" % (len(cov), len(pk))])
            n = mm["n_components"]
            note = ("%d component(s) required; %d equally small solution(s)"
                    % (n, len(mm["solutions"]))) if n else \
                   "no combination up to 3 components explains every peak"
            return hdr, rows, [], note
        self.run(job, "searching mixtures")

    # ---- io ------------------------------------------------------------

    def load_peaks(self):
        fn = filedialog.askopenfilename(
            title="Load peak list",
            filetypes=[("Peak lists", "*.csv *.txt *.tsv *.mgf"),
                       ("All files", "*.*")])
        if not fn:
            return
        keep = []
        with open(fn, "r", errors="ignore") as fh:
            for line in fh:
                s = line.replace(",", " ").replace("\t", " ").strip()
                if not s or not (s[0].isdigit() or s[0] == "."):
                    continue
                keep.append(s)
        self.peaks.delete("1.0", "end")
        self.peaks.insert("1.0", "\n".join(keep))
        self.status.config(text="loaded %d peaks from %s"
                           % (len(keep), os.path.basename(fn)))

    def export(self):
        if not self.rows:
            messagebox.showwarning(APP, "Nothing to export yet.")
            return
        fn = filedialog.asksaveasfilename(defaultextension=".csv",
                                          filetypes=[("CSV", "*.csv")])
        if not fn:
            return
        with open(fn, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(self.headers)
            w.writerows(self.rows)
        self.status.config(text="wrote %s" % os.path.basename(fn))


def main():
    App().mainloop()


if __name__ == "__main__":
    sys.exit(main())
