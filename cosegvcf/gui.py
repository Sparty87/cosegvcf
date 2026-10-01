"""Tkinter graphical interface for CoSegVCF."""

import copy
import os
import queue
import threading
import traceback
from collections import OrderedDict

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from . import __version__
from .core import (DEFAULT_AF_FIELDS, MODEL_HELP, MODELS, PAR_REGIONS, ROLE_SEX, ROLES, SEXES,
                   STATUSES, CoSegError, deduce_relationships, load_samples, read_ped,
                   run_analysis, validate, write_ped)


class ScrollFrame:
    """A vertically scrollable frame holding the pedigree table."""

    def __init__(self, parent, height=230):
        self.frame = ttk.Frame(parent)
        self.canvas = tk.Canvas(self.frame, height=height, highlightthickness=0)
        vsb = ttk.Scrollbar(self.frame, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(
            scrollregion=self.canvas.bbox("all")))
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(
            self._win, width=e.width))
        self.canvas.configure(yscrollcommand=vsb.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")


class App:
    HEAD = ["Sample", "File", "Relationship", "Sex", "Father", "Mother", "Status"]
    WIDTHS = [22, 26, 20, 5, 20, 20, 12]

    def __init__(self, root):
        self.root = root
        root.title("CoSegVCF %s - family variant segregation" % __version__)
        root.geometry("1180x920")
        self.files = []
        self.samples = OrderedDict()
        self.rowvars = {}
        self.q = queue.Queue()
        self.worker = None
        self._build()
        root.after(100, self._poll)

    # --- layout -------------------------------------------------------------------
    def _build(self):
        pad = {"padx": 6, "pady": 4}
        main = ttk.Frame(self.root)
        main.pack(fill="both", expand=True, padx=8, pady=6)

        f1 = ttk.LabelFrame(main, text="1. VCF / gVCF files")
        f1.pack(fill="x", **pad)
        self.lb = tk.Listbox(f1, height=4)
        self.lb.pack(side="left", fill="x", expand=True, padx=6, pady=6)
        b1 = ttk.Frame(f1)
        b1.pack(side="right", padx=6)
        ttk.Button(b1, text="Add VCF...", command=self._add_files).pack(fill="x", pady=2)
        ttk.Button(b1, text="Remove selected", command=self._remove_file).pack(fill="x", pady=2)

        f2 = ttk.LabelFrame(main, text="2. Samples, relationships and disease status")
        f2.pack(fill="both", expand=True, **pad)
        self.sf = ScrollFrame(f2)
        self.sf.frame.pack(fill="both", expand=True, padx=6, pady=4)
        b2 = ttk.Frame(f2)
        b2.pack(fill="x", padx=6, pady=4)
        ttk.Button(b2, text="Fill parents and sex from relationships",
                   command=self._deduce).pack(side="left")
        ttk.Button(b2, text="Import PED...", command=self._import_ped).pack(side="left", padx=6)
        ttk.Button(b2, text="Export PED...", command=self._export_ped).pack(side="left")
        ttk.Label(b2, text="Relationships are relative to the proband. Status 'unknown' "
                           "does not constrain filtering.",
                  foreground="gray").pack(side="left", padx=12)
        self._render_rows()

        f3 = ttk.LabelFrame(main, text="3. Inheritance model and filters")
        f3.pack(fill="x", **pad)
        ttk.Label(f3, text="Suspected model:").grid(row=0, column=0, sticky="w", **pad)
        self.model_var = tk.StringVar(value=MODELS["AD"])
        cb = ttk.Combobox(f3, textvariable=self.model_var, values=list(MODELS.values()),
                          state="readonly", width=56)
        cb.grid(row=0, column=1, columnspan=3, sticky="w", **pad)
        cb.bind("<<ComboboxSelected>>", lambda e: self._update_help())
        self.help_lbl = ttk.Label(f3, wraplength=1080, foreground="#1f4e79", justify="left")
        self.help_lbl.grid(row=1, column=0, columnspan=6, sticky="w", **pad)
        self._update_help()

        ttk.Label(f3, text="Genome build (PAR):").grid(row=2, column=0, sticky="w", **pad)
        self.build_var = tk.StringVar(value="GRCh38")
        ttk.Combobox(f3, textvariable=self.build_var, values=list(PAR_REGIONS),
                     state="readonly", width=12).grid(row=2, column=1, sticky="w", **pad)
        ttk.Label(f3, text="Min DP:").grid(row=2, column=2, sticky="e", **pad)
        self.dp_var = tk.StringVar(value="8")
        ttk.Spinbox(f3, from_=0, to=1000, textvariable=self.dp_var, width=6).grid(
            row=2, column=3, sticky="w", **pad)
        ttk.Label(f3, text="Min GQ:").grid(row=2, column=4, sticky="e", **pad)
        self.gq_var = tk.StringVar(value="20")
        ttk.Spinbox(f3, from_=0, to=99, textvariable=self.gq_var, width=6).grid(
            row=2, column=5, sticky="w", **pad)

        self.pass_var = tk.BooleanVar(value=True)
        self.absent_var = tk.BooleanVar(value=True)
        self.miss_var = tk.BooleanVar(value=False)
        self.pen_var = tk.BooleanVar(value=False)
        self.obl_var = tk.BooleanVar(value=True)
        checks = [
            (self.pass_var, "FILTER=PASS calls only"),
            (self.absent_var, "Absent from a VCF without gVCF blocks = 0/0"),
            (self.miss_var, "Missing / low-quality genotypes do not reject"),
            (self.pen_var, "Incomplete penetrance (unaffected may carry)"),
            (self.obl_var, "Check obligate carriers (parents for AR, mothers for XLR)"),
        ]
        for i, (var, txt) in enumerate(checks):
            ttk.Checkbutton(f3, text=txt, variable=var).grid(
                row=3 + i // 3, column=(i % 3) * 2, columnspan=2, sticky="w", **pad)

        ttk.Label(f3, text="Max population AF:").grid(row=5, column=0, sticky="w", **pad)
        self.af_var = tk.StringVar(value="")
        ttk.Entry(f3, textvariable=self.af_var, width=8).grid(row=5, column=1, sticky="w", **pad)
        ttk.Label(f3, text="INFO/CSQ fields:").grid(row=5, column=2, sticky="e", **pad)
        self.af_fields_var = tk.StringVar(value=DEFAULT_AF_FIELDS)
        ttk.Entry(f3, textvariable=self.af_fields_var, width=60).grid(
            row=5, column=3, columnspan=3, sticky="w", **pad)

        ttk.Label(f3, text="Gene BED (optional, compound het):").grid(
            row=6, column=0, columnspan=2, sticky="w", **pad)
        self.bed_var = tk.StringVar()
        ttk.Entry(f3, textvariable=self.bed_var, width=60).grid(
            row=6, column=2, columnspan=3, sticky="we", **pad)
        ttk.Button(f3, text="Browse...", command=self._browse_bed).grid(
            row=6, column=5, sticky="w", **pad)

        f4 = ttk.LabelFrame(main, text="4. Output")
        f4.pack(fill="both", expand=True, **pad)
        r = ttk.Frame(f4)
        r.pack(fill="x", padx=6, pady=4)
        ttk.Label(r, text="Output VCF:").pack(side="left")
        self.out_var = tk.StringVar()
        ttk.Entry(r, textvariable=self.out_var).pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(r, text="Save as...", command=self._browse_out).pack(side="left")
        r2 = ttk.Frame(f4)
        r2.pack(fill="x", padx=6, pady=4)
        self.run_btn = ttk.Button(r2, text="Run", command=self._run)
        self.run_btn.pack(side="left")
        self.pb = ttk.Progressbar(r2, mode="indeterminate", length=260)
        self.pb.pack(side="left", padx=10)
        self.log_txt = tk.Text(f4, height=11, wrap="word", font=("Courier", 9))
        self.log_txt.pack(fill="both", expand=True, padx=6, pady=4)

    def _model_code(self):
        for k, v in MODELS.items():
            if v == self.model_var.get():
                return k
        return "AD"

    def _update_help(self):
        self.help_lbl.configure(text=MODEL_HELP[self._model_code()])

    # --- files ------------------------------------------------------------------------
    def _add_files(self):
        paths = filedialog.askopenfilenames(
            title="Select VCF/gVCF files",
            filetypes=[("VCF", "*.vcf *.vcf.gz *.g.vcf *.g.vcf.gz *.gvcf *.gvcf.gz"),
                       ("All files", "*")])
        self._sync()
        for p in paths:
            if p in self.files:
                continue
            try:
                added = load_samples([p], existing=self.samples.values())
            except (CoSegError, OSError) as e:
                messagebox.showerror("Error", str(e))
                continue
            self.files.append(p)
            self.lb.insert("end", p)
            for s in added:
                self.samples[s.name] = s
            self._log("Added %s (%d samples)" % (os.path.basename(p), len(added)))
        if not self.out_var.get() and self.files:
            self.out_var.set(os.path.join(os.path.dirname(self.files[0]),
                                          "cosegvcf_output.vcf"))
        self._render_rows()

    def _remove_file(self):
        sel = self.lb.curselection()
        if not sel:
            return
        self._sync()
        p = self.files.pop(sel[0])
        self.lb.delete(sel[0])
        gone = [n for n, s in self.samples.items() if s.path == p]
        for n in gone:
            del self.samples[n]
        for s in self.samples.values():
            if s.father in gone:
                s.father = "0"
            if s.mother in gone:
                s.mother = "0"
        self._render_rows()

    # --- pedigree table ---------------------------------------------------------------
    def _render_rows(self):
        for w in self.sf.inner.winfo_children():
            w.destroy()
        self.rowvars = {}
        for c, h in enumerate(self.HEAD):
            ttk.Label(self.sf.inner, text=h, font=("TkDefaultFont", 9, "bold")).grid(
                row=0, column=c, sticky="w", padx=3, pady=2)
        if not self.samples:
            ttk.Label(self.sf.inner, text="Add one or more VCF files to fill this table.",
                      foreground="gray").grid(row=1, column=0, columnspan=7, sticky="w", padx=3)
            return
        names = list(self.samples)
        w = self.WIDTHS
        for r, (n, s) in enumerate(self.samples.items(), start=1):
            v = {k: tk.StringVar(value=getattr(s, k))
                 for k in ("role", "sex", "father", "mother", "status")}
            self.rowvars[n] = v
            others = ["0"] + [x for x in names if x != n]
            ttk.Label(self.sf.inner, text=n, width=w[0]).grid(row=r, column=0, sticky="w", padx=3)
            ttk.Label(self.sf.inner, text=os.path.basename(s.path), width=w[1],
                      foreground="gray").grid(row=r, column=1, sticky="w", padx=3)
            for c, (key, values) in enumerate([("role", ROLES), ("sex", SEXES),
                                                ("father", others), ("mother", others),
                                                ("status", STATUSES)], start=2):
                ttk.Combobox(self.sf.inner, textvariable=v[key], values=values,
                             state="readonly", width=w[c]).grid(row=r, column=c, padx=3, pady=1)
            v["role"].trace_add("write", lambda *a, name=n: self._role_changed(name))

    def _role_changed(self, name):
        v = self.rowvars.get(name)
        if v and v["sex"].get() == "?" and v["role"].get() in ROLE_SEX:
            v["sex"].set(ROLE_SEX[v["role"].get()])

    def _sync(self):
        for n, v in self.rowvars.items():
            s = self.samples.get(n)
            if s:
                for k, var in v.items():
                    setattr(s, k, var.get())

    def _deduce(self):
        self._sync()
        try:
            warns = deduce_relationships(list(self.samples.values()))
        except CoSegError as e:
            messagebox.showwarning("Relationships", str(e))
            return
        self._render_rows()
        for w in warns:
            self._log("WARNING: " + w)
        self._log("Parents and sex filled from relationships (existing values kept).")

    def _import_ped(self):
        p = filedialog.askopenfilename(title="Import PED",
                                       filetypes=[("PED", "*.ped *.fam *.txt"), ("All", "*")])
        if not p:
            return
        self._sync()
        n = read_ped(p, list(self.samples.values()))
        self._render_rows()
        self._log("PED imported: %d samples updated." % n)

    def _export_ped(self):
        self._sync()
        p = filedialog.asksaveasfilename(title="Export PED", defaultextension=".ped",
                                         filetypes=[("PED", "*.ped")])
        if p:
            write_ped(p, list(self.samples.values()))
            self._log("PED exported: " + p)

    def _browse_bed(self):
        p = filedialog.askopenfilename(filetypes=[("BED", "*.bed *.bed.gz"), ("All", "*")])
        if p:
            self.bed_var.set(p)

    def _browse_out(self):
        p = filedialog.asksaveasfilename(title="Output VCF", defaultextension=".vcf",
                                         filetypes=[("VCF", "*.vcf"), ("VCF gzip", "*.vcf.gz")])
        if p:
            self.out_var.set(p)

    # --- run -------------------------------------------------------------------------------
    def _collect_cfg(self):
        try:
            min_dp, min_gq = int(self.dp_var.get()), int(self.gq_var.get())
        except ValueError:
            raise CoSegError("Min DP and Min GQ must be integers.")
        af = self.af_var.get().strip().replace(",", ".")
        max_af = None
        if af:
            try:
                max_af = float(af)
            except ValueError:
                raise CoSegError("Invalid maximum AF (e.g. 0.01).")
        out = self.out_var.get().strip()
        if not out:
            raise CoSegError("Choose an output VCF.")
        bed = self.bed_var.get().strip()
        if bed and not os.path.exists(bed):
            raise CoSegError("BED file not found: " + bed)
        return {"model": self._model_code(), "build": self.build_var.get(),
                "min_dp": min_dp, "min_gq": min_gq, "pass_only": self.pass_var.get(),
                "absent_as_ref": self.absent_var.get(), "ignore_missing": self.miss_var.get(),
                "incomplete_penetrance": self.pen_var.get(),
                "obligate_carriers": self.obl_var.get(), "max_af": max_af,
                "af_fields": self.af_fields_var.get(), "gene_bed": bed, "output": out}

    def _run(self):
        if self.worker and self.worker.is_alive():
            return
        self._sync()
        if not self.samples:
            messagebox.showerror("Error", "No VCF loaded.")
            return
        try:
            cfg = self._collect_cfg()
        except CoSegError as e:
            messagebox.showerror("Error", str(e))
            return
        samples = [copy.copy(s) for s in self.samples.values()]
        errs, warns = validate(samples, cfg)
        if errs:
            messagebox.showerror("Pedigree / settings", "\n".join(errs))
            return
        if warns and not messagebox.askyesno("Warnings", "\n\n".join(warns) + "\n\nContinue?"):
            return
        self.run_btn.state(["disabled"])
        self.pb.start(12)
        self.worker = threading.Thread(target=self._work, args=(samples, cfg), daemon=True)
        self.worker.start()

    def _work(self, samples, cfg):
        try:
            n = run_analysis(samples, cfg, log=lambda m: self.q.put(("log", m)))
            self.q.put(("done", n))
        except Exception as e:  # reported in the GUI, never lost in a thread
            self.q.put(("error", "%s\n%s" % (e, traceback.format_exc())))

    def _poll(self):
        try:
            while True:
                kind, msg = self.q.get_nowait()
                if kind == "log":
                    self._log(msg)
                    continue
                self.pb.stop()
                self.run_btn.state(["!disabled"])
                if kind == "done":
                    messagebox.showinfo("Done", "Compatible variants: %d" % msg)
                else:
                    self._log("ERROR: " + msg)
                    messagebox.showerror("Error", msg.splitlines()[0])
        except queue.Empty:
            pass
        self.root.after(100, self._poll)

    def _log(self, msg):
        self.log_txt.insert("end", msg + "\n")
        self.log_txt.see("end")


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    App(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    main()
