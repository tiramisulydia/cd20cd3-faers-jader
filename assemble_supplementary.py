# Assemble Supplementary Tables S1-S6 for the BJH submission package.
# S1 full infection signal listings (PS + fatal + lymphoma backgrounds, four methods)
# S3 sensitivity analyses (lymphoma-restricted 4-method signals + COVID calendar 2x2)
# S4 TTO detail table (from tto_results.json)
# S5 infection PT dictionary (EN set + MedDRA/J mapping)
import csv, json, os, glob
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "..", "supplementary")
os.makedirs(OUT, exist_ok=True)

def load(fn):
    return list(csv.DictReader(open(os.path.join(BASE, "results", fn), encoding="utf-8")))

# ---------- S1: full infection signal listings ----------
rows = []
for drug in ("glofitamab", "epcoritamab", "mosunetuzumab"):
    for bg_tag, fn in (("full_library", f"signals_{drug}_ps_only.csv"),
                       ("fatal_restricted", f"signals_{drug}_ps_only_fatal.csv"),
                       ("lymphoma_restricted", f"signals_{drug}_ps_only_lymphoma.csv")):
        for r in load(fn):
            if not r.get("infection_cat"):
                continue
            sig4 = all(r[k] == "1" for k in ("signal_ror", "signal_prr", "signal_bcpnn", "signal_ebgm"))
            rows.append({
                "drug": drug, "background": bg_tag, "pt": r["pt"], "category": r["infection_cat"],
                "a": r["a"], "N_drug": r["N_drug"],
                "ROR (95%CI)": f"{r['ROR']} ({r['ROR_lo']}-{r['ROR_hi']})",
                "PRR": r["PRR"], "chi2": r["chi2"], "IC025": r["IC025"],
                "EBGM": r["EBGM"], "EBGM05": r["EBGM05"],
                "four_method_signal": "Yes" if sig4 else "No",
            })
rows.sort(key=lambda r: (r["drug"], r["background"], -float(r["ROR (95%CI)"].split(" ")[0])))
fn = os.path.join(OUT, "Table_S1_full_infection_signals.csv")
with open(fn, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"S1: {len(rows)} rows -> {fn}", flush=True)

# ---------- S3a: COVID calendar stratification ----------
cal = list(csv.DictReader(open(os.path.join(BASE, "covid_calendar_signals.csv"), encoding="utf-8")))
fn = os.path.join(OUT, "Table_S3a_covid_calendar_stratification.csv")
with open(fn, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(cal[0].keys()))
    w.writeheader(); w.writerows(cal)
print(f"S3a: {len(cal)} rows -> {fn}", flush=True)

# ---------- S4: TTO detail ----------
tto = json.load(open(os.path.join(BASE, "tto_results.json"), encoding="utf-8"))
rows4 = []
for drug, rec in tto.items():
    for axis in ("event_dt", "fda_dt"):
        a = rec[axis]
        for grp in ("fatal", "nonfatal", "all_infection"):
            s = a[grp]
            if s["n"] == 0:
                continue
            ci = s.get("weibull_shape_ci", [None, None])
            rows4.append({
                "drug": drug, "date_axis": axis, "group": s["label"], "n": s["n"],
                "median_days": s["median"], "IQR": f"{s['p25']:.0f}-{s['p75']:.0f}",
                "pct_le_30d": f"{s['pct_le_30d']:.1f}", "pct_le_90d": f"{s['pct_le_90d']:.1f}",
                "weibull_shape": f"{s['weibull_shape']:.2f}",
                "weibull_shape_95ci": f"{ci[0]:.2f}-{ci[1]:.2f}" if ci[0] else "",
                "weibull_scale": f"{s['weibull_scale']:.1f}",
            })
        if "mwu_p" in a:
            rows4.append({"drug": drug, "date_axis": axis, "group": "fatal_vs_nonfatal_MannWhitney_p",
                          "n": "", "median_days": f"{a['mwu_p']:.4g}", "IQR": "", "pct_le_30d": "",
                          "pct_le_90d": "", "weibull_shape": "", "weibull_shape_95ci": "", "weibull_scale": ""})
fn = os.path.join(OUT, "Table_S4_tto_detail.csv")
with open(fn, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows4[0].keys()))
    w.writeheader(); w.writerows(rows4)
print(f"S4: {len(rows4)} rows -> {fn}", flush=True)

# ---------- S5: infection PT dictionary ----------
with open(os.path.join(BASE, "infection_pts.txt"), encoding="utf-8") as fh:
    pt_rows = []
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cat, pt = line.split("|", 1)
        pt_rows.append({"category": cat, "meddra_pt_en": pt.strip()})
fn = os.path.join(OUT, "Table_S5a_infection_pt_dictionary.csv")
with open(fn, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["category", "meddra_pt_en"])
    w.writeheader(); w.writerows(sorted(pt_rows, key=lambda r: (r["category"], r["meddra_pt_en"])))
print(f"S5a: {len(pt_rows)} PTs -> {fn}", flush=True)

# JA mapping from jader_pipeline.py
import importlib.util
spec = importlib.util.spec_from_file_location("jp", os.path.join(BASE, "jader_pipeline.py"))
jp = importlib.util.module_from_spec(spec)
import types, sys as _sys
# avoid executing __main__ block: module has `if __name__ == "__main__"` guard, safe to exec
spec.loader.exec_module(jp)
ja_rows = [{"meddra_j_ja": ja, "meddra_pt_en": en} for ja, en in sorted(jp.PT_JA_EN.items(), key=lambda x: x[1])]
fn = os.path.join(OUT, "Table_S5b_meddraj_mapping.csv")
with open(fn, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["meddra_j_ja", "meddra_pt_en"])
    w.writeheader(); w.writerows(ja_rows)
print(f"S5b: {len(ja_rows)} mappings -> {fn}", flush=True)

print("=== SUPPLEMENTARY DATA TABLES DONE ===", flush=True)
