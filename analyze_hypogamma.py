# Hypogammaglobulinemia / neutropenia / CRS co-reporting analysis (PS scope)
# Questions: (a) how often do these markers co-occur with infection reports;
# (b) does marker co-reporting associate with FATAL infection outcomes;
# (c) among marker+ vs marker- cases, infection and fatal-infection rates;
# (d) pathogen-category composition shift in hypogamma+ infection cases.
import csv, glob, json, os
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact

BASE = os.path.dirname(os.path.abspath(__file__))

HYPOGAMMA = {"HYPOGAMMAGLOBULINAEMIA", "BLOOD IMMUNOGLOBULIN G DECREASED",
             "IMMUNOGLOBULINS DECREASED", "BLOOD IMMUNOGLOBULIN M DECREASED",
             "BLOOD IMMUNOGLOBULIN A DECREASED"}
NEUTRO = {"NEUTROPENIA", "FEBRILE NEUTROPENIA", "NEUTROPHIL COUNT DECREASED",
          "NEUTROPENIC SEPSIS"}
CRS = {"CYTOKINE RELEASE SYNDROME", "CYTOKINE STORM"}

def or_ci(a, b, c, d):
    # OR with Haldane correction if zero cells
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    import math
    OR = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return OR, OR * math.exp(-1.96 * se), OR * math.exp(1.96 * se)

ds = list(csv.DictReader(open(os.path.join(BASE, "analysis_dataset.csv"), encoding="utf-8")))
pids_all = {int(r["primaryid"]) for r in ds}
case_pts = defaultdict(set)
for f in glob.glob(os.path.join(BASE, "filtered", "reactions_*.csv")):
    for x in csv.DictReader(open(f, encoding="utf-8")):
        pid = int(x["primaryid"])
        if pid in pids_all and x.get("pt"):
            case_pts[pid].add(x["pt"].upper())

drugs = ["glofitamab", "epcoritamab", "mosunetuzumab"]
out = {}
for drug in drugs:
    cases = [r for r in ds if drug in r["target_drugs"].split(";")
             and "PS" in r["role_cods"].split(";")]
    def flags(r):
        pts = case_pts.get(int(r["primaryid"]), set())
        return {
            "hg": bool(pts & HYPOGAMMA), "neu": bool(pts & NEUTRO),
            "crs": bool(pts & CRS), "inf": int(r["n_infection_pt"]) > 0,
            "fatal": r["fatal"] == "1",
            "fatal_inf": int(r["n_infection_pt"]) > 0 and r["fatal"] == "1",
        }
    rec = {"n_ps": len(cases)}
    F = [(r, flags(r)) for r in cases]

    # (a) marker prevalence overall and among infection cases
    for mk in ("hg", "neu", "crs"):
        rec[f"n_{mk}"] = sum(1 for _, f in F if f[mk])
    rec["n_inf"] = sum(1 for _, f in F if f["inf"])
    rec["n_fatal_inf"] = sum(1 for _, f in F if f["fatal_inf"])
    for mk in ("hg", "neu", "crs"):
        rec[f"{mk}_in_inf"] = sum(1 for _, f in F if f["inf"] and f[mk])
        rec[f"{mk}_in_fatal_inf"] = sum(1 for _, f in F if f["fatal_inf"] and f[mk])
        rec[f"{mk}_in_nonfatal_inf"] = sum(1 for _, f in F if f["inf"] and not f["fatal"] and f[mk])

    # (b) marker -> fatal-infection association among infection cases (2x2)
    #     rows: marker+/-; cols: fatal inf / nonfatal inf
    rec["marker_vs_fatal_inf"] = {}
    for mk in ("hg", "neu", "crs"):
        a = rec[f"{mk}_in_fatal_inf"]
        b = rec[f"{mk}_in_nonfatal_inf"]
        c = rec["n_fatal_inf"] - a
        d = (rec["n_inf"] - rec["n_fatal_inf"]) - b
        if a + b + c + d > 0:
            OR, lo, hi = or_ci(a, b, c, d)
            _, p = fisher_exact([[a, b], [c, d]])
            rec["marker_vs_fatal_inf"][mk] = {
                "a_marker_fatal": a, "b_marker_nonfatal": b,
                "c_nomarker_fatal": c, "d_nomarker_nonfatal": d,
                "OR": round(OR, 2), "OR_lo": round(lo, 2), "OR_hi": round(hi, 2),
                "p_fisher": round(p, 4),
                "pct_marker_in_fatal_inf": round(100 * a / max(rec["n_fatal_inf"], 1), 1),
                "pct_marker_in_nonfatal_inf": round(100 * b / max(rec["n_inf"] - rec["n_fatal_inf"], 1), 1),
            }

    # (c) infection rate and fatal-infection rate by marker status (all PS cases)
    rec["rates_by_marker"] = {}
    for mk in ("hg", "neu", "crs"):
        m_plus = [f for _, f in F if f[mk]]
        m_minus = [f for _, f in F if not f[mk]]
        def rate(xs, key):
            return round(100 * sum(1 for x in xs if x[key]) / max(len(xs), 1), 1)
        rec["rates_by_marker"][mk] = {
            "n_plus": len(m_plus), "n_minus": len(m_minus),
            "inf_rate_plus": rate(m_plus, "inf"), "inf_rate_minus": rate(m_minus, "inf"),
            "fatal_rate_plus": rate(m_plus, "fatal"), "fatal_rate_minus": rate(m_minus, "fatal"),
            "fatal_inf_rate_plus": rate(m_plus, "fatal_inf"),
            "fatal_inf_rate_minus": rate(m_minus, "fatal_inf"),
        }

    # (d) hypogamma+ vs hypogamma- infection cases: pathogen category composition
    comp = {"hg_plus": defaultdict(int), "hg_minus": defaultdict(int)}
    for r, f in F:
        if not f["inf"]:
            continue
        key = "hg_plus" if f["hg"] else "hg_minus"
        for x in (r.get("infection_cats") or "").split(";"):
            if x:
                comp[key][x] += 1
    rec["hypogamma_cat_composition"] = {
        "hg_plus": dict(comp["hg_plus"]), "hg_minus": dict(comp["hg_minus"]),
        "n_hg_plus_inf": sum(1 for _, f in F if f["inf"] and f["hg"]),
        "n_hg_minus_inf": sum(1 for _, f in F if f["inf"] and not f["hg"]),
    }

    # (e) co-report network overlaps (2-way/3-way among infection cases)
    inf_cases = [(r, f) for r, f in F if f["inf"]]
    rec["overlap_inf"] = {
        "inf&hg": sum(1 for _, f in inf_cases if f["hg"]),
        "inf&neu": sum(1 for _, f in inf_cases if f["neu"]),
        "inf&crs": sum(1 for _, f in inf_cases if f["crs"]),
        "inf&hg&neu": sum(1 for _, f in inf_cases if f["hg"] and f["neu"]),
        "inf&hg&crs": sum(1 for _, f in inf_cases if f["hg"] and f["crs"]),
        "inf&neu&crs": sum(1 for _, f in inf_cases if f["neu"] and f["crs"]),
        "inf&hg&neu&crs": sum(1 for _, f in inf_cases if f["hg"] and f["neu"] and f["crs"]),
        "inf_all3_fatal": sum(1 for _, f in inf_cases if f["hg"] and f["neu"] and f["crs"] and f["fatal"]),
    }
    out[drug] = rec

with open(os.path.join(BASE, "hypogamma_coreport.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

for drug in drugs:
    r = out[drug]
    print(f"\n===== {drug} (PS n={r['n_ps']}, infection n={r['n_inf']}, fatal-inf n={r['n_fatal_inf']}) =====", flush=True)
    print(f"  markers: hypogamma {r['n_hg']} ({100*r['n_hg']/r['n_ps']:.1f}%), neutropenia {r['n_neu']} ({100*r['n_neu']/r['n_ps']:.1f}%), CRS {r['n_crs']} ({100*r['n_crs']/r['n_ps']:.1f}%)", flush=True)
    for mk, nm in (("hg", "hypogamma"), ("neu", "neutropenia"), ("crs", "CRS")):
        m = r["marker_vs_fatal_inf"].get(mk)
        if m:
            print(f"  {nm:10s} in fatal-inf: {m['pct_marker_in_fatal_inf']}% ({m['a_marker_fatal']}/{r['n_fatal_inf']}) vs nonfatal-inf {m['pct_marker_in_nonfatal_inf']}% ({m['b_marker_nonfatal']}/{r['n_inf']-r['n_fatal_inf']})  OR={m['OR']} ({m['OR_lo']}-{m['OR_hi']}) p={m['p_fisher']}", flush=True)
    for mk, nm in (("hg", "hypogamma"), ("neu", "neutropenia"), ("crs", "CRS")):
        rr = r["rates_by_marker"][mk]
        print(f"  {nm:10s}+ (n={rr['n_plus']}): inf {rr['inf_rate_plus']}% vs marker- {rr['inf_rate_minus']}% | fatal-inf {rr['fatal_inf_rate_plus']}% vs {rr['fatal_inf_rate_minus']}% | fatal {rr['fatal_rate_plus']}% vs {rr['fatal_rate_minus']}%", flush=True)
    cc = r["hypogamma_cat_composition"]
    print(f"  pathogen cats in hg+ infections (n={cc['n_hg_plus_inf']}): {cc['hg_plus']}", flush=True)
    print(f"  pathogen cats in hg- infections (n={cc['n_hg_minus_inf']}): {cc['hg_minus']}", flush=True)
    print(f"  overlaps (infection cases): {r['overlap_inf']}", flush=True)
print("\n=== HYPOGAMMA COREPORT DONE ===", flush=True)
