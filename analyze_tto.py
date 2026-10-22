# TTO (time-to-onset) analysis for infection reports of the three CD20xCD3 bispecifics.
# TTO = event onset (DEMO event_dt, primary axis) - latest drug start_dt (THER, PS row) <= event.
# Robustness axis: FDA receipt date (fda_dt, 100% complete) - start_dt.
# Outputs: tto_results.json + console tables; analysis set = PS cases with >=1 infection PT.
import csv, json, os, statistics
from collections import defaultdict
from datetime import date

import numpy as np
from scipy.stats import weibull_min, mannwhitneyu

BASE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261019)

def parse_date(s):
    s = (s or "").strip()
    if len(s) == 8 and s.isdigit():
        try:
            return date(int(s[:4]), int(s[4:6]), int(s[6:8]))
        except ValueError:
            return None
    return None

def tto_for(starts, anchor):
    """TTO in days vs the latest start <= anchor. Returns (tto, flag)."""
    valid = [s for s in starts if s]
    if not valid or anchor is None:
        return None, "missing"
    prior = [s for s in valid if s <= anchor]
    if prior:
        return (anchor - max(prior)).days, "ok"
    return (anchor - min(valid)).days, "negative"  # event before therapy start

def weibull_fit(x):
    c, loc, scale = weibull_min.fit(x, floc=0)
    return float(c), float(scale)

def weibull_boot_ci(x, n_boot=500):
    if len(x) < 5:
        return (float("nan"), float("nan"))
    ks = []
    x = np.asarray(x, dtype=float)
    for _ in range(n_boot):
        s = rng.choice(x, size=len(x), replace=True)
        if np.all(s == s[0]):
            continue
        try:
            c, _, _ = weibull_min.fit(s, floc=0)
            ks.append(c)
        except Exception:
            continue
    if len(ks) < 50:
        return (float("nan"), float("nan"))
    return (float(np.percentile(ks, 2.5)), float(np.percentile(ks, 97.5)))

def ecdf_at(x_sorted, t):
    import bisect
    return bisect.bisect_right(x_sorted, t) / len(x_sorted) if len(x_sorted) else float("nan")

def summarize(x, label):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n == 0:
        return {"label": label, "n": 0}
    k, lam = weibull_fit(x)
    lo, hi = weibull_boot_ci(x)
    return {
        "label": label, "n": n,
        "median": float(np.median(x)),
        "p25": float(np.percentile(x, 25)), "p75": float(np.percentile(x, 75)),
        "mean": float(np.mean(x)), "min": float(x[0]), "max": float(x[-1]),
        "pct_le_30d": float(100 * ecdf_at(x, 30)),
        "pct_le_90d": float(100 * ecdf_at(x, 90)),
        "weibull_shape": k, "weibull_scale": lam,
        "weibull_shape_ci": [lo, hi],
        "ecdf": {str(t): float(ecdf_at(x, t)) for t in (7, 14, 30, 60, 90, 180, 365)},
    }

if __name__ == "__main__":
    ds = list(csv.DictReader(open(os.path.join(BASE, "analysis_dataset.csv"), encoding="utf-8")))
    ther = list(csv.DictReader(open(os.path.join(BASE, "ther_ps_target.csv"), encoding="utf-8")))

    starts_by_pid = defaultdict(list)
    for r in ther:
        if r["role_cod"] == "PS":
            starts_by_pid[int(r["primaryid"])].append(parse_date(r["start_dt"]))

    drugs = ["glofitamab", "epcoritamab", "mosunetuzumab"]
    results = {}
    detail_rows = []
    for drug in drugs:
        cases = [r for r in ds if drug in r["target_drugs"].split(";")
                 and "PS" in r["role_cods"].split(";")
                 and int(r["n_infection_pt"]) > 0]
        rec = {"n_ps_infection_cases": len(cases)}
        for axis, dfield in (("event_dt", "event_dt"), ("fda_dt", "fda_dt")):
            fatal_x, nonfatal_x, neg, miss = [], [], 0, 0
            covid_x = []
            for r in cases:
                anchor = parse_date(r.get(dfield))
                starts = starts_by_pid.get(int(r["primaryid"]), [])
                tto, flag = tto_for(starts, anchor)
                if flag == "missing":
                    miss += 1
                    continue
                if flag == "negative":
                    neg += 1
                    continue
                if r["fatal"] == "1":
                    fatal_x.append(tto)
                else:
                    nonfatal_x.append(tto)
                if "COVID-19 PNEUMONIA" in (r.get("infection_cats") or "") or True:
                    pass
            rec[axis] = {
                "n_fatal_inf": len(fatal_x), "n_nonfatal_inf": len(nonfatal_x),
                "excluded_negative": neg, "excluded_missing_date": miss,
                "fatal": summarize(fatal_x, "fatal_infection"),
                "nonfatal": summarize(nonfatal_x, "nonfatal_infection"),
                "all_infection": summarize(fatal_x + nonfatal_x, "all_infection"),
            }
            if len(fatal_x) >= 3 and len(nonfatal_x) >= 3:
                try:
                    u, p = mannwhitneyu(fatal_x, nonfatal_x, alternative="two-sided")
                    rec[axis]["mwu_p"] = float(p)
                except Exception:
                    pass
        results[drug] = rec

        # COVID-19 pneumonia TTO (both axes, all fatal+nonfatal)
        covid_rec = {}
        for axis, dfield in (("event_dt", "event_dt"), ("fda_dt", "fda_dt")):
            xs = []
            for r in cases:
                anchor = parse_date(r.get(dfield))
                tto, flag = tto_for(starts_by_pid.get(int(r["primaryid"]), []), anchor)
                if flag == "ok" and "COVID-19 PNEUMONIA" in (r.get("infection_cats") or ""):
                    xs.append(tto)
            covid_rec[axis] = summarize(xs, "covid19_pneumonia")
        rec["covid19_pneumonia_tto"] = covid_rec

        # per-case detail for fatal infection cases (event axis)
        for r in cases:
            anchor = parse_date(r.get("event_dt"))
            tto, flag = tto_for(starts_by_pid.get(int(r["primaryid"]), []), anchor)
            if flag == "ok":
                detail_rows.append({"drug": drug, "primaryid": r["primaryid"],
                                    "fatal": int(r["fatal"]), "tto_days": tto,
                                    "infection_cats": r["infection_cats"]})

    with open(os.path.join(BASE, "tto_results.json"), "w", encoding="utf-8") as fh:
        json.dump(results, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(BASE, "tto_detail.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["drug", "primaryid", "fatal", "tto_days", "infection_cats"])
        w.writeheader(); w.writerows(detail_rows)

    # ---- console report ----
    for drug in drugs:
        r = results[drug]
        print(f"\n===== {drug} (PS, infection cases n={r['n_ps_infection_cases']}) =====", flush=True)
        for axis in ("event_dt", "fda_dt"):
            a = r[axis]
            print(f"  [{axis}] analyzable: fatal_inf={a['n_fatal_inf']}, nonfatal_inf={a['n_nonfatal_inf']}, "
                  f"excluded(neg)={a['excluded_negative']}, excluded(missing date)={a['excluded_missing_date']}", flush=True)
            for grp in ("fatal", "nonfatal", "all_infection"):
                s = a[grp]
                if s["n"] == 0:
                    continue
                ci = s.get("weibull_shape_ci", [float('nan'), float('nan')])
                print(f"    {s['label']:18s} n={s['n']:4d} median={s['median']:7.1f}d IQR=({s['p25']:.0f},{s['p75']:.0f}) "
                      f"<=30d={s['pct_le_30d']:5.1f}% <=90d={s['pct_le_90d']:5.1f}% "
                      f"Weibull k={s['weibull_shape']:.2f} ({ci[0]:.2f}-{ci[1]:.2f}) scale={s['weibull_scale']:.1f}", flush=True)
            if "mwu_p" in a:
                print(f"    fatal vs nonfatal Mann-Whitney p={a['mwu_p']:.4g}", flush=True)
        c = r.get("covid19_pneumonia_tto", {})
        for axis in ("event_dt", "fda_dt"):
            s = c.get(axis, {})
            if s.get("n"):
                print(f"  [COVID-19 pneumonia TTO, {axis}] n={s['n']} median={s['median']:.0f}d "
                      f"IQR=({s['p25']:.0f},{s['p75']:.0f}) k={s['weibull_shape']:.2f}", flush=True)
    print("\n=== TTO ANALYSIS DONE ===", flush=True)
