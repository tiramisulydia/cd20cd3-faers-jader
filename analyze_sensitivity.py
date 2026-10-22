# Sensitivity analyses for the CD20xCD3 bispecific fatal-infection FAERS study
# (2) Disease-restricted control (lymphoma background)  -- replicates Zhong 2026 approach
# (1) Fatal-restricted analysis (fatal background)       -- primary-outcome-adjacent
# 2x2: a=drug+PT, b=drug+other, c=bg+PT, d=bg+other (case level, deduped)
# Backgrounds reused from background_counts.json: pt_counts / lymphoma_pt_counts / fatal_pt_counts
import csv, glob, json, os
from collections import defaultdict
import analyze_main as M

BASE = os.path.dirname(os.path.abspath(__file__))

def run_scope(cases, pt_a_fn, nD, pt_bg, N_bg, tag):
    """Four-method signals for one drug-scope-background combination."""
    pt_a = pt_a_fn()
    pairs = [(a, nD * pt_bg.get(pt, 0) / N_bg) for pt, a in pt_a.items()
             if pt_bg.get(pt, 0) > 0]
    if not pairs:
        return [], None, None
    alpha, beta = M.ebgm_fit(pairs)
    rows = []
    for pt, a in sorted(pt_a.items()):
        c_bg = pt_bg.get(pt, 0)
        if c_bg < 1 or a < 1:
            continue
        c = c_bg - a
        b = nD - a
        d = N_bg - a - b - c
        if b <= 0 or c <= 0 or d <= 0:
            continue
        E = nD * c_bg / N_bg
        ror, ror_lo, ror_hi = M.ror_ci(a, b, c, d)
        prr, chi2 = M.prr_chi2(a, b, c, d)
        ic, ic025 = M.bcpnn(a, b, c, d)
        ebgm, ebgm05 = M.ebgm_post(a, E, alpha, beta)
        rows.append({
            "background": tag, "pt": pt, "infection_cat": M.load_pts().get(pt, ""),
            "a": a, "N_drug": nD, "N_bg": N_bg, "pt_bg": c_bg,
            "ROR": round(ror, 2), "ROR_lo": round(ror_lo, 2), "ROR_hi": round(ror_hi, 2),
            "PRR": round(prr, 2), "chi2": round(chi2, 2),
            "IC": round(ic, 2), "IC025": round(ic025, 2),
            "EBGM": round(ebgm, 2), "EBGM05": round(ebgm05, 2),
            "signal_ror": int(ror_lo > 1 and a >= 3),
            "signal_prr": int(prr >= 2 and chi2 >= 4 and a >= 3),
            "signal_bcpnn": int(ic025 > 0),
            "signal_ebgm": int(ebgm05 >= 2),
        })
    return rows, alpha, beta

if __name__ == "__main__":
    ds = list(csv.DictReader(open(os.path.join(BASE, "analysis_dataset.csv"), encoding="utf-8")))
    bg = json.load(open(os.path.join(BASE, "background_counts.json"), encoding="utf-8"))
    pid2case = {r["primaryid"]: r for r in ds}
    case_pts = defaultdict(set)
    for f in glob.glob(os.path.join(BASE, "filtered", "reactions_*.csv")):
        for x in csv.DictReader(open(f, encoding="utf-8")):
            if x["primaryid"] in pid2case and x.get("pt"):
                case_pts[x["primaryid"]].add(x["pt"].upper())

    drugs = ["glofitamab", "epcoritamab", "mosunetuzumab"]
    outdir = os.path.join(BASE, "results")
    os.makedirs(outdir, exist_ok=True)

    for scope in ("ps_only", "all_roles"):
        for drug in drugs:
            if scope == "all_roles":
                cases = [r for r in ds if drug in r["target_drugs"].split(";")]
            else:
                cases = [r for r in ds if drug in r["target_drugs"].split(";")
                         and "PS" in r["role_cods"].split(";")]
            nD = len(cases)
            pids = [r["primaryid"] for r in cases]
            fatal_cases = [r for r in cases if r["fatal"] == "1"]
            fatal_pids = [r["primaryid"] for r in fatal_cases]

            def pt_a_all(pids=pids):
                cnt = defaultdict(int)
                for pid in pids:
                    for pt in case_pts.get(pid, ()):
                        cnt[pt] += 1
                return cnt

            # ---- (2) lymphoma-restricted ----
            rows, alpha, beta = run_scope(
                cases, lambda: pt_a_all(pids), nD,
                bg["lymphoma_pt_counts"], bg["N_lymphoma"], "lymphoma")
            fn = os.path.join(outdir, f"signals_{drug}_{scope}_lymphoma.csv")
            with open(fn, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader(); w.writerows(rows)
            n4 = sum(1 for r in rows if r["infection_cat"] and all(
                r[k] for k in ("signal_ror", "signal_prr", "signal_bcpnn", "signal_ebgm")))
            print(f"{drug} [{scope}] LYMPHOMA: nD={nD}, N_bg={bg['N_lymphoma']}, "
                  f"4-method infection signals={n4}, prior a={alpha:.3f} b={beta:.3f}", flush=True)

            # ---- (1) fatal-restricted ----
            rows_f, alpha_f, beta_f = run_scope(
                fatal_cases, lambda: pt_a_all(fatal_pids), len(fatal_cases),
                bg["fatal_pt_counts"], bg["N_fatal"], "fatal")
            if rows_f:
                fn = os.path.join(outdir, f"signals_{drug}_{scope}_fatal.csv")
                with open(fn, "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(rows_f[0].keys()))
                    w.writeheader(); w.writerows(rows_f)
                n4f = sum(1 for r in rows_f if r["infection_cat"] and all(
                    r[k] for k in ("signal_ror", "signal_prr", "signal_bcpnn", "signal_ebgm")))
                print(f"{drug} [{scope}] FATAL: nD={len(fatal_cases)}, N_bg={bg['N_fatal']}, "
                      f"4-method infection signals={n4f}, prior a={alpha_f:.3f} b={beta_f:.3f}", flush=True)
            else:
                print(f"{drug} [{scope}] FATAL: no pairs", flush=True)

    # ---- fatal infection category composition (case level, PS scope) ----
    print("\n=== FATAL INFECTION CATEGORY COMPOSITION (PS scope) ===", flush=True)
    comp = {}
    for drug in drugs:
        cases = [r for r in ds if drug in r["target_drugs"].split(";")
                 and "PS" in r["role_cods"].split(";")]
        fatal = [r for r in cases if r["fatal"] == "1"]
        nonfatal = [r for r in cases if r["fatal"] != "1"]
        def catcount(rs):
            c = defaultdict(int)
            for r in rs:
                cats = [x for x in (r.get("infection_cats") or "").split(";") if x]
                for x in cats:
                    c[x] += 1
            return c
        cf, cn = catcount(fatal), catcount(nonfatal)
        n_fatal_inf = sum(1 for r in fatal if r.get("has_fatal_infection") == "1")
        n_nonfatal_inf = sum(1 for r in nonfatal if r.get("has_fatal_infection") == "1")
        comp[drug] = {
            "n_cases": len(cases), "n_fatal": len(fatal), "n_nonfatal": len(nonfatal),
            "fatal_with_infection": n_fatal_inf,
            "nonfatal_with_infection": n_nonfatal_inf,
            "fatal_cats": dict(cf), "nonfatal_cats": dict(cn),
        }
        print(f"{drug}: fatal={len(fatal)} (with infection PT {n_fatal_inf}, "
              f"{100*n_fatal_inf/max(len(fatal),1):.1f}%), "
              f"nonfatal={len(nonfatal)} (with infection PT {n_nonfatal_inf}, "
              f"{100*n_nonfatal_inf/max(len(nonfatal),1):.1f}%)", flush=True)
        print(f"   fatal cats: {dict(cf)}", flush=True)
        print(f"   nonfatal cats: {dict(cn)}", flush=True)
    with open(os.path.join(BASE, "fatal_infection_composition.json"), "w", encoding="utf-8") as fh:
        json.dump(comp, fh, ensure_ascii=False, indent=1)
    print("=== SENSITIVITY + FATAL ANALYSIS DONE ===", flush=True)
