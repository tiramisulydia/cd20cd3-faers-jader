# Main disproportionality analysis: ROR / PRR / BCPNN(IC) / EBGM
# for the three CD20xCD3 bispecifics x infection PTs
# 2x2: a=drug+PT, b=drug+other, c=other+PT, d=other+other (case level, deduped)
import csv, glob, json, math, os
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))

def load_pts():
    pt2cat = {}
    with open(os.path.join(BASE, "infection_pts.txt"), encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cat, pt = line.split("|", 1)
            pt2cat[pt.strip().upper()] = cat
    return pt2cat

def ror_ci(a, b, c, d):
    v = (a * d) / (b * c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return v, v * math.exp(-1.96 * se), v * math.exp(1.96 * se)

def prr_chi2(a, b, c, d):
    prr = (a / (a + b)) / (c / (c + d))
    n = a + b + c + d
    chi2 = (a * d - b * c) ** 2 * n / ((a + b) * (c + d) * (a + c) * (b + d))
    return prr, chi2

def bcpnn(a, b, c, d):
    # Bate 1998 / Noren: alpha1=alpha2=1, beta1=beta2=1
    n = a + b + c + d
    ab, ac = a + b, a + c
    ic = math.log2(((a + 1) * (n + 2)) / ((ab + 1) * (ac + 1)))
    var = (1 / math.log(2) ** 2) * (
        (n - a + 1) / ((a + 1) * (n + 2)) +
        (n - ab + 1) / ((ab + 1) * (n + 2)) +
        (n - ac + 1) / ((ac + 1) * (n + 2)))
    return ic, ic - 2 * math.sqrt(var)

def ebgm_fit(pairs):
    # single-component Gamma-Poisson (GPS): ML fit of prior (alpha, beta)
    # pairs: list of (a, E). NB likelihood: a ~ NB(alpha, p=beta/(beta+E))
    from scipy.optimize import minimize
    from scipy.special import gammaln
    import numpy as np
    a_arr = np.array([p[0] for p in pairs], dtype=float)
    e_arr = np.array([p[1] for p in pairs], dtype=float)
    def nll(th):
        alpha, beta = np.exp(th[0]), np.exp(th[1])
        p = beta / (beta + e_arr)
        ll = (gammaln(a_arr + alpha) - gammaln(alpha) - gammaln(a_arr + 1)
              + alpha * np.log(p) + a_arr * np.log(1 - p))
        return -ll.sum()
    res = minimize(nll, x0=[0.0, -2.0], method="Nelder-Mead",
                   options={"maxiter": 5000})
    return float(np.exp(res.x[0])), float(np.exp(res.x[1]))

def ebgm_post(a, e, alpha, beta):
    from scipy.stats import gamma as gdist
    return (a + alpha) / (e + beta), gdist.ppf(0.05, a + alpha, scale=1 / (e + beta))

if __name__ == "__main__":
    pt2cat = load_pts()
    ds = list(csv.DictReader(open(os.path.join(BASE, "analysis_dataset.csv"), encoding="utf-8")))
    bg = json.load(open(os.path.join(BASE, "background_counts.json"), encoding="utf-8"))
    N = bg["N_total_cases"]
    pt_counts = bg["pt_counts"]

    # per-case PT sets for target cases
    pid2case = {r["primaryid"]: r for r in ds}
    case_pts = defaultdict(set)
    for f in glob.glob(os.path.join(BASE, "filtered", "reactions_*.csv")):
        for x in csv.DictReader(open(f, encoding="utf-8")):
            if x["primaryid"] in pid2case and x.get("pt"):
                case_pts[x["primaryid"]].add(x["pt"].upper())

    drugs = ["glofitamab", "epcoritamab", "mosunetuzumab"]
    outdir = os.path.join(BASE, "results")
    os.makedirs(outdir, exist_ok=True)

    for scope in ("all_roles", "ps_only"):
        for drug in drugs:
            if scope == "all_roles":
                cases = [r for r in ds if drug in r["target_drugs"].split(";")]
            else:
                cases = [r for r in ds if drug in r["target_drugs"].split(";")
                         and "PS" in r["role_cods"].split(";")]
            nD = len(cases)
            pids = [r["primaryid"] for r in cases]
            pt_a = defaultdict(int)
            for pid in pids:
                for pt in case_pts.get(pid, ()):
                    pt_a[pt] += 1
            # EBGM prior fit on this drug's pairs (a>=1)
            pairs = [(a, nD * pt_counts.get(pt, 0) / N) for pt, a in pt_a.items()
                     if pt_counts.get(pt, 0) > 0]
            alpha, beta = ebgm_fit(pairs)
            rows_out = []
            for pt, a in sorted(pt_a.items()):
                if pt not in pt_counts or a < 1:
                    continue
                c = pt_counts[pt] - a
                b = nD - a
                d = N - a - b - c
                if b <= 0 or c <= 0:
                    continue
                E = nD * pt_counts[pt] / N
                ror, ror_lo, ror_hi = ror_ci(a, b, c, d)
                prr, chi2 = prr_chi2(a, b, c, d)
                ic, ic025 = bcpnn(a, b, c, d)
                ebgm, ebgm05 = ebgm_post(a, E, alpha, beta)
                rows_out.append({
                    "drug": drug, "scope": scope, "pt": pt,
                    "infection_cat": pt2cat.get(pt, ""),
                    "a": a, "N_drug": nD, "N_pt_bg": pt_counts[pt],
                    "ROR": round(ror, 2), "ROR_lo": round(ror_lo, 2), "ROR_hi": round(ror_hi, 2),
                    "PRR": round(prr, 2), "chi2": round(chi2, 2),
                    "IC": round(ic, 2), "IC025": round(ic025, 2),
                    "EBGM": round(ebgm, 2), "EBGM05": round(ebgm05, 2),
                    "signal_ror": int(ror_lo > 1 and a >= 3),
                    "signal_prr": int(prr >= 2 and chi2 >= 4 and a >= 3),
                    "signal_bcpnn": int(ic025 > 0),
                    "signal_ebgm": int(ebgm05 >= 2),
                })
            fn = os.path.join(outdir, f"signals_{drug}_{scope}.csv")
            with open(fn, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
                w.writeheader(); w.writerows(rows_out)
            n_sig = sum(r["signal_ror"] for r in rows_out)
            n_inf = sum(1 for r in rows_out if r["infection_cat"])
            n_inf_sig = sum(1 for r in rows_out if r["infection_cat"] and r["signal_ror"] and r["signal_prr"] and r["signal_bcpnn"] and r["signal_ebgm"])
            print(f"{drug} [{scope}]: cases={nD}, PTs={len(rows_out)}, "
                  f"ROR-signals={n_sig}, infection PTs={n_inf}, 4-method infection signals={n_inf_sig}, "
                  f"EBGM prior a={alpha:.3f} b={beta:.3f}", flush=True)
    print("=== MAIN ANALYSIS DONE ===")
