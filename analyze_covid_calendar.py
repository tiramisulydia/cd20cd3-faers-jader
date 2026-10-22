# COVID-19 calendar-stratified analysis (2020-2022 pandemic vs 2023+ post-pandemic)
# Needs period-restricted background counts -> targeted re-scan of the 32 quarterly zips:
#   pass1 DEMO: dedup map (caseid -> max primaryid) + event year of kept version
#   pass2 REAC: COVID PT counts by period; OUTC: fatal flag by period
# Drug-side counts come from analysis_dataset.csv (event_dt) + filtered/reactions_*.csv
import csv, glob, json, os, zipfile
from collections import defaultdict
import analyze_main as M

BASE = os.path.dirname(os.path.abspath(__file__))
COVID_PTS = {"COVID-19", "COVID-19 PNEUMONIA"}
EXTRA_PTS = {"RESPIRATORY SYNCYTIAL VIRUS INFECTION", "INFLUENZA"}

def iter_table(zf, prefix):
    hits = [n for n in zf.namelist()
            if os.path.basename(n).upper().startswith(prefix) and n.upper().endswith(".TXT")]
    if not hits:
        return
    with zf.open(hits[0]) as fh:
        text = (line.decode("latin-1", errors="replace") for line in fh)
        yield from csv.DictReader(text, delimiter="$")

def event_year(s):
    s = (s or "").strip()
    if len(s) >= 4 and s[:4].isdigit():
        y = int(s[:4])
        if 2015 <= y <= 2030:
            return y
    return None

def period_of(y):
    if y is None:
        return "unknown"
    if y <= 2019:
        return "pre2019"
    if y <= 2022:
        return "pandemic2020_2022"
    return "post2023"

zips = sorted(f for f in os.listdir(BASE) if f.startswith("faers_ascii_") and f.endswith(".zip"))

# ---- pass 1: dedup + event year (both event_dt onset and fda_dt receipt) ----
best = {}          # caseid -> (primaryid, year_event, year_fda)
for f in zips:
    with zipfile.ZipFile(os.path.join(BASE, f)) as zf:
        for row in iter_table(zf, "DEMO"):
            try:
                cid = int(row["caseid"]); pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            if pid > best.get(cid, (0,))[0]:
                best[cid] = (pid, event_year(row.get("event_dt")),
                             event_year(row.get("fda_dt")))
    print(f"[pass1] {f} done, cases={len(best)}", flush=True)

keep = {pid: (ye, yf) for pid, ye, yf in best.values()}
N = len(keep)
period_N = defaultdict(int)
period_N_fda = defaultdict(int)
for ye, yf in keep.values():
    period_N[period_of(ye)] += 1
    period_N_fda[period_of(yf)] += 1

# ---- pass 2: COVID/extra PT counts + fatal, by period (event_dt and fda_dt axes) ----
pt_period = defaultdict(lambda: defaultdict(int))   # pt -> period -> n cases (event_dt)
pt_period_fatal = defaultdict(lambda: defaultdict(int))
pt_period_fda = defaultdict(lambda: defaultdict(int))  # pt -> period -> n (fda_dt)
for f in zips:
    with zipfile.ZipFile(os.path.join(BASE, f)) as zf:
        fatal_pids = set()
        for row in iter_table(zf, "OUTC"):
            try:
                pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            if pid in keep and (row.get("outc_cod") or "") == "DE":
                fatal_pids.add(pid)
        for row in iter_table(zf, "REAC"):
            try:
                pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            if pid not in keep:
                continue
            pt = (row.get("pt") or "").strip().upper()
            if pt in COVID_PTS or pt in EXTRA_PTS:
                ye, yf = keep[pid]
                pe, pf = period_of(ye), period_of(yf)
                pt_period[pt][pe] += 1
                pt_period_fda[pt][pf] += 1
                if pid in fatal_pids:
                    pt_period_fatal[pt][pe] += 1
    print(f"[pass2] {f} done", flush=True)

bg = {
    "N_total": N,
    "period_N": dict(period_N),
    "period_N_fda": dict(period_N_fda),
    "pt_period": {pt: dict(v) for pt, v in pt_period.items()},
    "pt_period_fda": {pt: dict(v) for pt, v in pt_period_fda.items()},
    "pt_period_fatal": {pt: dict(v) for pt, v in pt_period_fatal.items()},
}
with open(os.path.join(BASE, "covid_calendar_background.json"), "w", encoding="utf-8") as fh:
    json.dump(bg, fh, ensure_ascii=False, indent=1)
print(f"=== COVID CALENDAR BACKGROUND DONE: N={N} ===", flush=True)
print(f"  event_dt periods: {dict(period_N)}", flush=True)
print(f"  fda_dt  periods: {dict(period_N_fda)}", flush=True)
for pt in sorted(pt_period):
    print(f"  {pt}: event {dict(pt_period[pt])} | fda {dict(pt_period_fda[pt])} (fatal: {dict(pt_period_fatal[pt])})", flush=True)

# ---- drug-side: COVID PT counts by period, PS scope (event_dt + fda_dt axes) ----
ds = list(csv.DictReader(open(os.path.join(BASE, "analysis_dataset.csv"), encoding="utf-8")))
pid_year = {}
pid_year_fda = {}
for r in ds:
    pid_year[r["primaryid"]] = event_year(r.get("event_dt"))
    pid_year_fda[r["primaryid"]] = event_year(r.get("fda_dt"))

drugs = ["glofitamab", "epcoritamab", "mosunetuzumab"]
rows = []
for drug in drugs:
    cases = [r for r in ds if drug in r["target_drugs"].split(";")
             and "PS" in r["role_cods"].split(";")]
    # PT sets per case from filtered reactions
    pids = {r["primaryid"] for r in cases}
    case_pts = defaultdict(set)
    for f in glob.glob(os.path.join(BASE, "filtered", "reactions_*.csv")):
        for x in csv.DictReader(open(f, encoding="utf-8")):
            if x["primaryid"] in pids and x.get("pt"):
                case_pts[x["primaryid"]].add(x["pt"].upper())
    period_cases = defaultdict(int)
    period_cases_fda = defaultdict(int)
    period_pt = defaultdict(lambda: defaultdict(int))
    period_pt_fda = defaultdict(lambda: defaultdict(int))
    period_fatal = defaultdict(int)
    for r in cases:
        pe = period_of(pid_year.get(r["primaryid"]))
        pf = period_of(pid_year_fda.get(r["primaryid"]))
        period_cases[pe] += 1
        period_cases_fda[pf] += 1
        if r["fatal"] == "1":
            period_fatal[pe] += 1
        for pt in case_pts.get(r["primaryid"], ()):
            if pt in COVID_PTS or pt in EXTRA_PTS:
                period_pt[pt][pe] += 1
                period_pt_fda[pt][pf] += 1
    print(f"\n{drug} PS: event-axis cases={dict(period_cases)}, fda-axis cases={dict(period_cases_fda)}, fatal={dict(period_fatal)}", flush=True)
    for pt in sorted(period_pt):
        print(f"  {pt}: event {dict(period_pt[pt])} | fda {dict(period_pt_fda[pt])}", flush=True)
    # 2x2 per period per PT vs period-restricted background, both axes
    for axis, pcases, ppts, bg_pt, bg_N in (
            ("event_dt", period_cases, period_pt, bg["pt_period"], bg["period_N"]),
            ("fda_dt", period_cases_fda, period_pt_fda, bg["pt_period_fda"], bg["period_N_fda"])):
        for pt in sorted(COVID_PTS | EXTRA_PTS):
            for p in ("pandemic2020_2022", "post2023"):
                a = ppts.get(pt, {}).get(p, 0)
                nD = pcases.get(p, 0)
                c_bg = bg_pt.get(pt, {}).get(p, 0)
                N_bg = bg_N.get(p, 0)
                if a == 0 or c_bg == 0 or nD == 0:
                    continue
                b = nD - a
                c = c_bg - a
                d = N_bg - a - b - c
                if b <= 0 or c <= 0 or d <= 0:
                    continue
                ror, lo, hi = M.ror_ci(a, b, c, d)
                prr, chi2 = M.prr_chi2(a, b, c, d)
                rows.append({"drug": drug, "pt": pt, "period": p, "axis": axis,
                             "a": a, "nD_period": nD, "c_bg": c_bg, "N_bg": N_bg,
                             "ROR": round(ror, 2), "ROR_lo": round(lo, 2), "ROR_hi": round(hi, 2),
                             "PRR": round(prr, 2), "chi2": round(chi2, 2)})
with open(os.path.join(BASE, "covid_calendar_signals.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("\n=== COVID CALENDAR SIGNALS ===", flush=True)
for r in rows:
    print(f"{r['drug']:14s} {r['pt'][:38]:40s} {r['period']:18s} [{r['axis']:8s}] a={r['a']:3d}/{r['nD_period']:4d} ROR={r['ROR']:8.2f} ({r['ROR_lo']}-{r['ROR_hi']})", flush=True)
