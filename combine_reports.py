# Combine parsed quarters -> deduplicated analysis dataset
# Rules: group by caseid, keep row with max caseversion (latest version of the case)
import csv, glob, os, re, sys
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
FILT = os.path.join(BASE, "filtered")

def load_pt_set():
    pt2cat = {}
    with open(os.path.join(BASE, "infection_pts.txt"), encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cat, pt = line.split("|", 1)
            pt2cat[pt.strip().upper()] = cat
    return pt2cat

def read_all(name):
    rows = []
    for f in sorted(glob.glob(os.path.join(FILT, f"{name}_*.csv"))):
        with open(f, newline="", encoding="utf-8") as fh:
            rows.extend(csv.DictReader(fh))
    return rows

if __name__ == "__main__":
    pt2cat = load_pt_set()
    reports = read_all("reports")
    reactions = read_all("reactions")
    outcomes = read_all("outcomes")
    drugrows = read_all("drugrows")

    # dedup reports by caseid -> keep max caseversion
    best = {}
    for r in reports:
        try:
            cv = int(r.get("caseversion") or 0)
        except ValueError:
            cv = 0
        key = r.get("caseid") or r["primaryid"]
        if key not in best or cv > best[key][0]:
            best[key] = (cv, r)
    keep_pid = {r["primaryid"] for _, r in best.values()}
    dup_removed = len(reports) - len(keep_pid)

    outc = defaultdict(set)
    for o in outcomes:
        if o["primaryid"] in keep_pid and o.get("outc_cod"):
            outc[o["primaryid"]].add(o["outc_cod"])

    reac = defaultdict(list)
    for x in reactions:
        if x["primaryid"] in keep_pid and x.get("pt"):
            reac[x["primaryid"]].append(x["pt"].upper())

    rep = {r["primaryid"]: r for _, r in best.values()}

    # case-level: one row per unique case; target drugs as a set (handles switch/combination)
    case_drugs = defaultdict(set)
    case_roles = defaultdict(set)
    for dr in drugrows:
        if dr["primaryid"] in keep_pid:
            case_drugs[rep[dr["primaryid"]].get("caseid") or dr["primaryid"]].add(dr["target_drug"])
            case_roles[rep[dr["primaryid"]].get("caseid") or dr["primaryid"]].add(dr.get("role_cod") or "")

    out_path = os.path.join(BASE, "analysis_dataset.csv")
    n_written = 0
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["primaryid", "caseid", "quarter", "target_drugs", "role_cods",
                    "fda_dt", "event_dt", "age", "age_cod", "age_grp", "sex",
                    "reporter_country", "occr_country", "report_type",
                    "fatal", "hosp", "life_threat", "n_pt",
                    "n_infection_pt", "infection_cats", "has_fatal_infection"])
        for pid, r in rep.items():
            cid = r.get("caseid") or pid
            if cid not in case_drugs:
                continue
            pts = reac.get(pid, [])
            inf_pts = [p for p in pts if p in pt2cat]
            cats = sorted({pt2cat[p] for p in inf_pts})
            ocs = outc.get(pid, set())
            fatal = "DE" in ocs
            w.writerow([pid, cid, "", ";".join(sorted(case_drugs[cid])),
                        ";".join(sorted(case_roles[cid])), r.get("fda_dt"),
                        r.get("event_dt"), r.get("age"), r.get("age_cod"),
                        r.get("age_grp"), r.get("sex"), r.get("reporter_country"),
                        r.get("occr_country"), r.get("report_type"),
                        int(fatal), int("HO" in ocs), int("LT" in ocs),
                        len(set(pts)), len(set(inf_pts)), ";".join(cats),
                        int(fatal and bool(inf_pts))])
            n_written += 1
    print(f"reports(raw)={len(reports)} dedup_removed={dup_removed} unique_cases={n_written}")
    print(f"dataset -> {out_path}")
