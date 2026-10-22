# Build full-database background counts from all parsed quarters (for ROR/PRR/BCPNN/EBGM denominators)
# Dedup rule: keep highest primaryid per caseid (latest report version; primaryid is sequential)
import csv, json, os, re, sys, zipfile
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
LYMPHOMA = re.compile(r"lymphoma|CLL|SLL|waldenstr|richter|macroglobulin", re.I)

def iter_table(zf, prefix):
    hits = [n for n in zf.namelist()
            if os.path.basename(n).upper().startswith(prefix) and n.upper().endswith(".TXT")]
    if not hits:
        return
    with zf.open(hits[0]) as fh:
        text = (line.decode("latin-1", errors="replace") for line in fh)
        yield from csv.DictReader(text, delimiter="$")

zips = sorted(f for f in os.listdir(BASE) if f.startswith("faers_ascii_") and f.endswith(".zip"))

# Pass 1: dedup map caseid -> max primaryid
best = {}
for f in zips:
    with zipfile.ZipFile(os.path.join(BASE, f)) as zf:
        for row in iter_table(zf, "DEMO"):
            try:
                cid = int(row["caseid"]); pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            if pid > best.get(cid, 0):
                best[cid] = pid
    print(f"[pass1] {f} done, cases so far={len(best)}", flush=True)
keep_pid = set(best.values())
N = len(keep_pid)

# Pass 2: REAC (PT counts), OUTC (fatal), INDI (lymphoma flag) for kept pids
pt_cases = defaultdict(set)      # pt -> set(caseid)
fatal_cases = set()
lymph_cases = set()
for f in zips:
    with zipfile.ZipFile(os.path.join(BASE, f)) as zf:
        for row in iter_table(zf, "REAC"):
            try:
                pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            if pid in keep_pid:
                pt = (row.get("pt") or "").strip().upper()
                if pt:
                    pt_cases[pt].add(pid)  # pid unique per case after dedup
        for row in iter_table(zf, "OUTC"):
            try:
                pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            if pid in keep_pid and (row.get("outc_cod") or "") == "DE":
                fatal_cases.add(pid)
        for row in iter_table(zf, "INDI"):
            try:
                pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            if pid in keep_pid and LYMPHOMA.search(row.get("indi_pt") or ""):
                lymph_cases.add(pid)
    print(f"[pass2] {f} done, pts={len(pt_cases)}, fatal={len(fatal_cases)}, lymph={len(lymph_cases)}", flush=True)

out = {
    "N_total_cases": N,
    "N_fatal": len(fatal_cases),
    "N_lymphoma": len(lymph_cases),
    "pt_counts": {pt: len(s) for pt, s in pt_cases.items()},
    "fatal_pt_counts": {pt: len(s & fatal_cases) for pt, s in pt_cases.items()},
    "lymphoma_pt_counts": {pt: len(s & lymph_cases) for pt, s in pt_cases.items()},
}
with open(os.path.join(BASE, "background_counts.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh)
print(f"=== BACKGROUND DONE: N={N}, fatal={len(fatal_cases)}, lymphoma={len(lymph_cases)}, PTs={len(pt_cases)} ===")
