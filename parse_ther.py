# Parse THER (therapy dates) for target-drug PS rows across all 32 quarterly zips.
# THER is keyed by dsg_drug_seq -> must join with DRUG (primaryid + drug_seq) to identify the drug.
# Output: ther_ps_target.csv  (primaryid, target_drug, drug_seq, start_dt, end_dt, dur, dur_cod)
import csv, os, re, zipfile
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
TARGETS = {
    "glofitamab": re.compile(r"glofitamab|columvi|ro[\s\-]?7082859", re.I),
    "epcoritamab": re.compile(r"epcoritamab|epkinly|gen[\s\-]?3013|duobody[\s\-]?cd3", re.I),
    "mosunetuzumab": re.compile(r"mosunetuzumab|lunsumio|ro[\s\-]?7030816|rg[\s\-]?7828", re.I),
}

def iter_table(zf, prefix):
    hits = [n for n in zf.namelist()
            if os.path.basename(n).upper().startswith(prefix) and n.upper().endswith(".TXT")]
    if not hits:
        return
    with zf.open(hits[0]) as fh:
        text = (line.decode("latin-1", errors="replace") for line in fh)
        yield from csv.DictReader(text, delimiter="$")

zips = sorted(f for f in os.listdir(BASE) if f.startswith("faers_ascii_") and f.endswith(".zip"))

# pass 1: DRUG -> (primaryid, drug_seq) -> (target_drug, role) for target rows
seq_info = {}   # (pid, seq) -> (target_drug, role)
for f in zips:
    with zipfile.ZipFile(os.path.join(BASE, f)) as zf:
        for row in iter_table(zf, "DRUG"):
            name = " ".join(filter(None, [row.get("drugname"), row.get("prod_ai")]))
            for drug, rx in TARGETS.items():
                if rx.search(name):
                    try:
                        pid = int(row["primaryid"]); seq = row.get("drug_seq") or ""
                    except (ValueError, KeyError):
                        continue
                    seq_info[(pid, seq)] = (drug, row.get("role_cod") or "")
                    break
    print(f"[pass1] {f} done, target drug rows={len(seq_info)}", flush=True)

ps_pids = {pid for (pid, seq), (d, role) in seq_info.items() if role == "PS"}
print(f"PS target rows: {sum(1 for v in seq_info.values() if v[1]=='PS')}, unique PS pids={len(ps_pids)}", flush=True)

# pass 2: THER -> therapy dates for those (pid, seq)
out_rows = []
seen = set()
for f in zips:
    with zipfile.ZipFile(os.path.join(BASE, f)) as zf:
        for row in iter_table(zf, "THER"):
            try:
                pid = int(row["primaryid"])
            except (ValueError, KeyError):
                continue
            seq = row.get("dsg_drug_seq") or ""
            key = (pid, seq)
            if key not in seq_info:
                continue
            drug, role = seq_info[key]
            start = (row.get("start_dt") or "").strip()
            end = (row.get("end_dt") or "").strip()
            dedup_key = (pid, seq, start, end)
            if dedup_key in seen:
                continue
            seen.add(dedup_key)
            out_rows.append({
                "primaryid": pid, "target_drug": drug, "role_cod": role,
                "drug_seq": seq, "start_dt": start, "end_dt": end,
                "dur": (row.get("dur") or "").strip(),
                "dur_cod": (row.get("dur_cod") or "").strip(),
            })
    print(f"[pass2] {f} done, ther rows={len(out_rows)}", flush=True)

fn = os.path.join(BASE, "ther_ps_target.csv")
with open(fn, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
    w.writeheader(); w.writerows(out_rows)
n_ps = sum(1 for r in out_rows if r["role_cod"] == "PS")
n_start = sum(1 for r in out_rows if r["role_cod"] == "PS" and len(r["start_dt"]) == 8)
print(f"=== THER PARSE DONE: total rows={len(out_rows)} (PS {n_ps}, PS with valid start_dt {n_start}) -> {fn} ===", flush=True)
