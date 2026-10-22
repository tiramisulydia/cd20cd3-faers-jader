# FAERS quarterly ASCII parser: extract target-drug reports (report-level + reaction-level)
# Streams directly from zip files; no full extraction needed.
import csv, glob, os, re, sys, zipfile

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "filtered")
os.makedirs(OUT, exist_ok=True)

TARGETS = {
    # 通用名 + 商品名 + 研发代号（试验期报告用代号：glofitamab=RO7082859、mosunetuzumab=RO7030816/RG7828、epcoritamab=GEN3013）
    "glofitamab": re.compile(r"glofitamab|columvi|ro[\s\-]?7082859", re.I),
    "epcoritamab": re.compile(r"epcoritamab|epkinly|gen[\s\-]?3013|duobody[\s\-]?cd3", re.I),
    "mosunetuzumab": re.compile(r"mosunetuzumab|lunsumio|ro[\s\-]?7030816|rg[\s\-]?7828", re.I),
}

def iter_table(zf, name_prefix):
    """Yield dict rows from ascii/<PREFIX>YYQq.txt inside zip; header auto-detected."""
    hits = [n for n in zf.namelist()
            if os.path.basename(n).upper().startswith(name_prefix) and n.upper().endswith(".TXT")]
    if not hits:
        return
    with zf.open(hits[0]) as fh:
        text = (line.decode("latin-1", errors="replace") for line in fh)
        reader = csv.DictReader(text, delimiter="$")
        for row in reader:
            yield row

def parse_quarter(zpath):
    qtag = os.path.basename(zpath).replace("faers_ascii_", "").replace(".zip", "")
    with zipfile.ZipFile(zpath) as zf:
        # 1) DRUG: find target-drug rows
        hits, ids = [], set()
        for row in iter_table(zf, "DRUG"):
            dn = (row.get("drugname") or "") + "|" + (row.get("prod_ai") or "")
            for drug, pat in TARGETS.items():
                if pat.search(dn):
                    pid = row["primaryid"]
                    ids.add(pid)
                    hits.append({
                        "quarter": qtag, "primaryid": pid, "caseid": row.get("caseid"),
                        "target_drug": drug, "role_cod": row.get("role_cod"),
                        "drugname": row.get("drugname"), "prod_ai": row.get("prod_ai"),
                        "route": row.get("route"), "dose_amt": row.get("dose_amt"),
                        "dose_unit": row.get("dose_unit"),
                    })
                    break
        if not ids:
            return qtag, 0, 0, 0
        # 2) DEMO / REAC / OUTC filtered to target primaryids
        demo_rows, reac_rows, outc_rows = [], [], []
        for row in iter_table(zf, "DEMO"):
            if row["primaryid"] in ids:
                demo_rows.append({k: row.get(k) for k in
                    ("primaryid", "caseid", "caseversion", "i_f_code", "event_dt",
                     "fda_dt", "rept_dt", "age", "age_cod", "age_grp", "sex",
                     "reporter_country", "occr_country", "occp_cod", "report_type",
                     "lit_ref", "mfr_sndr")})
        for row in iter_table(zf, "REAC"):
            if row["primaryid"] in ids:
                reac_rows.append({"primaryid": row["primaryid"], "pt": row.get("pt")})
        for row in iter_table(zf, "OUTC"):
            if row["primaryid"] in ids:
                outc_rows.append({"primaryid": row["primaryid"], "outc_cod": row.get("outc_cod")})
        for name, rows in (("drugrows", hits), ("reports", demo_rows),
                           ("reactions", reac_rows), ("outcomes", outc_rows)):
            with open(os.path.join(OUT, f"{name}_{qtag}.csv"), "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader(); w.writerows(rows)
        return qtag, len(ids), len(reac_rows), len(hits)

if __name__ == "__main__":
    zips = sorted(glob.glob(os.path.join(BASE, "faers_ascii_*.zip")))
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for z in zips:
        if only and only not in z:
            continue
        try:
            qtag, n, nreac, ndrug = parse_quarter(z)
            print(f"[{qtag}] target reports={n}, reaction rows={nreac}, drug rows={ndrug}", flush=True)
        except Exception as e:
            print(f"[{os.path.basename(z)}] ERROR: {type(e).__name__}: {e}", flush=True)
    print("=== PARSING DONE ===")
