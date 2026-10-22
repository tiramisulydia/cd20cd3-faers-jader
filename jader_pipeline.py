# JADER cross-validation pipeline for the CD20xCD3 bispecific fatal-infection study.
# Prereq: user manually downloads JADER bulk zip (CAPTCHA-gated) from
#   https://www.info.pmda.go.jp/fukusayoudb/CsvDownload.jsp
# and places it (or the extracted folder) under  outputs/jader_raw/
# Tables: demo/drug/reac/hist, Shift-JIS (cp932), comma-separated.
# PS-equivalent = 疑わしい医薬品 (suspect). Fatal = reac row outcome 死亡.
import csv, glob, json, math, os, re, sys, zipfile
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "..", "jader_raw")

DRUG_RX = {
    "glofitamab": re.compile(r"グロフィタマブ|コールムヴィ|コルムヴィ|glofitamab", re.I),
    "epcoritamab": re.compile(r"エプコリタマブ|エポリタマブ|エプキンリ|epcoritamab", re.I),
    "mosunetuzumab": re.compile(r"モスネツズマブ|ルンスミオ|mosunetuzumab", re.I),
}

# MedDRA/J (Japanese) -> English PT, curated for infection PTs of interest
# (spellings verified against JADER 2026-10 real data; full-width COVID variants included)
PT_JA_EN = {
    "サイトメガロウイルス感染": "CMV infection",
    "サイトメガロウイルス感染再活性化": "CMV infection reactivation",
    "サイトメガロウイルス感染再燃": "CMV infection reactivation",
    "サイトメガロウイルスウイルス血症": "CMV viraemia",
    "サイトメガロウイルス血症": "CMV viraemia",
    "サイトメガロウイルス腸炎": "CMV enterocolitis",
    "サイトメガロウイルス性腸炎": "CMV enterocolitis",
    "サイトメガロウイルス大腸炎": "CMV enterocolitis",
    "サイトメガロウイルス検査陽性": "CMV test positive",
    "サイトメガロウイルス性肺炎": "CMV pneumonia",
    "サイトメガロウイルス性脈絡網膜炎": "CMV chorioretinitis",
    "サイトメガロウイルス性脳炎": "CMV encephalitis",
    "サイトメガロウイルス脳炎": "CMV encephalitis",
    "COVID-19": "COVID-19",
    "ＣＯＶＩＤ－１９": "COVID-19",
    "COVID-19肺炎": "COVID-19 pneumonia",
    "ＣＯＶＩＤ－１９肺炎": "COVID-19 pneumonia",
    "敗血症": "Sepsis",
    "敗血症性ショック": "Septic shock",
    "肺炎": "Pneumonia",
    "誤嚥性肺炎": "Pneumonia aspiration",
    "細菌性肺炎": "Bacterial pneumonia",
    "ニューモシスチス・イロベチイ肺炎": "Pneumocystis jirovecii pneumonia",
    "ニューモシスチス・イロベチ肺炎": "Pneumocystis jirovecii pneumonia",
    "ニューモシスチス・イロベツィ肺炎": "Pneumocystis jirovecii pneumonia",
    "ニューモシスチス肺炎": "Pneumocystis jirovecii pneumonia",
    "ニューモシスチス・カリーニ肺炎": "Pneumocystis jirovecii pneumonia",
    "帯状疱疹": "Herpes zoster",
    "水痘帯状疱疹ウイルス感染": "Varicella zoster virus infection",
    "水痘・帯状疱疹ウイルス感染": "Varicella zoster virus infection",
    "真菌感染": "Fungal infection",
    "真菌性肺炎": "Fungal pneumonia",
    "アスペルギルス感染": "Aspergillus infection",
    "全身性カンジダ症": "Systemic candida",
    "進行性多巣性白質脳症": "Progressive multifocal leukoencephalopathy",
    "結核": "Tuberculosis",
    "肺結核": "Tuberculosis",
    "播種性結核": "Disseminated tuberculosis",
    "非定型マイコバクテリア感染": "Atypical mycobacterial infection",
    "B型肝炎再活性化": "Hepatitis B reactivation",
    "Ｂ型肝炎再活性化": "Hepatitis B reactivation",
    "呼吸器合胞体ウイルス感染": "Respiratory syncytial virus infection",
    "緑膿菌敗血症": "Pseudomonal sepsis",
    "緑膿菌感染": "Pseudomonas infection",
    "菌血症": "Bacteraemia",
    "黄色ブドウ球菌菌血症": "Staphylococcal bacteraemia",
    "ブドウ球菌感染": "Staphylococcal infection",
    "ブドウ球菌敗血症": "Staphylococcal sepsis",
    "髄膜炎": "Meningitis",
    "蜂窩織炎": "Cellulitis",
    "感染症": "Infection",
    "感染": "Infection",
    "細菌感染": "Bacterial infection",
    "ウイルス感染": "Viral infection",
    "低ガンマグロブリン血症": "Hypogammaglobulinaemia",
    "発熱性好中球減少症": "Febrile neutropenia",
    "好中球減少症": "Neutropenia",
    "サイトカイン放出症候群": "Cytokine release syndrome",
    "尿路感染": "Urinary tract infection",
    "大腸菌感染": "Escherichia infection",
    "腸球菌感染": "Enterococcal infection",
    "ヘルペスウイルス感染": "Herpes virus infection",
    "単純ヘルペス": "Herpes simplex",
    "インフルエンザ": "Influenza",
    "気管支肺アスペルギルス症": "Bronchopulmonary aspergillosis",
    "気管支肺アスペルギルス感染": "Bronchopulmonary aspergillosis",
    # ---- round 2 additions (spelling variants from real JADER 2026-10) ----
    "シュードモナス感染": "Pseudomonas infection",
    "シュードモナス敗血症": "Pseudomonal sepsis",
    "ブドウ球菌性敗血症": "Staphylococcal sepsis",
    "カンジダ感染": "Candida infection",
    "エプスタイン・バーウイルス感染": "Epstein-Barr virus infection",
    "エプスタイン・バーウイルス感染再燃": "Epstein-Barr virus infection reactivation",
    "コロナウイルス感染": "Coronavirus infection",
    "肺炎球菌性肺炎": "Pneumococcal pneumonia",
    "肺炎球菌感染": "Pneumococcal infection",
    "肺炎球菌性菌血症": "Pneumococcal bacteraemia",
    "クリプトコッカス性肺炎": "Cryptococcal pneumonia",
    "マイコバクテリウム・アビウムコンプレックス感染": "Mycobacterium avium complex infection",
    "播種性帯状疱疹": "Disseminated herpes zoster",
    "好中球数減少": "Neutrophil count decreased",
    "白血球数減少": "White blood cell count decreased",
    "免疫エフェクター細胞関連神経毒性症候群": "ICANS",
    "血球貪食性リンパ組織球症": "Haemophagocytic lymphohistiocytosis",
    "腫瘍フレア": "Tumour flare",
    # ---- round 3 additions ----
    "低γグロブリン血症": "Hypogammaglobulinaemia",
    "低ガンマグロブリン血症": "Hypogammaglobulinaemia",
    "ＲＳウイルス感染": "Respiratory syncytial virus infection",
    "RSウイルス感染": "Respiratory syncytial virus infection",
    "シュードモナス性肺炎": "Pseudomonas pneumonia",
    "細菌性敗血症": "Bacterial sepsis",
    "全身性真菌症": "Systemic mycosis",
    "ＢＫウイルス感染": "BK virus infection",
    "ブドウ球菌性肺炎": "Staphylococcal pneumonia",
    "結核性胸膜炎": "Tuberculous pleurisy",
    "骨結核": "Bone tuberculosis",
    "ヘルペス脳炎": "Herpes encephalitis",
}

LYMPHOMA_RX = re.compile(r"リンパ腫")

def find_tables():
    """Return dict table -> path or (zip, member)."""
    out = {}
    if not os.path.isdir(RAW):
        return None
    # case 1: zip file present
    zips = glob.glob(os.path.join(RAW, "*.zip"))
    for z in zips:
        with zipfile.ZipFile(z) as zf:
            for n in zf.namelist():
                low = os.path.basename(n).lower()
                for t in ("demo", "drug", "reac", "hist"):
                    if low.startswith(t) and low.endswith(".csv"):
                        out.setdefault(t, (z, n))
    # case 2: extracted csvs
    for f in glob.glob(os.path.join(RAW, "**", "*.csv"), recursive=True):
        low = os.path.basename(f).lower()
        for t in ("demo", "drug", "reac", "hist"):
            if low.startswith(t):
                out.setdefault(t, f)
    return out if len(out) == 4 else out

def iter_rows(loc):
    """Yield dict rows with cp932 decoding; loc = path or (zip, member)."""
    if isinstance(loc, tuple):
        z, member = loc
        with zipfile.ZipFile(z) as zf, zf.open(member) as fh:
            text = (b.decode("cp932", errors="replace") for b in fh)
            yield from csv.DictReader(text)
    else:
        with open(loc, encoding="cp932", errors="replace", newline="") as fh:
            yield from csv.DictReader(fh)

def col(row, *cands):
    # exact match first (avoid 有害事象連番 shadowing 有害事象), then substring
    # (excluding serial-number columns like 原疾患連番)
    for c in cands:
        for k in row:
            if k and k.strip() == c:
                return row[k]
    for c in cands:
        for k in row:
            if k and c in k and "連番" not in k:
                return row[k]
    return ""

def ror_ci(a, b, c, d):
    v = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return v, v * math.exp(-1.96 * se), v * math.exp(1.96 * se)

if __name__ == "__main__":
    tabs = find_tables()
    if not tabs or len(tabs) != 4:
        print("!! JADER tables not found. Put the bulk zip or extracted CSVs into outputs/jader_raw/", flush=True)
        print("   found:", list(tabs) if tabs else "nothing", flush=True)
        sys.exit(1)
    print("tables:", {k: (v if isinstance(v, str) else v[1]) for k, v in tabs.items()}, flush=True)

    # ---- pass 1: DRUG -> target cases (suspect role) ----
    case_drugs = defaultdict(set)      # caseid -> {drugs}
    case_drugrole = defaultdict(set)   # caseid -> {(drug, role)}
    for row in iter_rows(tabs["drug"]):
        cid = (col(row, "識別番号") or "").strip()
        name = " ".join([col(row, "一般名"), col(row, "販売名")])
        role = (col(row, "関与") or "").strip()
        for drug, rx in DRUG_RX.items():
            if rx.search(name):
                case_drugrole[cid].add((drug, role))
                if "被疑" in role or "疑わしい" in role:  # 被疑薬 = suspect (PS-equivalent)
                    case_drugs[cid].add(drug)
                break
    n_suspect = len(case_drugs)
    print(f"[drug] suspect-role target cases={n_suspect} "
          f"(glofitamab={sum(1 for s in case_drugs.values() if 'glofitamab' in s)}, "
          f"epcoritamab={sum(1 for s in case_drugs.values() if 'epcoritamab' in s)}, "
          f"mosunetuzumab={sum(1 for s in case_drugs.values() if 'mosunetuzumab' in s)})", flush=True)

    # ---- pass 2: REAC -> all-case PT counts (background) + target PTs + fatal ----
    N = 0
    bg_pt_cases = defaultdict(set)     # PT(en or ja) -> set(caseid)
    tgt_pt_cases = defaultdict(lambda: defaultdict(set))  # drug -> PT -> set(caseid)
    fatal_cases = set()
    lymph_cases = set()
    unmatched = defaultdict(int)
    seen_cases = set()
    for row in iter_rows(tabs["reac"]):
        cid = (col(row, "識別番号") or "").strip()
        pt = (col(row, "有害事象") or "").strip()
        outc = (col(row, "転帰") or "").strip()
        if not cid or not pt:
            continue
        seen_cases.add(cid)
        en = PT_JA_EN.get(pt)
        bg_pt_cases[en or pt].add(cid)
        if not en and any(x in pt for x in ("感染", "肺炎", "敗血", "サイトメガロ", "ヘルペス",
                                            "帯状疱疹", "真菌", "結核", "COVID", "ニューモシスチス")):
            unmatched[pt] += 1
        if "死亡" in outc:
            fatal_cases.add(cid)
        for drug in case_drugs.get(cid, ()):
            tgt_pt_cases[drug][en or pt].add(cid)
    N = len(seen_cases)
    print(f"[reac] total cases={N}, fatal cases={len(fatal_cases)}, mapped PTs={len(bg_pt_cases)}", flush=True)
    if unmatched:
        print("[reac] unmatched infection-like JA PTs (top 20, review for map additions):", flush=True)
        for pt, n in sorted(unmatched.items(), key=lambda x: -x[1])[:20]:
            print(f"   {pt}: {n}", flush=True)

    # ---- pass 3: HIST -> lymphoma flag ----
    for row in iter_rows(tabs["hist"]):
        cid = (col(row, "識別番号") or "").strip()
        dis = (col(row, "原疾患") or "").strip()
        if LYMPHOMA_RX.search(dis):
            lymph_cases.add(cid)
    print(f"[hist] lymphoma cases={len(lymph_cases)}", flush=True)

    # ---- ROR per drug per mapped PT (suspect role) vs full JADER background ----
    rows_out = []
    for drug in ("glofitamab", "epcoritamab", "mosunetuzumab"):
        dset = {c for c, s in case_drugs.items() if drug in s}
        nD = len(dset)
        for pt, cset in tgt_pt_cases[drug].items():
            a = len(cset)
            if a < 1:
                continue
            c_all = len(bg_pt_cases.get(pt, set()))
            c = c_all - a
            b = nD - a
            d = N - a - b - c
            if b <= 0 or c <= 0 or d <= 0:
                continue
            ror, lo, hi = ror_ci(a, b, c, d)
            n_fatal = len(cset & fatal_cases)
            n_lymph = len(cset & lymph_cases)
            rows_out.append({
                "drug": drug, "pt": pt, "a": a, "N_drug": nD, "N_bg": N,
                "bg_count": c_all, "fatal_among": n_fatal, "lymphoma_among": n_lymph,
                "ROR": round(ror, 2), "ROR_lo": round(lo, 2), "ROR_hi": round(hi, 2),
                "signal_ror": int(lo > 1 and a >= 3),
            })
    fn = os.path.join(BASE, "jader_signals.csv")
    if not rows_out:
        print("!! no target-drug rows found — check drug name regex / role field", flush=True)
        sys.exit(1)
    with open(fn, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        w.writeheader(); w.writerows(sorted(rows_out, key=lambda r: (r["drug"], -r["ROR"])))
    print(f"\n=== JADER ROR DONE -> {fn} ===", flush=True)
    for drug in ("glofitamab", "epcoritamab", "mosunetuzumab"):
        sig = [r for r in rows_out if r["drug"] == drug and r["signal_ror"]]
        print(f"{drug}: suspect cases={sum(1 for s in case_drugs.values() if drug in s)}, "
              f"ROR signals={len(sig)}", flush=True)
        for r in sorted(sig, key=lambda x: -x["ROR"])[:15]:
            print(f"   {r['pt']}: a={r['a']} ROR={r['ROR']} ({r['ROR_lo']}-{r['ROR_hi']}) fatal={r['fatal_among']}", flush=True)
