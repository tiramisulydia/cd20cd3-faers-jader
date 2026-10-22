# Phase 4: publication-grade figures (BJH) + main tables
# Figures: F1 flow diagram, F2 fatal pathogen spectrum, F3 FAERS-vs-JADER forest, F4 TTO KM curves
# Style: Okabe-Ito colorblind-safe palette; 300 DPI PNG + vector PDF; 8-12pt fonts
import csv, json, os
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "..", "figures")
os.makedirs(OUT, exist_ok=True)

DRUGS = ["glofitamab", "epcoritamab", "mosunetuzumab"]
LBL = {"glofitamab": "Glofitamab", "epcoritamab": "Epcoritamab", "mosunetuzumab": "Mosunetuzumab"}
COL = {"glofitamab": "#0072B2", "epcoritamab": "#E69F00", "mosunetuzumab": "#009E73"}  # Okabe-Ito
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.linewidth": 0.8})

def save(fig, name):
    fig.savefig(os.path.join(OUT, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"{name} saved", flush=True)

# ---------- Figure 1: flow diagram ----------
fig, ax = plt.subplots(figsize=(7.6, 5.0))
ax.axis("off")
def box(x, y, w, h, text, fc="#F2F6FA", ec="#345", fs=8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012", fc=fc, ec=ec, lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)
def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13, color="#345", lw=1.1))

box(0.01, 0.74, 0.30, 0.20, "FAERS quarterly raw ASCII\n2018Q3\u20132026Q2 (32 packages)\n9,546 target-drug report\nrecords", fs=7.8)
box(0.36, 0.74, 0.26, 0.20, "Deduplication\n(case ID, latest version)\n6,973 unique cases", fs=8)
box(0.68, 0.74, 0.31, 0.20, "Primary-suspect (PS) set\n4,006 cases: glofitamab 1,477,\nepcoritamab 1,867,\nmosunetuzumab 662", fs=7.8)
box(0.68, 0.38, 0.31, 0.22, "Fatal reports (PS): 491 / 657 / 89\nFatal reports with \u22651 infection PT:\n111 / 188 / 51", fs=8)
box(0.36, 0.38, 0.26, 0.22, "Disproportionality analysis\nROR / PRR / BCPNN / EBGM\nfull-library, lymphoma- and\nfatal-restricted backgrounds", fs=7.8)
box(0.01, 0.38, 0.30, 0.22, "JADER (Oct 2026)\n1,060,670 cases\nepcoritamab 1,507,\nmosunetuzumab 214\n(glofitamab: none)", fs=7.8)
box(0.36, 0.02, 0.26, 0.22, "Sensitivity analyses\nCOVID-19 calendar strata\nTTO (Weibull/KM)\nco-reporting \u00b7 label audit", fs=7.8)
box(0.68, 0.02, 0.31, 0.22, "Cross-database\nvalidation (ROR,\nsuspect role)", fs=7.8)
arrow(0.31, 0.84, 0.36, 0.84); arrow(0.62, 0.84, 0.68, 0.84)
arrow(0.835, 0.74, 0.835, 0.60); arrow(0.68, 0.49, 0.62, 0.49)
arrow(0.31, 0.49, 0.36, 0.49); arrow(0.49, 0.38, 0.49, 0.24)
arrow(0.68, 0.13, 0.62, 0.13)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
save(fig, "Figure_1_flow")

# ---------- Figure 2: fatal infection pathogen spectrum (grouped bars) ----------
comp = json.load(open(os.path.join(BASE, "fatal_infection_composition.json"), encoding="utf-8"))
cats = ["viral", "bacterial", "site", "fungal", "opportunistic", "unspecified"]
cat_lbl = {"viral": "Viral", "bacterial": "Bacterial", "site": "Site-specific", "fungal": "Fungal", "opportunistic": "Opportunistic", "unspecified": "Unspecified"}
fig, ax = plt.subplots(figsize=(7.0, 4.2))
x = np.arange(len(cats)); w = 0.26
data = {}
for d in DRUGS:
    rec = comp[f"{d}_ps_only"]
    fc = rec["fatal_cats"]
    data[d] = [100 * fc.get(c, 0) / rec["fatal_inf"] for c in cats]
for j, d in enumerate(DRUGS):
    ax.bar(x + (j - 1) * w, data[d], w * 0.92, color=COL[d], label=f"{LBL[d]} (n={comp[f'{d}_ps_only']['fatal_inf']})", edgecolor="white", lw=0.4)
ax.set_xticks(x); ax.set_xticklabels([cat_lbl[c] for c in cats], fontsize=8.5)
ax.set_ylabel("% of fatal reports with \u22651 infection PT")
ax.set_ylim(0, max(max(v) for v in data.values()) * 1.2)
ax.legend(fontsize=8, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Pathogen categories among fatal reports with infection (PS set; multi-category per case possible)", fontsize=8.8)
save(fig, "Figure_2_spectrum")

# ---------- Figure 3: FAERS vs JADER forest ----------
rows = [
    # (label, drug, faers_ror, lo, hi, jader_ror, jlo, jhi or None)
    ("CMV infection reactivation", "epcoritamab", 94.16, 72.82, 121.74, 63.88, 48.02, 84.99),
    ("CMV test positive", "epcoritamab", 171.73, 115.12, 256.17, None, None, None),
    ("CMV enterocolitis", "epcoritamab", 114.22, 58.98, 221.18, 6.87, 3.55, 13.26),
    ("CMV viraemia", "epcoritamab", 46.22, 30.86, 69.24, 19.19, 13.05, 28.22),
    ("CMV infection", "epcoritamab", 26.75, 19.69, 36.34, 2.45, 1.35, 4.43),
    ("COVID-19 pneumonia", "epcoritamab", 10.15, 6.46, 15.96, 22.56, 14.85, 34.28),
    ("PML", "epcoritamab", 11.65, 5.82, 23.35, 7.68, 3.43, 17.20),
    ("PJP", "epcoritamab", 8.82, 4.74, 16.43, 2.53, 1.57, 4.09),
    ("Hypogammaglobulinaemia", "epcoritamab", None, None, None, 62.43, 45.69, 85.31),
    ("COVID-19 pneumonia", "mosunetuzumab", 55.20, 39.25, 77.62, 48.23, 22.63, 102.79),
    ("CMV infection reactivation", "mosunetuzumab", 25.17, 11.26, 56.25, None, None, None),
    ("RSV infection", "mosunetuzumab", 7.84, 2.52, 24.37, None, None, None),
    ("Disseminated tuberculosis", "glofitamab", 19.41, 6.24, 60.34, None, None, None),
    ("Varicella zoster virus infection", "glofitamab", 18.89, 7.07, 50.47, None, None, None),
    ("CMV infection reactivation", "glofitamab", 18.77, 10.07, 34.99, None, None, None),
    ("PML", "glofitamab", 7.35, 2.75, 19.61, None, None, None),
    ("Fungal infection", "glofitamab", 7.07, 4.32, 11.58, None, None, None),
    ("PJP", "glofitamab", 6.68, 3.00, 14.90, None, None, None),
]
fig, ax = plt.subplots(figsize=(7.4, 6.2))
ypos = np.arange(len(rows))[::-1]
for (label, drug, fr, flo, fhi, jr, jlo, jhi), y in zip(rows, ypos):
    c = COL[drug]
    if fr:
        ax.errorbar(fr, y, xerr=[[fr - flo], [fhi - fr]], fmt="o", ms=5.5, color=c, ecolor=c, elinewidth=1.2, capsize=2.5, zorder=3)
    if jr:
        ax.errorbar(jr, y - 0.18, xerr=[[jr - jlo], [jhi - jr]], fmt="s", ms=5, color=c, ecolor=c, elinewidth=1.2, capsize=2.5, mfc="white", zorder=3)
    ax.text(0.062, y, f"{label}", ha="right", va="center", fontsize=8.3)
    ax.text(0.064, y - 0.42, LBL[drug], ha="right", va="center", fontsize=6.8, color=c, style="italic")
ax.axvline(1, color="#666", lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_xlim(0.8, 320)
ax.set_xticks([1, 3, 10, 30, 100, 300]); ax.set_xticklabels(["1", "3", "10", "30", "100", "300"])
ax.set_yticks([])
ax.set_xlabel("Reporting odds ratio (95% CI), log scale")
ax.spines[["top", "right", "left"]].set_visible(False)
from matplotlib.lines import Line2D
handles = [Line2D([], [], marker="o", color="#333", ls="", ms=6, label="FAERS (PS, full library)"),
           Line2D([], [], marker="s", color="#333", mfc="white", ls="", ms=6, label="JADER (suspect)")]
ax.legend(handles=handles, fontsize=8, frameon=False, loc="lower right")
save(fig, "Figure_3_forest")

# ---------- Figure 4: TTO KM curves ----------
tto = list(csv.DictReader(open(os.path.join(BASE, "tto_detail.csv"), encoding="utf-8")))
fig, ax = plt.subplots(figsize=(6.6, 4.4))
for d in DRUGS:
    for fatal, ls, lab in (("0", "-", "all"),):
        xs = sorted(float(r["tto_days"]) for r in tto if r["drug"] == d)
        if not xs:
            continue
        xs = np.array(xs); y = np.arange(1, len(xs) + 1) / len(xs)
        ax.step(xs, y, where="post", color=COL[d], lw=1.8, label=f"{LBL[d]} (n={len(xs)})")
    xs_f = sorted(float(r["tto_days"]) for r in tto if r["drug"] == d and r["fatal"] == "1")
    if xs_f:
        xs_f = np.array(xs_f); y_f = np.arange(1, len(xs_f) + 1) / len(xs_f)
        ax.step(xs_f, y_f, where="post", color=COL[d], lw=1.1, ls="--")
ax.set_xscale("log")
ax.set_xlim(1, 1500)
ax.set_xticks([1, 7, 30, 90, 180, 365, 1500])
ax.set_xticklabels(["1", "7", "30", "90", "180", "365", "1500"])
ax.set_xlabel("Time to onset (days, log scale)")
ax.set_ylabel("Cumulative proportion of infection onsets")
ax.set_ylim(0, 1.02)
from matplotlib.lines import Line2D
h1 = [Line2D([], [], color=COL[d], lw=1.8, label=f"{LBL[d]}") for d in DRUGS]
h1.append(Line2D([], [], color="#555", lw=1.1, ls="--", label="Fatal subset (same colour)"))
ax.legend(handles=h1, fontsize=8, frameon=False, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Cumulative time-to-onset of infection reports (event-date axis, PS set)", fontsize=9.5)
save(fig, "Figure_4_tto")

print("=== ALL FIGURES DONE ===", flush=True)
