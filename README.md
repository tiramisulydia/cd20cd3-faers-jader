# Fatal infection pathogen spectra of CD20×CD3 bispecific antibodies — analysis code

Analysis code and data dictionaries accompanying the manuscript:

> **Fatal infection pathogen spectra of CD20×CD3 bispecific antibodies in B-cell lymphoma: a disproportionality analysis of FAERS with JADER cross-validation**

## Data sources (publicly available)

- **FAERS** quarterly ASCII files, 2018Q3–2026Q2: https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
- **JADER** (Pharmaceuticals and Medical Devices Agency, Japan), October 2026 snapshot: https://www.pmda.go.jp/english/

## Repository contents

| File | Role |
|---|---|
| `parse_quarters.py` | FAERS quarterly ASCII parser — target-drug reports (report-level + reaction-level) |
| `combine_reports.py` | Combine parsed quarters and deduplicate by case identifier |
| `build_background.py` | Full-database background counts (denominators for disproportionality) |
| `analyze_main.py` | Main four-method disproportionality analysis (ROR / PRR / BCPNN / EBGM) |
| `analyze_sensitivity.py` | Lymphoma-restricted and fatal-restricted sensitivity analyses |
| `analyze_covid_calendar.py` | COVID-19 calendar stratification (event-date and receipt-date axes) |
| `parse_ther.py` | THER (therapy date) table parsing |
| `analyze_tto.py` | Time-to-onset analysis with Weibull shape modelling |
| `analyze_hypogamma.py` | Hypogammaglobulinaemia / neutropenia / CRS co-reporting analysis |
| `jader_pipeline.py` | JADER cross-validation pipeline (incl. MedDRA/J mapping) |
| `make_figures.py` | Publication-grade figures |
| `assemble_supplementary.py` | Supplementary tables S1–S6 |
| `Table_S5a_infection_pt_dictionary.csv` | Prespecified infection Preferred Term (PT) dictionary |
| `Table_S5b_meddraj_mapping.csv` | MedDRA/J → English PT mapping used for JADER |

## Requirements

Python 3.13 (numpy 2.5.2, scipy 1.18.1, matplotlib for figures). Random seeds are fixed throughout for reproducibility.

## Notes

- FAERS and JADER are de-identified public databases; no ethical approval was required.
- Signal detection required concordance of all four methods (ROR lower 95% CI > 1 with ≥ 3 cases; PRR ≥ 2 with χ² ≥ 4; IC025 > 0; EBGM05 ≥ 2), per the READUS-PV reporting guideline.

## License

Code: MIT License (see `LICENSE`). Data dictionaries: CC BY 4.0.

## Citation

If you use this code, please cite the manuscript and this repository (Zenodo DOI: [to be added on first release]).
