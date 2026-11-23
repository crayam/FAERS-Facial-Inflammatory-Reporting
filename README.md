# FAERS Analysis of Facial Inflammatory Reporting in Atopic Dermatitis

This repository contains the analytical code for a retrospective active-comparator disproportionality study using FDA Adverse Event Reporting System (FAERS) quarterly ASCII data from **2022 Q2 to 2026 Q2 (17 quarters)**.

The study compares **dupilumab** with **tralokinumab** (primary active comparator) and **upadacitinib** (secondary active comparator).

## Study overview

Reports were deduplicated globally by CASEID, retaining the report version selected by CASEVERSION, FDA_DT, PRIMARYID and source quarter, followed by DELETE-record handling and the study date restriction. Eligibility required a primary-suspect (PS) drug role and a drug-linked indication exactly matching `Dermatitis atopic`, with matching drug and indication sequence numbers in the retained report.

The final eligible cohort comprised 164,515 dupilumab, 3,637 tralokinumab and 6,254 upadacitinib reports (total N = 174,406). Main adjusted analyses excluded 65 reports with missing geography (N = 174,341).

Outcomes were defined using the frozen Phenotype Dictionary v1.0. Primary11, Narrow6 and BroadSX23 were calculated as CASEID-level unions of their constituent preferred terms. Secondary outcomes included rosacea, Demodex-related conditions, blepharitis and seborrhoeic dermatitis.

## Repository contents

| File or directory | Purpose |
|---|---|
| `recompute_primary11.py` | Report selection, drug-specific cohort construction and Primary11 union |
| `recompute_narrow_broad.py` | Narrow6 and BroadSX23 outcome flags |
| `recompute_all_outcomes.py` | All seven prespecified outcome flags |
| `validate_primary11.py`, `validate_narrow_broad.py`, `validate_all_outcomes.py` | Independent cohort and outcome verification |
| `crude_comparative_analysis_v1_0.py` | Case-level 2×2 tables and historical crude ROR calculations |
| `table1/` | Report characteristics, missingness and standardized mean differences |
| `sparsity/` | Outcome-specific covariate sparsity assessment |
| `firth_primary11/` | Primary11 main Firth model and numerical validation |
| `firth_three/` | Narrow6 and BroadSX23 main-model preparation and validation |
| `final_stats/` | Calendar-year and age sensitivity analyses, secondary models and result integration |
| `restriction/` | HCP-only and US-only sensitivity analyses |
| `prr_pt/` | Supplementary PRR, Primary11 PT composition and overlap |
| `tto_primary11/` | Primary11 date extraction, date audit and descriptive time-to-onset analysis |
| `tto_correction_v1_1/` | Independent TTO verification, including R type-7 quantile checks |

## Statistical implementation

The main adjusted model was:

```text
Outcome ~ drug + sex + HCP_status + US_nonUS + calendar_quarter_index
```

Firth penalized logistic regression was implemented in **R using custom code**, rather than `logistf` or another ready-made Firth package. The penalized likelihood includes the Jeffreys-information penalty. Optimization uses base R BFGS with analytic adjusted scores and numerical-Hessian Newton refinement. Inference uses 95% profile penalized-likelihood confidence intervals and two-sided penalized likelihood-ratio tests. Independent Python/NumPy adjusted Fisher scoring provides numerical validation.

Drug is a three-level variable with dupilumab as the reference. Reported dupilumab-versus-comparator estimates invert the comparator coefficients. Sex and HCP missingness are separate categories. The continuous calendar-quarter index runs from 0 (2022 Q2) to 16 (2026 Q2).

The completed analyses comprise 21 three-drug Firth models: three main models, nine calendar-year/age sensitivity models, six restriction models and three exploratory secondary models. Demodex outcomes are descriptive only. P values are nominal, without multiplicity adjustment.

The historical crude script includes exact intervals for zero-cell tables. Final reporting marks non-finite conventional RORs as not estimable, without continuity correction; Demodex remains descriptive. The final integration logic is retained in `final_stats/`.

TTO analysis is exploratory and descriptive. It uses complete retained DEMO.EVENT_DT and index-drug-linked THER.START_DT dates, excludes ambiguous distinct nonblank start dates and negative intervals, and performs no date imputation or between-drug testing. NumPy linear quantiles were independently checked with R `quantile(type = 7)`. The v1.1 verification confirmed the original numerical TTO summaries.

## Software and dependencies

- **R 4.2.2** is recorded in the saved model results; `jsonlite` is used for JSON input/output.
- **Python 3.12.14** is documented for crude analysis and independent TTO verification; versions were not recorded for every historical stage.
- Python dependencies used across the scripts include `numpy`, `pandas`, `python-docx` and `openpyxl`, together with standard-library modules such as `sqlite3`.
- Historical package versions were not comprehensively recorded. NumPy 2.3.5 is documented for independent TTO verification, not necessarily for every earlier analysis.

## Inputs and execution

This is a **code-only research repository**. Raw FAERS ZIPs, the original phenotype dictionary DOCX, report-level intermediate datasets, saved results and stage-specific checksum files are not included.

The scripts retain the original project-relative layout, in which these repository contents were located in a `work/` directory alongside `inputs/` and `outputs/`. To use that layout, place this checkout at `PROJECT/work/` and run the relevant scripts from `PROJECT/`. Some scripts determine the project root from their own location; others use the current working directory.

Raw extraction requires the 17 source ZIPs and the frozen dictionary. Later stages require the corresponding retained-report masters, outcome flags, characteristics, validation records and checksum manifests. Certain historical validators also use prior workbooks. Inspect each script's input paths before execution. Do not bypass checksum assertions or change frozen definitions to accommodate missing inputs.

The approximate processing order is cohort construction → outcome unions → characteristics and crude analyses → model preparation, fitting and validation → sensitivity/secondary analyses → PRR/PT composition and TTO. This repository does not provide a self-contained, one-command reconstruction without the required inputs.

## Interpretation

FAERS analyses describe reporting odds, reporting patterns and disproportionate reporting. They do not estimate adverse-event incidence or establish causal effects. Reports may contain multiple preferred terms, and DEMO.EVENT_DT is a report-level event date rather than a validated onset date for a specific outcome term.
