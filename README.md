<div align="center">

# Cool the Heat, Not Just the Air

**A framework for making site-specific cooling decisions in outdoor and semi-outdoor public spaces**

HiMCM 2026 · Problem B · Team #19386

![Python](https://img.shields.io/badge/python-3.13-blue)
![Contest](https://img.shields.io/badge/HiMCM-2026_Problem_B-1f3a5f)
![Team](https://img.shields.io/badge/Team-%2319386-555555)

[Overview](#overview) · [Results](#results) · [Method](#method) · [Robustness](#robustness-analyses) · [Repository structure](#repository-structure) · [Reproduce](#reproducing-the-official-results) · [Limitations](#limitations-and-caveats)

</div>

---

## Overview

Installing more powerful cooling equipment does not by itself reduce the heat that people actually experience in open or semi-open spaces. Cooled air disperses, solar load and wind change the local environment, and an evaporative cooler that works well in dry air does little in humid air. This project asks which cooling technology, and at what scale, reduces *occupant heat exposure* at a given site.

We compare a no-cooling baseline (`Baseline-0`) with three systems at four venues:

| System | Description |
| :--- | :--- |
| **A** | Recirculating fans. They raise local air speed and do not lower dry-bulb temperature. |
| **B** | Medium-capacity evaporative cooling. |
| **C** | High-capacity evaporative cooling. |

| Venue | Climate and layout |
| :--- | :--- |
| Hanoi Bus Shelter | Hot and humid. Small roadside structure, roof of about 6.4 m². |
| Hanoi Railway Station | Hot and humid. Large semi-enclosed hall with strong natural ventilation. |
| Honolulu International Airport | Most enclosed and shaded of the four. Modelled exposure is far lower than at the other venues. |
| Sun Devil Stadium, Phoenix | Hot and dry. Fully open air with intense solar load. |

A configuration `X-n` means `n` units of system `X`. The official comparison allows at most 40 units per configuration (`cap40`). The Hanoi sites share one weather series, so three hourly weather files serve four venues.

The submitted solution paper is in [`paper/19386.pdf`](paper/19386.pdf).

## Results

Official deterministic results of model V2.4.0 under `cap40`, evaluated over 7,728 operating hours per venue (seven summers, June to August 2020–2026, 09:00–21:00). The selected configuration is the top-ranked option after Pareto screening and AHP.

| Venue | Selected | Danger share, baseline → selected | EDH38 reduction | η<sub>EDH38</sub> | Electricity (kWh/season) | Annualised cost (USD/season) | AHP score |
| :--- | :---: | :---: | :---: | :---: | ---: | ---: | :---: |
| Hanoi Bus Shelter | `C-7` | 45.22% → 35.27% | 43.63% | 0.185 | 19,320 | 10,271 | 0.8590 |
| Hanoi Railway Station | `B-40` | 56.84% → 49.50% | 27.02% | 0.820 | 52,992 | 25,334 | 0.8101 |
| Honolulu Airport | `Baseline-0` | 0.00% → 0.00% | not applicable | undefined (0/0) | 0 | 0 | 0.7730 |
| Sun Devil Stadium | `C-40` | 66.31% → 50.25% | 30.03% | 15.59 | 110,400 | 60,948 | 0.6507 |

Reading the table:

- **Climate decides what works.** In humid Hanoi the wet-bulb depression is small, which caps how far evaporative units can cool supply air. In dry Phoenix it is large, and System C stays top-ranked at the stadium in all 56 one-at-a-time tests.
- **Cooling hardware is not enough.** Even the best configuration leaves about half of person-hours above the danger threshold at the stadium and at the railway station. The paper pairs equipment with shade, drinking water, heat monitoring and medical response.
- **Bus Shelter: `C-7` is a mathematical optimum, not a physical one.** One System C unit has an effective service area of about 20 m², roughly three times the 6.4 m² floor, and the model does not limit the unit count. We recommend a feasibility-checked pilot of 1–2 units rather than seven.
- **Honolulu:** modelled EDH38 is already zero in the operating window, so EDH38 gives no case for added equipment and η is undefined rather than zero.
- AHP scores are min–max normalised within each venue's Pareto set, so they are **not comparable across venues**. Prices come from standardised workbook assumptions, not supplier quotes.

## Method

The model is a reduced-order, one-zone, quasi-steady hourly heat-and-moisture balance. It is a comparison tool, not a CFD model or a validated building model.

```text
Hourly weather (Open-Meteo)  +  occupancy profile  +  venue geometry and envelope
                                   │
                                   ▼
            baseline zone state: T_base, humidity ratio, ventilation
            (envelope conduction, sol-air gains, internal heat, wind + stack flow)
                                   │
                                   ▼
        cooling response, per system A / B / C and unit count n
        (fan air speed · evaporative supply-air drop · coverage · mixing)
                                   │
                                   ▼
        served / unserved occupants scored with AT_mod  →  danger share, comfort share, EDH38
                                   │
                                   ▼
        electricity, water, annualised cost  →  Pareto filter (4 objectives)  →  AHP ranking
```

**Key relations** (symbols as in the paper and in `PhysicalModel_Team19386_14Nov_v2.4.0.py`):

Evaporative supply-air drop, bounded by the wet-bulb depression:

$$\Delta T_{\mathrm{supply}} = \max\!\big[0,\ \min\!\big(\Delta T_{\mathrm{rated}},\ \varepsilon_{\mathrm{evap}}\,(T_{\mathrm{amb}} - T_{\mathrm{wb}})\big)\big]$$

Coverage and mixing with background ventilation:

$$f_{\mathrm{cov}} = \min\!\Big(1,\ \frac{n\,A_{\mathrm{cover}}}{A_{\mathrm{floor}}}\Big), \qquad \varphi = \frac{\mathrm{mix\_eff}\,\dot V_{\mathrm{sys}}}{\dot V_{\mathrm{sys}} + f_{\mathrm{cov}}\,\dot V_{\mathrm{inf}}}, \qquad T_{\mathrm{served}} = T_{\mathrm{base}} - \varphi\,\Delta T_{\mathrm{supply}}$$

Team-modified apparent temperature (an internal comparison index, not the Steadman index, UTCI or the NWS Heat Index):

$$AT_{\mathrm{mod}} = T + 0.348\,e - 0.70\,v + \frac{0.70\,Q_g}{v + 10} - 4.25, \qquad Q_g = 0.1333\, f_{\mathrm{sun}}\, I_{\mathrm{sol}}$$

Exposure and efficiency metrics:

$$EDH_{38} = \frac{1}{7}\sum_{y=2020}^{2026}\sum_{t\in y}\Big[N^{\mathrm{served}}_t \max(0, AT^{\mathrm{served}}_t - 38) + N^{\mathrm{unserved}}_t \max(0, AT^{\mathrm{unserved}}_t - 38)\Big]\Delta t, \qquad \eta_{EDH_{38}} = \frac{EDH_{38}^{\mathrm{base}} - EDH_{38}^{\mathrm{config}}}{E_{\mathrm{kWh}}}$$

- **Danger share**: occupancy-weighted fraction of person-hours with AT<sub>mod</sub> ≥ 38 °C. **Comfort share**: the same with AT<sub>mod</sub> ≤ 28 °C.
- **EDH38** is in person·°C·h per average summer. **E<sub>kWh</sub>** assumes constant rated power for 92 days × 12 h.
- The semi-open railway station reuses the same conservation framework with venue-specific areas, openings, U-values, sheltering and occupancy.
- **Decision rule.** Candidates are filtered on four objectives (comfort share up; danger share, annualised cost and electricity down), then ranked with AHP weights Safety 0.4236, Comfort 0.2270, Cost 0.2270, Energy 0.1223 (consistency ratio CR = 0.003838).

> The submitted PDF prints 0.33 for the vapour-pressure coefficient in Eq. (22). The V2.4.0 code, and therefore every reported number, uses 0.348.

## Robustness analyses

These runs supplement the official `cap40` results and do not replace them.

| Analysis | Setup | Outcome |
| :--- | :--- | :--- |
| One-at-a-time sensitivity | 56 perturbed cases per venue | Top-ranked configuration changes in 13 (Bus Shelter), 5 (Railway Station), 2 (Honolulu) and 0 (Stadium) of 56 cases. |
| Heatwave stress test | Ambient temperature +3 °C, all else fixed | Baseline danger share rises to 81.54%, 88.38%, 0.95% and 84.79% at the Bus Shelter, Railway Station, Honolulu and Stadium. |
| Monte Carlo | 10,000 draws per venue; seven-summer block resampling; independent multipliers in [0.8, 1.2] on ten parameter groups; seed 42 | Mean danger share, baseline → selected: Bus Shelter 46.11% → 36.06%, Railway Station 56.91% → 49.63%, Honolulu 0.00% → 0.00%, Stadium 66.42% → 50.41%. |

Monte Carlo means differ from the deterministic values by 0.05–0.90 percentage points, and the draw-to-draw standard deviation of danger share is 2.4–3.3 percentage points. Upper clipping at 1.0 on parameters such as `alpha_shelter`, `eps_evap`, `mix_eff` and `f_dung` is a computational bound, so the tails of the distributions should be read with care.

**Exploratory extension (not official).** Raising the cap to 120 units moves the AHP leaders to `C-42` (Bus Shelter), `B-40` (Railway Station), `B-80` (Honolulu) and `C-104` (Stadium). Margins between rank 1 and rank 2 are only 0.0005–0.0023 AHP points, so neighbouring unit counts are practically tied and the exact count is not robustly identified.

## Repository structure

```text
.
├── paper/
│   └── 19386.pdf                      submitted solution (summary sheet, report, memo, blog post, references, AI-use report)
├── 01_Official_Model/                 official V2.4.0 files, byte-identical copies (SHA-256 verified)
│   ├── Original_V2.4.0/               PhysicalModel_…_v2.4.0.py, parameter workbook, three hourly weather files
│   └── Official_Results/              ConfigurationResults, ParetoFrontier, AHP ranking, sensitivities, audits (CSV)
├── 02_Extended_Analysis/              EXPLORATORY: cap-120 scan, Monte Carlo, heatwave, rank tracking, figures
│   ├── CSV_Tables/   Monte_Carlo/   Sensitivity_Analysis/   Audit_Tables/   Figures_300dpi/
├── 03_Literature/                     contest statement, device specifications, literature verification register
├── 04_Audit_and_Corrections/          technical audit, C7 coefficient audit, corrected efficiency table, correction note
├── 05_Code_and_Reproducibility/
│   ├── Tests/                         reproduction test (13 checks)
│   ├── Source_Code/                   extended-analysis, Monte Carlo, audit and plotting scripts
│   ├── Historical_Logs/               logs of the 15 Nov runs (not produced by the test)
│   ├── Reproduction_Log.txt           raw output of the 16 Nov test run
│   └── requirements.txt
├── MANIFEST_SHA256.txt   99_Archive_HASHES.txt   CHANGELOG.md   RELEASE_STATUS.md
└── README.md
```

Venue identifiers in code and CSV files: `Tram_Xe_Buyt` (Bus Shelter), `Ga_Ha_Noi` (Railway Station), `San_Bay_Honolulu` (Honolulu Airport), `Sun_Devil_Stadium`. File names carry the author and date of the version, for example `…_LeMinh_15Nov_v1…`.

## Reproducing the official results

**Requirements.** Python 3.13 (tested on Linux) and the packages pinned in `05_Code_and_Reproducibility/requirements.txt`: `numpy`, `pandas`, `matplotlib`, `seaborn`, `reportlab`, `openpyxl`. No network access is needed. Windows and macOS have not been tested.

```bash
git clone https://github.com/lenghiemgiaphu/HiMCM_2026_ProblemB_Team19386.git
cd HiMCM_2026_ProblemB_Team19386

python3 -m venv .venv && source .venv/bin/activate
pip install -r 05_Code_and_Reproducibility/requirements.txt

python3 05_Code_and_Reproducibility/Tests/Test_Reproduction_LeMinh_16Nov_v1.py
```

The test copies the inputs to a temporary folder, runs model V2.4.0 with `--max-units 40`, compares the eight output CSVs with the official ones, re-computes the efficiency table, checks that the cap-120 table reproduces the official rows for n ≤ 40, and confirms by SHA-256 that no official file was modified. Expected result: **13 of 13 checks pass in about 6 seconds** (see `Reproduction_Log.txt`).

<details>
<summary>Running the model manually</summary>

V2.4.0 expects weather files with the suffixes `(1)`, `(1)` and `(2)`. Copy the three weather CSVs from `01_Official_Model/Original_V2.4.0/Weather_Inputs/` and the workbook into one folder, rename the weather files to

```text
WeatherData_Ga_Ha_Noi_GiaPhu_08Nov_v1(1).csv
WeatherData_Honolulu_International_Airport_GiaPhu_08Nov_v1(1).csv
WeatherData_Sun_Devil_Stadium_GiaPhu_08Nov_v1(2).csv
```

and run, quoting the file names in the shell:

```bash
python3 PhysicalModel_Team19386_14Nov_v2.4.0.py --data-dir <folder> \
  --workbook THONG_SO_THONG_NHAT_CHUAN_100_Team19386_14Nov_v2.4.0.xlsx --max-units 40
```

The workbook must be inside `--data-dir`. A `FileNotFoundError` on the weather files means the names do not match.
</details>

The extended-analysis scripts in `05_Code_and_Reproducibility/Source_Code/` read the environment variables `HIMCM_MODEL`, `HIMCM_RUN_DIR`, `HIMCM_WORKBOOK` and `HIMCM_OUT` (plus `HIMCM_TABLES` and `HIMCM_PLOTS` for figures). They were run on 15 November and were not re-run for the release, except for one Monte Carlo chunk (500 draws per venue), which reproduced identically. Monte Carlo seeding is `default_rng([42, site, chunk])`, with 80 chunks of 500 draws per venue.

## Data

Hourly weather comes from the [Open-Meteo Historical Weather API](https://open-meteo.com/en/docs/historical-weather-api): three series (Hanoi, Honolulu, Phoenix), each with 54,792 hourly records and 11 columns from 1 June 2020 to 31 August 2026. The audit finds no duplicate timestamps, no missing values and no imputed rows. The model uses June to August, hourly labels 10 to 21, which gives 12 operating hours × 92 days × 7 summers = 7,728 hours per venue. These are gridded model values, not on-site measurements. Please follow Open-Meteo's terms of use and attribution requirements when reusing the data.

## Limitations and caveats

- **Not calibrated or validated.** No venue was compared with field measurements. Results are conditional decision support, not certified safety predictions.
- **One-zone approximation.** The model resolves no spatial gradients, supply jets or occupant-level microclimates. Ventilation, sheltering, mixing and coverage coefficients are coarse.
- **Internal index.** AT<sub>mod</sub> is a comparison metric. Its 38 °C threshold is the team's choice and is not the NWS danger level.
- **Placeholder economics.** Costs use standardised workbook values, and electricity assumes constant rated power with no part-load behaviour.
- **No ablation study.** The release contains no controlled component-ablation runs, so no ablation percentages are reported.
- **Bus Shelter physics.** The model has no cap on unit count, so the Bus Shelter result needs an engineering feasibility check (roof loading, airflow, drainage, electrical protection, noise, clearance) before any purchase.
- **Name collision.** In the extended materials, "C7" also denotes a proposed railway-station wind-break case (`alpha_shelter` 0.35 → 0.20). It is unrelated to the configuration `C-7`, was never approved by the team, and does not enter any official result.

**Release status: CONDITIONAL PASS** (16 November 2026). Official files are byte-identical to the 14 November outputs and the executed checks pass. Open items recorded in [`RELEASE_STATUS.md`](RELEASE_STATUS.md): the origin of seven mismatched EDH38 cells in an internal AHP summary table (superseded by the CSV-derived values reported here), six of seven literature sources not read in full, unmodelled physical limits at the Bus Shelter, and no Windows test.

## Possible extensions

- Controlled ablation runs for stack ventilation, envelope conduction and latent moisture generation.
- A physical unit-count limit for the Bus Shelter, plus an instrumented pilot.
- Calibration against multi-sensor field data at one or more venues.
- Part-load and duty-cycle electricity models, and supplier-quoted prices.

## Team

Team #19386: Le Nghiem Gia Phu (team leader), Le Minh, Do Le Chi Hung and Tran Hung Phat.

## Citation

```bibtex
@misc{team19386_himcm2026,
  title        = {Cool the Heat, Not Just the Air: A Framework for Making Site-Specific
                  Decisions About Outdoor and Semi-Outdoor Public Spaces},
  author       = {{Team \#19386}},
  year         = {2026},
  note         = {HiMCM 2026, Problem B. Submitted solution paper},
  howpublished = {\url{https://github.com/lenghiemgiaphu/HiMCM_2026_ProblemB_Team19386}}
}
```

## License

A `LICENSE` file has not been added yet. Until one is, the code and documents in this repository are covered by default copyright. Weather data are derived from Open-Meteo and remain subject to its terms. The contest statement and device specifications in `03_Literature/` belong to their original owners (COMAP and the cited sources).
