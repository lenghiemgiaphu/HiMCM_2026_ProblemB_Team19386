<div align="center">

# Cool the Heat, Not Just the Air

**A framework for making site-specific decisions about outdoor and semi-outdoor public spaces**

COMAP High School Mathematical Contest in Modeling (HiMCM) 2026 · Problem B · Team #19386

[![Python](https://img.shields.io/badge/Python-3.13-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Model](https://img.shields.io/badge/Model-Deterministic%20V2.4.0%20(cap40)-0F766E?style=flat-square)](#2-physical-formulation)
[![Verification](https://img.shields.io/badge/Reproducibility%20Checks-13%2F13%20Passed-16A34A?style=flat-square)](./02_Code_and_Reproducibility/Team19386_Final_Technical_Package_16Nov2026_v1_LITE/05_Code_and_Reproducibility/Reproduction_Log.txt)
[![AHP Consistency](https://img.shields.io/badge/AHP%20CR-0.0038%20%3C%200.10-2563EB?style=flat-square)](#24-multi-criteria-ranking-pareto--ahp)
[![Monte Carlo](https://img.shields.io/badge/Monte%20Carlo-10%2C000%20draws%20(seed%2042)-7C3AED?style=flat-square)](#32-robustness-analyses)

*An hourly, quasi-steady, one-zone heat-and-moisture balance with a served sub-zone, scored with a team-modified apparent-temperature index, screened with a four-objective Pareto filter, and ranked with the Analytic Hierarchy Process (AHP).*

[Final report (PDF)](./01_Final_Report/19386.pdf) · [AHP ranking (CSV)](./04_Data_and_Results/AHP_Ranking_Official_Team19386_14Nov_v2.4.0.csv) · [Rank tracking (CSV)](./04_Data_and_Results/Sensitivity_RankTracking_LeMinh_15Nov_v1.csv) · [Validation checks (CSV)](./04_Data_and_Results/ValidationChecks_Team19386_14Nov_v2.4.0.csv)

</div>

---

## 1. Key findings

Open-air and semi-open public venues cool people while they are exposed to the surrounding atmosphere, so most conditioned air is lost to ambient exchange. A qualitative analogy is open-front refrigerated display cases, where infiltration is reported to account for the majority of the cooling load (NREL, 2024). The analogy motivates the question but is not used to validate the model. Rated capacity or kWh alone therefore says little about the relief people actually receive; what matters is the **reduction in occupant heat exposure**.

We evaluate **7,728 summer operating hours per venue** (June–August 2020–2026, 09:00–21:00) at four venues in three climates: the **Hanoi Bus Shelter**, the **Hanoi Railway Station**, **Honolulu International Airport** and **Sun Devil Stadium** in Phoenix. The two Hanoi venues share one weather series. Five conclusions follow.

1. **Efficiency is not effectiveness.** At the Bus Shelter, one fan (`A-1`) is about **3× more efficient per kWh** than one large evaporative cooler (`C-1`): η<sub>EDH38</sub> = 2.11 versus 0.70 person·°C·h/kWh. Yet it removes only **2.49 percentage points** of danger share, against **4.63 points** for `C-1`. Ranking by a single metric leads to the wrong purchase.
2. **Climate, not equipment, sets the value of evaporative cooling.** Supply-air cooling is bounded by the wet-bulb depression *D* = *T* − *T*<sub>wb</sub>. System C reaches its rated 8 °C drop in **93.1%** of Phoenix operating hours (mean *D* = 17.6 °C), versus **1.2%** in Hanoi (mean *D* = 3.8 °C) and **1.6%** in Honolulu (mean *D* = 5.6 °C). The efficiency η<sub>EDH38</sub> of the selected configuration is roughly **84× higher** at the stadium (`C-40`: 15.59) than at the Bus Shelter (`C-7`: 0.185).
3. **Efficiency saturates with coverage** (observation from the official configuration table, not in the submitted paper). With the coverage ratio ρ = *c*·*n*·*A*<sub>cov</sub> / *A*<sub>floor</sub>, η<sub>EDH38</sub> is constant while each added unit reaches new occupants (ρ ≤ 1) and then falls roughly as 1/ρ. For example, it is 15.59 for every System C count up to 104 at the stadium, and 0.820 for System B up to 40 units at the Railway Station, falling to 0.549 at 60 units. Beyond saturation the decline is slightly slower than 1/*n* because the mixing factor still rises with *n*.
4. **Site-specific optima** (484 configurations scanned under `cap40`, 290 Pareto-efficient):
   * **Sun Devil Stadium:** `C-40` (AHP rank 1, score 0.6507) cuts danger share from **66.31% to 50.25%** for USD 60,948 per season-year (about USD 3,795 per percentage point). An exploratory scan beyond the cap puts the ρ ≈ 1 point at `C-104` (danger share 24.55%).
   * **Hanoi Railway Station:** `B-40` (rank 1, score 0.8101) cuts danger share from **56.84% to 49.50%** for USD 25,334 per season-year.
   * **Hanoi Bus Shelter:** `C-7` is the mathematical AHP optimum (score 0.8590, danger share **45.22% to 35.27%**, USD 10,271), but seven units on a roof of about 6.4 m² is physically unproven. Two fans (`A-2`) already saturate coverage (danger share 41.23% for USD 308 per season-year). We recommend a feasibility-checked pilot of one or two units.
   * **Honolulu International Airport:** `Baseline-0` (score 0.7730). Modelled danger share is already 0.00%, so EDH38 gives no case for added equipment and η<sub>EDH38</sub> is undefined (0/0), not zero.
5. **Extreme heat needs operations, not just machines.** With ambient temperature raised by +3 °C, baseline danger share rises to **81.54%, 88.38% and 84.79%** at the Bus Shelter, Railway Station and Stadium. The selected configurations only bring this to 70.09%, 80.76% and 71.76%. Across the hottest 5% of operating hours (*T* ≥ 44.3 °C in Phoenix, 35.5 °C in Hanoi, 29.7 °C in Honolulu), 83–94% fall between hour labels 12 and 17. Cooling hardware must be paired with shade, drinking water, schedule changes, heat monitoring and medical response.

### Energy-to-exposure chain

```text
+-------------------+     +---------------------+     +------------------------------------+
|  Electricity E    | --> |  Cooling produced   | --> |  Reaches people                    |
|      (kWh)        |     |  (dT_supply, v_air) |     |  (f_people, mixing factor phi)     |
+-------------------+     +---------------------+     |  Avoided EDH38 / E = eta_EDH38     |
                                    |                 +------------------------------------+
                                    v
                          +---------------------+     +------------------------------------+
                          | Diluted by          | --> | Warm, moist ambient air re-enters  |
                          | ventilation V_inf   |     | (open-case analogy)                |
                          +---------------------+     +------------------------------------+
```

---

## 2. Physical formulation

The model is a reduced-order, one-zone, quasi-steady hourly balance. It is a comparison tool, not CFD and not a validated building model. Symbols follow the paper and `PhysicalModel_Team19386_14Nov_v2.4.0.py`.

### 2.1 Baseline heat and moisture balance

Each venue is one well-mixed background zone with no thermal mass, simulated for 7 × 92 × 12 = 7,728 hours (Δ*t* = 1 h). Moist-air properties follow Buck (saturation vapour pressure) and Stull (wet-bulb temperature). Heat gains combine roof and wall sol-air temperatures, glazing solar gain, internal loads (*q*<sub>s</sub> = 75 W/person, *q*<sub>lat</sub> = 55 W/person, *P*<sub>tb</sub> = 100 W) and wind-plus-buoyancy ventilation:

$$T_{\mathrm{sol,roof}} = T_{\mathrm{amb}} + \frac{\alpha_{\mathrm{roof}} I_{\mathrm{sol}}}{h_o} - K_{\mathrm{rad}}(1 - 0.8C), \qquad T_{\mathrm{sol,wall}} = T_{\mathrm{amb}} + \frac{\alpha_{\mathrm{wall}}\, f_{\mathrm{dung}}\, I_{\mathrm{sol}}}{h_o}$$

$$v_{\mathrm{ef}} = \max(k_{\mathrm{gio}}\,\alpha_{\mathrm{shelter}}\, v_{\mathrm{amb}},\ v_{\min}), \qquad \dot V_{\mathrm{inf}} = \sqrt{\dot V_{\mathrm{wind}}^2 + \dot V_{\mathrm{stack}}^2}$$

$$T_{\mathrm{base}} = \frac{H_{\mathrm{roof}} T_{\mathrm{sol,roof}} + H_{\mathrm{wall}} T_{\mathrm{sol,wall}} + H_{\mathrm{glass}} T_{\mathrm{amb}} + \dot Q_{\mathrm{solar,glass}} + \dot Q_{\mathrm{internal}} + H_{\mathrm{vent}} T_{\mathrm{amb}}}{H_{\mathrm{env}} + H_{\mathrm{vent}}}, \qquad H_{\mathrm{vent}} = \rho_{\mathrm{da}} c_{p,\mathrm{da}} \dot V_{\mathrm{inf}}$$

The background humidity ratio is capped at saturation in semi-enclosed venues. The open-air stadium branch sets *T*<sub>base</sub> = *T*<sub>amb</sub> and *w*<sub>zone</sub> = *w*<sub>amb</sub>:

$$w_{\mathrm{zone}} = \min\!\Big( w_{\mathrm{sat}}(T_{\mathrm{base}}),\ w_{\mathrm{amb}} + \frac{N_t\, q_{\mathrm{lat}}}{h_{fg}\rho_{\mathrm{da}} \max(\dot{V}_{\mathrm{inf}}, 0.05)} \Big)$$

### 2.2 Cooling systems

**Evaporative cooling** (System B: 5 °C, 1.2 kW; System C: 8 °C, 2.5 kW). Supply-air cooling is limited by both the rated drop and the wet-bulb depression, with effectiveness ε<sub>evap</sub> = 0.90:

$$\Delta T_{\mathrm{supply}} = \max\!\big[0,\ \min\!\big(\Delta T_{\mathrm{rated}},\ \varepsilon_{\mathrm{evap}}(T_{\mathrm{amb}} - T_{\mathrm{wb}})\big)\big]$$

$$f_{\mathrm{cov}} = \min\!\Big(1, \frac{n A_{\mathrm{cover}}}{A_{\mathrm{floor}}}\Big), \quad f_{\mathrm{people}} = \min(1,\ c\, f_{\mathrm{cov}}), \quad \varphi = \frac{\mathrm{mix\_eff}\,\dot V_{\mathrm{sys}}}{\dot V_{\mathrm{sys}} + f_{\mathrm{cov}} \dot V_{\mathrm{inf}}}, \quad T_{\mathrm{served}} = T_{\mathrm{base}} - \varphi\, \Delta T_{\mathrm{supply}}$$

Here *c* is the venue's occupant-clustering factor and mix_eff = 0.5. Evaporation also adds moisture, estimated from an adiabatic enthalpy balance, and the served air is capped at saturation.

**Recirculating fans** (System A: 0.25 kW). Fans raise local air speed by up to 1.2 m/s and do not change dry-bulb temperature or moisture. The effect fades linearly between 33 °C and 38 °C zone temperature:

$$f_{\mathrm{fan}} = \mathrm{clip}\!\Big(\frac{38 - T_{\mathrm{base}}}{38 - 33}, 0, 1\Big), \qquad v_{\mathrm{served}} = \sqrt{v_{\mathrm{background}}^2 + (1.2\, f_{\mathrm{fan}})^2}$$

### 2.3 Apparent temperature and exposure metrics

Exposure is scored with a team-modified apparent temperature AT<sub>mod</sub> (°C). It is an internal comparison index, not the Steadman index, UTCI or the NWS Heat Index. Person-hours are classed as comfort (AT<sub>mod</sub> ≤ 28 °C), stress, or danger (AT<sub>mod</sub> ≥ 38 °C):

$$AT_{\mathrm{mod}} = T + 0.348\, e - 0.70\, v + \frac{0.70\, Q_g}{v + 10} - 4.25, \qquad Q_g = 0.1333\, f_{\mathrm{sun}}\, I_{\mathrm{sol}}$$

> The submitted PDF prints 0.33 for the vapour-pressure coefficient in its Eq. (22). The V2.4.0 code, and therefore every reported number, uses 0.348.

Accumulated exposure and cooling efficiency, with electricity *E*<sub>kWh</sub> = *P*<sub>rated</sub> · *n* · 92 days · 12 h:

$$EDH_{38} = \frac{1}{7}\sum_{y=2020}^{2026}\sum_{t\in y}\Big[N^{\mathrm{served}}_t \max(0, AT^{\mathrm{served}}_t - 38) + N^{\mathrm{unserved}}_t \max(0, AT^{\mathrm{unserved}}_t - 38)\Big]\Delta t, \qquad \eta_{EDH_{38}} = \frac{EDH_{38}^{\mathrm{base}} - EDH_{38}^{\mathrm{config}}}{E_{\mathrm{kWh}}}$$

*Danger share* is the occupancy-weighted fraction of person-hours with AT<sub>mod</sub> ≥ 38 °C. EDH38 is in person·°C·h per average summer.

### 2.4 Multi-criteria ranking (Pareto + AHP)

Annualised cost combines capital recovery, maintenance, electricity and water:

$$C_{\mathrm{ann}} = n\,(C_{\mathrm{capex}} + C_{\mathrm{install}})(\mathrm{CRF} + r_{\mathrm{maint}}) + E_{\mathrm{kWh}}\, c_e + W_L\, c_w$$

with discount rate 6%, lifetime 10 years (CRF = 0.1359) and maintenance 5%. Prices are standardised workbook assumptions, not supplier quotes. The 290 non-dominated configurations of 484 (Pareto filter on comfort share, danger share, annualised cost and electricity) are min–max normalised **within each venue** and scored with AHP weights (consistency ratio CR = 0.003838 < 0.10):

1. **Safety** (minimise danger share): 0.4236
2. **Comfort** (maximise comfort share): 0.2270
3. **Cost** (minimise annualised cost): 0.2270
4. **Energy** (minimise seasonal electricity): 0.1223

Because normalisation is per venue, AHP scores are **not comparable across venues**.

---

## 3. Results and robustness

### 3.1 Official recommendations (`cap40`) and exploratory alternatives

Rows marked † lie outside the official `cap40` ranking (taken from the exploratory cap-120 table) and are shown for context only.

| Venue | Option | AHP score | Danger share (base → option) | Comfort share | η<sub>EDH38</sub> | Energy (MWh) / CO<sub>2</sub> (t) | Water (m³) | Annualised cost (USD) | USD per danger point |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Sun Devil Stadium**<br>hot-arid, open air | `C-40` official<br>`C-104` † | 0.6507<br>— | 66.31% → **50.25%**<br>66.31% → 24.55% | 2.63%<br>5.65% | 15.59<br>15.59 | 110.4 / 35.2<br>287.0 / 91.6 | 1,591<br>4,137 | 60,948<br>158,465 | 3,795<br>3,795 |
| **Hanoi Railway Station**<br>hot-humid, semi-open | `B-40` official<br>`A-100` † | 0.8101<br>— | 56.84% → **49.50%**<br>56.84% → 50.34% | 0.00%<br>0.00% | 0.820<br>0.815 | 53.0 / 34.9<br>27.6 / 18.2 | 344<br>0 | 25,334<br>15,393 | 3,448<br>2,367 |
| **Hanoi Bus Shelter**<br>hot-humid, roof ≈ 6.4 m² | `C-7` official<br>`A-2` † | 0.8590<br>— | 45.22% → **35.27%**<br>45.22% → 41.23% | 0.02%<br>0.00% | 0.185<br>1.69 | 19.3 / 12.7<br>0.55 / 0.36 | 129<br>0 | 10,271<br>308 | 1,032<br>77 |
| **Honolulu Airport**<br>mild, shaded | `Baseline-0` official<br>`C-14` † | 0.7730<br>— | 0.00% → **0.00%**<br>0.00% → 0.00% | 15.79%<br>18.36% | undefined (0/0)<br>0 (no EDH38 to avoid) | 0 / 0<br>38.6 / 26.1 | 0<br>392 | 0<br>32,622 | —<br>— |

Notes:

* `C-7` is a mathematical optimum: the model has no cap on unit count, and one System C unit has an effective service area of about 20 m², roughly three times the shelter floor. `A-2` is a practical alternative, not an official result.
* The headline recommendations of the paper are `C-40`, `B-40`, no extra cooling at Honolulu on EDH38 alone, and a 1–2-unit pilot at the Bus Shelter.
* Costs use placeholder workbook prices. Even the best configuration leaves roughly half of person-hours above the danger threshold at the stadium and the Railway Station.

### 3.2 Robustness analyses

These runs supplement the official results and do not replace them.

| Venue (selected) | One-at-a-time rank changes (of 56) | AHP weight perturbation (±20%, 9 cases) | Monte Carlo mean danger share, baseline → selected (SD) | +3 °C warming, baseline → selected |
| :--- | :---: | :---: | :---: | :---: |
| **Sun Devil Stadium** (`C-40`) | **0 / 56** | rank 1 retained 9 / 9 | 66.42% → 50.41% (2.7 pp) | 84.79% → 71.76% |
| **Hanoi Railway Station** (`B-40`) | **5 / 56** | rank 1 retained 9 / 9 | 56.91% → 49.63% (2.8 pp) | 88.38% → 80.76% |
| **Hanoi Bus Shelter** (`C-7`) | **13 / 56** | rank 1 retained 9 / 9 | 46.11% → 36.06% (3.3 pp) | 81.54% → 70.09% |
| **Honolulu Airport** (`Baseline-0`) | **2 / 56** | rank 1 retained 9 / 9 | 0.00% → 0.00% (0.0 pp) | 0.95% → 0.95% |

* **Monte Carlo design:** 10,000 draws per venue; each draw resamples the seven summer blocks with replacement and applies independent uniform multipliers in [0.8, 1.2] to ten parameter groups; seed 42. Sixteen configurations are tracked (baseline plus the best System A, B and C configuration at each of the four venues). Upper clipping at 1.0 on `alpha_shelter`, `eps_evap`, `mix_eff` and `f_dung` is a computational bound, so distribution tails should be read with care. Monte Carlo means differ from the deterministic values by 0.05–0.90 percentage points.
* **Automated verification (13 / 13 passed):** the test re-runs model V2.4.0 with `--max-units 40`, compares the eight output CSVs with the official ones (deviation 0), re-computes the efficiency table, checks that the cap-120 table reproduces the official rows for *n* ≤ 40, and confirms by SHA-256 that no official file changed. It runs in about 6 seconds.

---

## 4. Figures and verification dossiers

Each primary figure has an audit dossier (PDF) with its underlying data tables. File numbering follows the repository; in the paper these are Figures 8, 7 and 9.

### Pareto frontiers under `cap40` (paper Figure 8) — [dossier](./03_Verification_and_Figures/Figure_1_Verification_Dossier_V2.4.0.pdf)
![Four-objective cap40 Pareto frontiers](./03_Verification_and_Figures/fig01_pareto_tradeoff.png)

### EDH38 reduction and efficiency, with AT<sub>mod</sub> distributions (paper Figure 7) — [dossier](./03_Verification_and_Figures/Figure_2_Verification_Dossier_V2.4.0.pdf)
![EDH38 efficiency and AT_mod distributions](./03_Verification_and_Figures/fig02_edh38_efficiency.png)

### AHP weight perturbation and operating-season sensitivity (paper Figure 9) — [dossier](./03_Verification_and_Figures/Figure_3_Verification_Dossier_V2.4.0.pdf)
![AHP weight and operating-season sensitivity](./03_Verification_and_Figures/fig03_ahp_season_sensitivity.png)

---

## 5. Operator guidelines

Model rankings are decision support, not installation instructions. Before procurement, verify structure, power, water, drainage, airflow, maintenance and public-safety requirements on site.

1. **Match the technology to the wet-bulb depression.** An evaporative unit reaches its rated drop only when ε<sub>evap</sub>·*D* ≥ Δ*T*<sub>rated</sub>, that is *D* ≥ 5.6 °C for System B and *D* ≥ 8.9 °C for System C. Share of operating hours meeting this: Phoenix 97% (B) and 93% (C), Hanoi 19% and 1.2%. Prefer evaporative cooling in dry climates; in humid climates rely on shade, ventilation and air movement.
2. **Size to coverage.** Beyond ρ ≈ 1 added units buy airflow without reaching new people; the saturation count is *n*\* = *A*<sub>floor</sub> / (*c*·*A*<sub>cov</sub>).
3. **Venue-specific actions from the paper:**
   * *Bus Shelter:* start with a 1–2-unit pilot after checking roof loading, mounting clearance, airflow direction, drainage, electrical protection, noise and pedestrian clearance, and log local temperature, humidity and air speed against the model.
   * *Railway Station:* deploy System B in phases focused on high-density waiting zones, keep egress clear, monitor humidity, water supply and drainage, and pair equipment with shade, ventilation and heat alerts.
   * *Honolulu:* keep baseline operations in the modelled window, keep monitoring, and reassess if hours, occupancy or extreme weather change.
   * *Stadium:* because 50.25% of person-hours stay above the threshold even with 40 units, combine cooling zones with overhead shade, free drinking water, schedule adjustments, real-time heat monitoring and medical response. Run an instrumented pilot before expanding.
4. **Replace placeholder prices** with local supplier quotes and verify utility requirements against supplier specifications.

---

## 6. Repository structure

```text
public-space-cooling-model/
├── 01_Final_Report/
│   └── 19386.pdf                                  # submitted solution (summary sheet, report, memo, blog post, references, AI-use report)
├── 02_Verification_and_Figures/
│   ├── fig01_pareto_tradeoff.png                  # + Figure_1_Verification_Dossier_V2.4.0.pdf
│   ├── fig02_edh38_efficiency.png                 # + Figure_2_Verification_Dossier_V2.4.0.pdf
│   └── fig03_ahp_season_sensitivity.png           # + Figure_3_Verification_Dossier_V2.4.0.pdf
├── 03_Final_Technical_Package.zip                 # self-contained package (95 files); contents below
├── 04_Data_and_Results/
│   ├── AHP_Ranking_Official_Team19386_14Nov_v2.4.0.csv        # 290 Pareto-efficient configurations with AHP scores
│   ├── Sensitivity_RankTracking_LeMinh_15Nov_v1.csv           # rank tracking across the 56 one-at-a-time cases
│   └── ValidationChecks_Team19386_14Nov_v2.4.0.csv            # automated invariant checks
└── README.md

03_Final_Technical_Package.zip  →  Team19386_Final_Technical_Package_16Nov2026/
├── 01_Official_Model/                             # V2.4.0 model, parameter workbook, three weather files, official result CSVs
├── 02_Extended_Analysis/                          # exploratory: cap-120 scan, Monte Carlo, heatwave, rank tracking, figures
├── 03_Literature/                                 # contest statement, device specifications, literature register
├── 04_Audit_and_Corrections/                      # technical audit, C7 coefficient audit, corrected efficiency table
├── 05_Code_and_Reproducibility/                   # Tests/, Source_Code/, Historical_Logs/, Reproduction_Log.txt, requirements.txt
└── MANIFEST_SHA256.txt  RELEASE_STATUS.md  CHANGELOG.md  README.md  99_Archive_HASHES.txt
```

Venue identifiers in code and CSV files: `Tram_Xe_Buyt` (Bus Shelter), `Ga_Ha_Noi` (Railway Station), `San_Bay_Honolulu` (Honolulu Airport), `Sun_Devil_Stadium`.

---

## 7. Installation and reproducibility

### 7.1 Environment

Tested on **Python 3.13 (Linux)**; Windows and macOS are untested. No network access is needed. From the package folder:

```bash
cd 02_Code_and_Reproducibility/Team19386_Final_Technical_Package_16Nov2026_v1_LITE
python3 -m venv .venv && source .venv/bin/activate
pip install -r 05_Code_and_Reproducibility/requirements.txt   # numpy, pandas, matplotlib, seaborn, reportlab, openpyxl (pinned)
```

### 7.2 Reproduce the official results

```bash
python3 05_Code_and_Reproducibility/Tests/Test_Reproduction_LeMinh_16Nov_v1.py
```

Expected: **13 of 13 checks pass in about 6 seconds** (compare with `05_Code_and_Reproducibility/Reproduction_Log.txt`). To check file integrity, run `sha256sum -c MANIFEST_SHA256.txt` from the package folder.

<details>
<summary>Running model V2.4.0 manually</summary>

The model expects weather files with the suffixes `(1)`, `(1)` and `(2)`. Copy the workbook and the three weather CSVs from `01_Official_Model/Original_V2.4.0/` into one folder, rename the weather files to

```text
WeatherData_Ga_Ha_Noi_GiaPhu_08Nov_v1(1).csv
WeatherData_Honolulu_International_Airport_GiaPhu_08Nov_v1(1).csv
WeatherData_Sun_Devil_Stadium_GiaPhu_08Nov_v1(2).csv
```

and run, quoting file names in the shell:

```bash
python3 PhysicalModel_Team19386_14Nov_v2.4.0.py --data-dir <folder> \
  --workbook THONG_SO_THONG_NHAT_CHUAN_100_Team19386_14Nov_v2.4.0.xlsx --max-units 40
```
</details>

The extended-analysis scripts in `05_Code_and_Reproducibility/Source_Code/` (including `MonteCarloV24_LeMinh_15Nov_v1.py`) read the environment variables `HIMCM_MODEL`, `HIMCM_RUN_DIR`, `HIMCM_WORKBOOK` and `HIMCM_OUT`. They were run on 15 November and were not re-run for the release, except one Monte Carlo chunk of 500 draws per venue, which reproduced identically. Monte Carlo seeding is `default_rng([42, site, chunk])`, in 80 chunks of 500 draws per venue.

### 7.3 Data

Hourly weather comes from the [Open-Meteo Historical Weather API](https://open-meteo.com/en/docs/historical-weather-api): three series (Hanoi, Honolulu, Phoenix), each with 54,792 records and 11 columns from 1 June 2020 to 31 August 2026, with no duplicate timestamps, missing values or imputed rows. The model uses June to August, hourly labels 10–21 (12 h × 92 days × 7 summers = 7,728 h). These are gridded model values, not on-site measurements. Follow Open-Meteo's terms of use and attribution requirements when reusing the data.

---

## 8. Scope, limitations and caveats

* **Not calibrated or validated.** No venue was compared with field measurements. Results are conditional decision support, not certified safety predictions.
* **Reduced-order one-zone physics.** No thermal mass, spatial gradients, supply jets or occupant-level microclimates. Ventilation, sheltering, mixing and coverage coefficients are coarse.
* **Open-air stadium branch.** It sets *T*<sub>base</sub> = *T*<sub>amb</sub>, so stadium reductions depend strongly on the assumed coverage and mixing parameters.
* **Internal index.** AT<sub>mod</sub> and its 28 °C / 38 °C thresholds are an internal comparison metric, not a clinical or regulatory standard.
* **Placeholder economics.** Prices, tariffs and equipment parameters are workbook assumptions, and electricity assumes constant rated power with no part-load behaviour.
* **No ablation study.** The release contains no controlled component-ablation runs, so no ablation percentages are reported.
* **Name collision.** In the extended materials, "C7" also denotes a proposed railway-station wind-break case (`alpha_shelter` 0.35 → 0.20). It is unrelated to the configuration `C-7` (System C, seven units) and does not enter any official result.

**Release status: CONDITIONAL PASS** (16 November 2026). The official files are byte-identical to the 14 November outputs and the executed checks pass. Open items are listed in the package's `RELEASE_STATUS.md`: six of seven literature sources not read in full, unmodelled physical limits at the Bus Shelter, and no Windows test.

---

## 9. Team and workstreams

**Team #19386**, High School for the Gifted, Vietnam National University Ho Chi Minh City (VNU-HCM). Workstreams as recorded in the project's planning and parameter documents (8–15 November 2026):

* **Lê Nghiêm Gia Phú (team leader):** coordination; weather data and references; Cooling-Efficiency and Semi-Open Rail Station model formulation; Monte Carlo workstream; public blog post (with Lê Minh).
* **Lê Minh:** Python implementation of the V2.4.0 solver, optimisation, sensitivity and extended analyses, the 10,000-draw Monte Carlo script, and the audit and reproducibility package.
* **Đỗ Lê Chí Hùng:** heat-balance, geometry and cooling-system parameters, and parameter audit.
* **Trần Hưng Phát:** cost and AHP inputs, summary sheet, and operations memorandum.

## Citation

```bibtex
@misc{team19386_himcm2026,
  title        = {Cool the Heat, Not Just the Air: A Framework for Making Site-Specific
                  Decisions About Outdoor and Semi-Outdoor Public Spaces},
  author       = {Le, Nghiem Gia Phu and Le, Minh and Do, Le Chi Hung and Tran, Hung Phat},
  year         = {2026},
  note         = {HiMCM 2026, Problem B, Team \#19386. Submitted solution paper}
}
```


