<div align="center">

# Cool the Heat, Not Just the Air: Site-Specific Thermal Exposure & Cooling Optimization Framework

**COMAP High School Mathematical Contest in Modeling (HiMCM) 2026 — Problem B | Team #19386**

[![Python](https://img.shields.io/badge/Python-3.13.18-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Model](https://img.shields.io/badge/Model-Deterministic%20V2.4.0%20(cap40)-0F766E?style=flat-square)](#2-methodology--physical-formulation)
[![Verification](https://img.shields.io/badge/Reproducibility%20Checks-13%2F13%20Passed-16A34A?style=flat-square)](./02_Code_and_Reproducibility/Reproduction_Log.txt)
[![AHP Consistency](https://img.shields.io/badge/AHP%20CR-0.003838%20%3C%200.10-2563EB?style=flat-square)](#23-exposure-metrics-efficiency--multi-criteria-ahp)
[![Monte Carlo](https://img.shields.io/badge/Monte%20Carlo-10%2C000%20Runs%20(seed%3D42)-7C3AED?style=flat-square)](#32-sensitivity-heatwave-stress--10000-run-monte-carlo-robustness)

*A reduced-order, one-zone, quasi-steady hourly heat-and-moisture balance and multi-criteria decision framework (4-Objective Pareto + AHP + 10,000-Run Block-Bootstrap Monte Carlo) for optimizing cooling technologies across outdoor and semi-outdoor public spaces.*

[📄 Final Report (PDF)](./01_Final_Report/19386.pdf) • [📊 AHP Rankings (CSV)](./04_Data_and_Results/AHP_Ranking_Official_Team19386_14Nov_v2.4.0.csv) • [🔍 Sensitivity Tracking (CSV)](./04_Data_and_Results/Sensitivity_RankTracking_LeMinh_15Nov_v5.0.csv) • [✅ Validation Audit (CSV)](./04_Data_and_Results/ValidationChecks_Team19386_14Nov_v2.4.0.csv)

</div>

---

## 1. Overview & System Architecture

Extreme heat threatens pedestrians, commuters, and spectators in outdoor and semi-outdoor public venues where conditioned air rapidly disperses and solar radiation, natural airflow, and crowd density strongly govern human thermal stress. Relying solely on rated device cooling capacity fails to capture real occupant heat exposure.

Developed by **Team #19386** (High School for the Gifted, VNU-HCM) for the **2026 COMAP High School Mathematical Contest in Modeling (HiMCM), Problem B**, this repository provides the official deterministic simulation engine (`V2.4.0`), 10,000-realization stochastic Monte Carlo robustness pipeline, SHA-256 verification manifests, figure audit dossiers, and final 30-page technical report evaluating three cooling technologies against a no-cooling baseline (`Baseline-0`) across four contrasting venues:

* **System A:** Recirculating air-movement fans (convective cooling via elevated local air speed without dry-bulb drop or moisture addition).
* **System B:** Medium-capacity direct evaporative coolers.
* **System C:** High-capacity direct evaporative coolers.

### End-to-End Modeling & Decision Pipeline

```text
+----------------------------------------------------------------------------+
| 1. METEOROLOGICAL & VENUE INPUTS (Open-Meteo, Jun-Aug 2020-2026, 09-21h)   |
|    * 3 Hourly Weather Series (54,792 rows each; 7,728 summer hrs/site)     |
|    * 4 Venues: Hanoi Bus Shelter | Hanoi Rail Station | Honolulu | Phoenix |
+----------------------------------------------------------------------------+
                                      |
                                      v
+----------------------------------------------------------------------------+
| 2. BASELINE ONE-ZONE QUASI-STEADY HEAT & MOISTURE BALANCE (V2.4.0 Solver)  |
|    * Sol-Air Envelope Conduction (H_env) + Glazing Solar + Metabolic Load  |
|    * Coupled Wind-Driven (V_wind) & Buoyancy Stack (V_stack) Ventilation   |
|    * Baseline Zone Dry-Bulb Temperature (T_base) & Humidity Ratio (w_zone) |
+----------------------------------------------------------------------------+
                                      |
                                      v
+----------------------------------------------------------------------------+
| 3. COOLING TECHNOLOGY RESPONSE & APPARENT TEMPERATURE EXPOSURE             |
|    * System A: Local air-speed amplification (v_served)                    |
|    * Systems B & C: Wet-bulb bounded drop (dT_supply), mixing (phi), m_ev  |
|    * Team-Modified Apparent Temperature (AT_mod) for Served & Unserved     |
+----------------------------------------------------------------------------+
                                      |
                                      v
+----------------------------------------------------------------------------+
| 4. MULTI-CRITERIA OPTIMIZATION & ROBUSTNESS DIAGNOSTICS                    |
|    * 4-Objective Pareto Screening under cap40 (n <= 40; 121 configs/site)  |
|    * Analytic Hierarchy Process (AHP) Ranking (CR = 0.003838 < 0.10)       |
|    * 56 OAT Sensitivity Tests | +3 C Heatwave Stress | 10,000-Run MC       |
+----------------------------------------------------------------------------+
```

---

## 2. Methodology & Physical Formulation

### 2.1 Baseline Thermal & Moisture Conservation

Each venue is modeled as a representative control volume resolved at hourly time steps ($\Delta t = 1\text{ h}$) across seven summers (June–August, 2020–2026, local operating hours 09:00–21:00, yielding $7 \times 92 \times 12 = 7,728$ operating hours per site). Total envelope conductance $H_{\text{env}}\ (\text{W/K})$ and natural infiltration volumetric flow $\dot{V}_{\text{inf}}\ (\text{m}^3/\text{s})$ couple radiative-conductive gains with wind-driven and buoyancy-driven stack ventilation:

$$H_{\text{env}} = U_{\text{roof}} A_{\text{roof}} + U_{\text{wall}} A_{\text{wall}} + U_{\text{glass}} A_{\text{glass}}$$

$$\dot{V}_{\text{wind}} = C_d A_{\text{open}} \max\left(k_{\text{gio}} \alpha_{\text{shelter}} v_{\text{amb}},\, v_{\text{min}}\right), \qquad \dot{V}_{\text{stack}} = C_d A_{\text{open}} \sqrt{\frac{g H \vert{}T_0 - T_{\text{amb}}\vert{}}{T_{\text{amb}} + 273.15}}$$

$$\dot{V}_{\text{inf}} = \sqrt{\dot{V}_{\text{wind}}^2 + \dot{V}_{\text{stack}}^2}, \qquad T_{\text{base}} = \frac{H_{\text{roof}} T_{\text{sol,roof}} + H_{\text{wall}} T_{\text{sol,wall}} + H_{\text{glass}} T_{\text{amb}} + Q_{\text{solar,glass}} + Q_{\text{internal}} + H_{\text{vent}} T_{\text{amb}}}{H_{\text{env}} + H_{\text{vent}}}$$

where $H_{\text{vent}} = \rho_{\text{da}} c_{p,\text{da}} \dot{V}_{\text{inf}}$. For enclosed and semi-open structures (such as the Hanoi Railway Station), background humidity ratio $w_{\text{zone}}$ incorporates occupant latent heat generation capped at saturation $w_{\text{sat}}(T_{\text{base}})$, whereas open-air Sun Devil Stadium sets $w_{\text{zone}} = w_{\text{amb}}$.

### 2.2 Cooling Technology Mechanics

* **System A (Recirculating Fans):** Elevates local air speed around served occupants ($v_{\text{served}} = \sqrt{v_{\text{background}}^2 + (v_{A,\text{person}} f_{\text{fan}})^2}$) without reducing dry-bulb temperature or altering zone moisture.
* **Systems B & C (Evaporative Coolers):** Bounded by rated device drop $\Delta T_{\text{rated}}$ and ambient wet-bulb depression $(T_{\text{amb}} - T_{\text{wb}})$:

$$\varphi = \frac{\text{mix\_eff} \cdot \dot{V}_{\text{sys}}}{\dot{V}_{\text{sys}} + f_{\text{cov}} \dot{V}_{\text{inf}}}, \qquad \Delta T_{\text{supply}} = \max\left[0,\, \min\left(\Delta T_{\text{rated}},\, \epsilon_{\text{evap}} (T_{\text{amb}} - T_{\text{wb}})\right)\right]$$

$$T_{\text{served}} = T_{\text{base}} - \varphi\,\Delta T_{\text{supply}}, \qquad \dot{m}_{\text{evap}} = \rho_{\text{da}} \dot{V}_{\text{sys}} \max(0,\, w_{\text{supply}} - w_{\text{amb}})$$

### 2.3 Exposure Metrics, Efficiency & Multi-Criteria AHP

Human thermal stress is evaluated using a team-modified apparent temperature index ($AT_{\text{mod}}$, in $^\circ\text{C}$) distinguishing served ($N_{\text{served},t}$) and unserved ($N_{\text{unserved},t}$) crowd fractions:

$$AT_{\text{mod}} = T + 0.33 e - 0.70 v + \frac{0.70\, Q_g}{v + 10} - 4.25, \qquad Q_g = 0.1333\, f_{\text{sun}} I_{\text{sol}}$$

with Comfort Threshold $AT_{\text{mod}} \le 28^\circ\text{C}$ and Danger Threshold $AT_{\text{mod}} \ge 38^\circ\text{C}$. Seven-summer mean occupancy-weighted Excess Degree-Hours above $38^\circ\text{C}$ ($\text{EDH}_{38}$, in $\text{person-}^\circ\text{C}\cdot\text{h}$), Danger Share ($S_{\text{danger}}$), Cooling Electricity Efficiency ($\eta_{\text{EDH}_{38}}$), and Equivalent Annualized Cost ($C_{\text{ann}}$) are computed as:

$$\text{EDH}_{38} = \frac{1}{7} \sum_{y=2020}^{2026} \sum_{t \in y} \left[ N_{\text{served},t} \max(0,\, AT_{\text{served},t} - 38) + N_{\text{unserved},t} \max(0,\, AT_{\text{unserved},t} - 38) \right] \Delta t$$

$$\eta_{\text{EDH}_{38}} = \frac{\text{EDH}_{38,\text{baseline}} - \text{EDH}_{38,\text{config}}}{E_{\text{kWh}}} \quad (\text{undefined } 0/0 \text{ when } \text{EDH}_{38,\text{baseline}} = 0)$$

$$C_{\text{ann}} = n (C_{\text{capex}} + C_{\text{install}}) \left( \frac{d(1+d)^L}{(1+d)^L - 1} + r_{\text{maint}} \right) + E_{\text{kWh}} c_e + W_L c_w$$

Candidates within `cap40` ($n \le 40$, 121 configurations per venue) are filtered via 4-objective Pareto dominance and ranked via Analytic Hierarchy Process (AHP) weights ($\text{CR} = 0.003838 < 0.10$):

1. **Thermal Heat-Exposure Reduction ($\text{EDH}_{38}$ & $S_{\text{danger}}$):** `0.4236`
2. **Equivalent Annualized Cost ($C_{\text{ann}}$):** `0.2270`
3. **Seasonal Electricity & Water Demand:** `0.2270`
4. **Operational & Spatial Deployment Feasibility:** `0.1223`

---

## 3. Key Results & Benchmarks

### 3.1 Official Deterministic V2.4.0 Recommendations (`cap40`)

| Venue (Geometry & Climate) | AHP Pick | AHP Score | Base $S_{\text{danger}}$ | New $S_{\text{danger}}$ | $\text{EDH}_{38}\ \Delta$ | $\eta_{\text{EDH}_{38}}$ | EAC (USD/yr) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Sun Devil Stadium, Phoenix**<br>($100 \times 25 \times 10\text{ m}$, open, hot-arid) | **C-40** | 0.6507 | 66.31% | 50.25% | -30.03% | 15.591 | 60,948.07 |
| **Hanoi Railway Station**<br>($40 \times 15 \times 6\text{ m}$, semi-open, humid) | **B-40** | 0.8101 | 56.84% | 49.50% | -27.02% | 0.8200 | 25,334.40 |
| **Hanoi Bus Shelter**<br>($\approx 6.4\ \text{m}^2$ roof, hot-humid) | **C-7\*** | 0.8590 | 45.22% | 35.27% | -43.63% | 0.1846 | 10,270.69 |
| **Honolulu Int'l Airport**<br>($60 \times 20 \times 7\text{ m}$, enclosed/shaded) | **Base-0** | 0.7730 | 0.00% | 0.00% | 0.00% | Undef (0/0) | 0.00 |

> **\*Operational Guidance:** Priority #1 is phased deployment of **C-40** at Sun Devil Stadium paired with shade and hydration protocols; phased **B-40** deployment in high-density platform waiting zones at Hanoi Railway Station; a **1–2-unit instrumented pilot** at Hanoi Bus Shelter due to the $6.4\ \text{m}^2$ roof footprint constraint; and retaining **Baseline-0** at Honolulu International Airport.

### 3.2 Sensitivity, Heatwave Stress & 10,000-Run Monte Carlo Robustness

| Venue & Selected Config | OAT Rank Shifts (of 56) | $+3^\circ\text{C}$ Heatwave Base $S_{\text{danger}}$ | MC Base Mean $S_{\text{danger}}$ | MC Config Mean $S_{\text{danger}}$ |
| :--- | :---: | :---: | :---: | :---: |
| **Hanoi Bus Shelter (C-7)** | 13/56 | 81.54% | 46.11% | 36.06% |
| **Hanoi Railway Station (B-40)** | 5/56 | 88.38% | 56.91% | 49.63% |
| **Honolulu Airport (Baseline-0)** | 2/56 | 0.95% | 0.00% | 0.00% |
| **Sun Devil Stadium (C-40)** | 0/56 | 84.79% | 66.42% | 50.41% |

* **Monte Carlo Design:** 10,000 realizations per configuration across 16 selected configurations (deterministic seed `42` via `default_rng([42, site, chunk])`, 80 chunks × 500 draws), combining 7-summer year-block bootstrap resampling (2020–2026) with independent $U(0.8, 1.2)$ perturbations across 10 parameter groups and computational upper-bound clipping at `1.0` for `alpha_shelter`, `eps_evap`, `mix_eff`, and `f_dung`.
* **Verification Integrity:** `13 / 13` automated reproducibility checks pass in ~5.7s (`0.0e+00` maximum absolute numerical deviation across all 484 configuration rows).

---

## 4. Visual Benchmarks & Verification Dossiers

Each primary figure is accompanied by a standalone PDF audit dossier in `./03_Verification_and_Figures/` documenting underlying numerical tables and invariance checks:

### Figure 1: Four-Objective `cap40` Pareto Frontiers ([Verification Dossier PDF](./03_Verification_and_Figures/Figure_1_Verification_Dossier_V2.4.0.pdf))
![Four-Objective cap40 Pareto Frontiers](./03_Verification_and_Figures/fig01_pareto_tradeoff.png)

### Figure 2: $\text{EDH}_{38}$ Electricity Efficiency & Occupancy-Weighted $AT_{\text{mod}}$ Distributions ([Verification Dossier PDF](./03_Verification_and_Figures/Figure_2_Verification_Dossier_V2.4.0.pdf))
![EDH38 Efficiency and AT_mod Distributions](./03_Verification_and_Figures/fig02_edh38_efficiency.png)

### Figure 3: AHP Weight Perturbation ($\pm 20\%$) & Operating-Season Sensitivity ([Verification Dossier PDF](./03_Verification_and_Figures/Figure_3_Verification_Dossier_V2.4.0.pdf))
![AHP Weight and Operating-Season Sensitivity](./03_Verification_and_Figures/fig03_ahp_season_sensitivity.png)

---

## 5. Repository Structure

```text
HiMCM_2026_ProblemB_19386/
├── 01_Final_Report/
│   └── 19386.pdf                                              # Official 30-page COMAP HiMCM submission
├── 02_Code_and_Reproducibility/
│   ├── PhysicalModel_LeMinh_13Nov_v1.0.py                   # Deterministic hourly solver, Pareto, AHP & OAT
│   ├── MonteCarloV24_GiaPhu_15Nov_v1.py                     # Parallel 10,000-run Monte Carlo engine (seed=42)
│   ├── MANIFEST_SHA256.txt                                  # SHA-256 checksums for inputs, code & tables
│   ├── Reproduction_Log.txt                                 # Python 3.13.18 verification log (13/13 passed)
│   └── Team19386_Final_Technical_Package_16Nov2026_v1_LITE/ # Weather CSVs, parameter workbook & test harness
├── 03_Verification_and_Figures/
│   ├── fig01_pareto_tradeoff.png                            # Fig 1: Four-objective cap40 Pareto frontiers
│   ├── Figure_1_Verification_Dossier_V2.4.0.pdf             # Audit verification dossier for Figure 1
│   ├── fig02_edh38_efficiency.png                           # Fig 2: EDH38 efficiency & AT_mod histograms
│   ├── Figure_2_Verification_Dossier_V2.4.0.pdf             # Audit verification dossier for Figure 2
│   ├── fig03_ahp_season_sensitivity.png                     # Fig 3: ±20% AHP weights & season scaling
│   └── Figure_3_Verification_Dossier_V2.4.0.pdf             # Audit verification dossier for Figure 3
└── 04_Data_and_Results/
    ├── AHP_Ranking_Official_Team19386_14Nov_v2.4.0.csv      # Official V2.4.0 AHP scores & rankings (484 rows)
    ├── Sensitivity_RankTracking_LeMinh_15Nov_v5.0.csv       # Rank reversal tracking across 56 OAT cases
    └── ValidationChecks_Team19386_14Nov_v2.4.0.csv          # Automated audit & invariant checks (10/10 True)
```

---

## 6. Installation & Reproducibility

### 6.1 Environment Setup

Tested on **Python 3.13.18 (Linux)**. Install the required scientific computing stack:

```bash
pip install numpy==2.5.3 pandas==3.0.5 matplotlib==3.11.2 seaborn==0.13.2 reportlab==5.0.1 openpyxl==3.1.5
```

### 6.2 Running the Deterministic & Monte Carlo Pipelines

1. **Deterministic V2.4.0 Solver (`cap40` sweep, Pareto screening, AHP & OAT sensitivity):**
   Ensure the standardized parameter workbook (`THONG_SO_THONG_NHAT_CHUAN_100_v2.3.2 (1).xlsx`) and the three Open-Meteo hourly CSV files (`Weather Data_Ga_Ha_Noi_GiaPhu_08Nov_v1(1).csv`, `WeatherData_Honolulu_International_Airport_GiaPhu_08Nov_v1(1).csv`, `WeatherData_Sun_Devil_Stadium_GiaPhu_08Nov_v1(2).csv`) are available inside `<DATA_DIR>` (bundled inside `02_Code_and_Reproducibility/Team19386_Final_Technical_Package_16Nov2026_v1_LITE/`), then run:

   ```bash
   python3 02_Code_and_Reproducibility/PhysicalModel_LeMinh_13Nov_v1.0.py \
     --data-dir <DATA_DIR> \
     --workbook "THONG_SO_THONG_NHAT_CHUAN_100_v2.3.2 (1).xlsx" \
     --max-units 40
   ```

2. **10,000-Run Stochastic Monte Carlo Simulation (`seed = 42`):**
   Ensure the extended model wrapper (`ModelExt_LeMinh_15Nov_v1.py`) and official V2.4.0 context variables (`HIMCM_MODEL`, `HIMCM_RUN_DIR`, `HIMCM_WORKBOOK`, `HIMCM_OUT`) are configured, then execute:

   ```bash
   export HIMCM_OUT=./mc_outputs
   python3 02_Code_and_Reproducibility/MonteCarloV24_GiaPhu_15Nov_v1.py
   ```

3. **SHA-256 Checksum Verification:**

   ```bash
   cd 02_Code_and_Reproducibility && sha256sum -c MANIFEST_SHA256.txt
   ```

---

## 7. Scope, Limitations & Caveats

* **Reduced-Order One-Zone Approximation:** Resolves venue-level hourly quasi-steady energy and moisture balances rather than 3D computational fluid dynamics (CFD), localized jet plumes, or train-induced piston airflow in the railway station.
* **Internal Comparative Index ($AT_{\text{mod}}$):** $AT_{\text{mod}}$ is an internal comparative exposure metric adapted from Steadman's apparent temperature formulation and is not interchangeable with clinical UTCI, NWS Heat Index, or field-calibrated physiological measurements.
* **Uncalibrated Site & Economic Assumptions:** While deterministic software reproducibility is verified (`13/13` checks), none of the four venues has been calibrated against multi-sensor field trials, and equipment costs rely on standardized workbook assumptions rather than binding vendor quotes.
* **Computational Upper Clipping in Monte Carlo:** Perturbed draws for `alpha_shelter`, `eps_evap`, `mix_eff`, and `f_dung` are clipped at `1.0`, introducing boundary skewness in upper-tail stochastic distributions.

---

## 8. Team & Authorship

**COMAP HiMCM 2026 Team #19386** — High School for the Gifted, Vietnam National University Ho Chi Minh City (VNU-HCM):

* **Lê Nghiêm Gia Phú (Team Leader):** Project coordination and task delegation; mathematical formulation of the core Cooling-Efficiency and Semi-Open Rail Station models; design and execution of the 10,000-scenario Monte Carlo simulation across 16 configurations (`MonteCarloV24_GiaPhu_15Nov_v1.py`); authorship of the public-facing Blog Post.
* **Lê Minh:** Python implementation of the deterministic V2.4.0 solver (`PhysicalModel_LeMinh_13Nov_v1.0.py`), extended `cap120` sweeps, OAT sensitivity scripts, and SHA-256 verification tables.
* **Đỗ Lê Chí Hùng:** Literature synthesis, equipment specification verification, and report drafting.
* **Trần Hưng Phát:** Data preprocessing, occupancy profile parameterization, and executive memo/guidelines drafting.

## Citation

```bibtex
@techreport{team19386_himcm2026,
  title       = {Cool the Heat, Not Just the Air: A Framework for Making Site-Specific Decisions About Outdoor and Semi-Outdoor Public Spaces},
  author      = {Le, Nghiem Gia Phu and Le, Minh and Do, Le Chi Hung and Tran, Hung Phat},
  institution = {Consortium for Mathematics and Its Applications (COMAP) High School Mathematical Contest in Modeling (HiMCM), Team \#19386},
  year        = {2026},
  month       = {November},
  note        = {Problem B Final Solution Paper}
}
```
