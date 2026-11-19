# HiMCM_2026_ProblemB_Team19386
# Cool the Heat, Not Just the Air: Site-Specific Thermal Exposure & Cooling Optimization Framework

[![Competition](https://img.shields.io/badge/COMAP_HiMCM_2026-Problem_B-1f4e79.svg)](https://www.comap.com/contests/himcm)
[![Team](https://img.shields.io/badge/Team-%2319386-2e7d32.svg)](#team--authorship)
[![Python](https://img.shields.io/badge/Python-3.13-3776ab.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Model Version](https://img.shields.io/badge/Deterministic_Model-V2.4.0_(cap40)-informational.svg)](#methodology--physical-formulation)
[![Verification](https://img.shields.io/badge/Audit_&_Reproducibility-13%2F13_Checks_Passed-success.svg)](#installation--reproducibility)

> **A reduced-order, one-zone, quasi-steady hourly heat-and-moisture balance and multi-criteria decision framework (Pareto + AHP + 10,000-run Monte Carlo) for optimizing cooling technologies across outdoor and semi-outdoor public spaces.**

---

## 1. Overview & System Architecture

Extreme heat threatens pedestrians, commuters, and spectators in outdoor and semi-outdoor public venues where conditioned air rapidly disperses and solar radiation, natural airflow, and crowd density strongly govern human thermal stress. Relying solely on rated device cooling capacity fails to capture real occupant heat exposure.

Developed by **Team #19386** for the **2026 COMAP High School Mathematical Contest in Modeling (HiMCM), Problem B**, this repository provides the official deterministic simulation engine (`V2.4.0`), 10,000-realization stochastic Monte Carlo robustness pipeline, verification dossiers, and final technical report evaluating three cooling technologies against a no-cooling baseline (`Baseline-0`) across four contrasting venues:

- **System A:** Recirculating air-movement fans (convective cooling via elevated local air speed without dry-bulb drop or moisture addition).
- **System B:** Medium-capacity direct evaporative coolers.
- **System C:** High-capacity direct evaporative coolers.

### End-to-End Modeling & Decision Pipeline

```text
+-----------------------------------------------------------------------------------+
| 1. METEOROLOGICAL & VENUE INPUTS (Open-Meteo API, Jun–Aug 2020–2026, 09:00–21:00) |
|    • 3 Hourly Weather Series (54,792 rows each; 7,728 summer operating hrs/site)  |
|    • 4 Venues: Hanoi Bus Shelter | Hanoi Rail Station | Honolulu Airport | Phoenix|
+-----------------------------------------+-----------------------------------------+
                                          |
                                          v
+-----------------------------------------------------------------------------------+
| 2. BASELINE ONE-ZONE QUASI-STEADY HEAT & MOISTURE BALANCE (V2.4.0 Solver)         |
|    • Sol-Air Envelope Conduction (H_env) + Glazing Solar Gain + Metabolic Loads   |
|    • Coupled Wind-Driven (V_wind) & Buoyancy Stack (V_stack) Natural Ventilation  |
|    • Baseline Zone Dry-Bulb Temperature (T_base) & Humidity Ratio (w_zone)        |
+-----------------------------------------+-----------------------------------------+
                                          |
                                          v
+-----------------------------------------------------------------------------------+
| 3. COOLING TECHNOLOGY RESPONSE & APPARENT TEMPERATURE EXPOSURE                    |
|    • System A: Local air-speed amplification (v_served)                           |
|    • Systems B & C: Wet-bulb bounded supply drop (dT_supply), mixing (phi), m_evap|
|    • Team-Modified Apparent Temperature (AT_mod) for Served & Unserved Occupants  |
+-----------------------------------------+-----------------------------------------+
                                          |
                                          v
+-----------------------------------------------------------------------------------+
| 4. MULTI-CRITERIA OPTIMIZATION & ROBUSTNESS DIAGNOSTICS                           |
|    • 4-Objective Pareto Screening under cap40 (n <= 40 units; 121 configs/venue)  |
|    • Analytic Hierarchy Process (AHP) Ranking (CR = 0.003838 < 0.10)              |
|    • 56 OAT Sensitivity Tests | +3°C Heatwave Stress | 10,000-Run Monte Carlo     |
+-----------------------------------------------------------------------------------+
2. Methodology & Physical Formulation2.1 Baseline Thermal & Moisture ConservationEach venue is modeled as a representative control volume resolved at hourly time steps ($\Delta t = 1\text{ h}$) across seven summers (June–August, 2020–2026, local operating hours 09:00–21:00, yielding $7 \times 92 \times 12 = 7,728$ operating hours per site).Total envelope conductance $H_{\text{env}}$ ($\text{W/K}$) and natural infiltration volumetric flow $\dot{V}_{\text{inf}}$ ($\text{m}^3/\text{s}$) couple radiative-conductive gains with wind-driven and buoyancy-driven stack ventilation:$$H_{\text{env}} = U_{\text{roof}} A_{\text{roof}} + U_{\text{wall}} A_{\text{wall}} + U_{\text{glass}} A_{\text{glass}}$$$$\dot{V}_{\text{wind}} = C_d A_{\text{open}} \max\left(k_{\text{gio}} \, \alpha_{\text{shelter}} \, v_{\text{amb}}, \, v_{\min}\right), \qquad \dot{V}_{\text{stack}} = C_d A_{\text{open}} \sqrt{\frac{g H \vert{}T_0 - T_{\text{amb}}\vert{}}{T_{\text{amb}} + 273.15}}$$$$\dot{V}_{\text{inf}} = \sqrt{\dot{V}_{\text{wind}}^2 + \dot{V}_{\text{stack}}^2}, \qquad T_{\text{base}} = \frac{H_{\text{roof}} T_{\text{sol,roof}} + H_{\text{wall}} T_{\text{sol,wall}} + H_{\text{glass}} T_{\text{amb}} + Q_{\text{solar,glass}} + Q_{\text{internal}} + H_{\text{vent}} T_{\text{amb}}}{H_{\text{env}} + H_{\text{vent}}}$$where $H_{\text{vent}} = \rho_{\text{da}} c_{p,\text{da}} \dot{V}_{\text{inf}}$. For enclosed and semi-open structures (such as the Hanoi Railway Station), background humidity ratio $w_{\text{zone}}$ incorporates occupant latent heat generation capped at saturation $w_{\text{sat}}(T_{\text{base}})$, whereas open-air Sun Devil Stadium sets $w_{\text{zone}} = w_{\text{amb}}$.
