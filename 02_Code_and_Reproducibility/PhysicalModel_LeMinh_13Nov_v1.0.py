#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HiMCM 2026 Problem B — Team 19386
PhysicalModel_LeMinh_13Nov_v2.3.3.py

Rebuilt from the V2.3.2 parameter workbook and the team's V2.1 heat-balance note.
Scope: deterministic hourly model, moisture balance, A/B/C spot-cooling comparison,
Pareto screening, provisional AHP and one-at-a-time sensitivity. Monte Carlo is
intentionally excluded (Gia Phu's workstream).

IMPORTANT MODEL STATUS
- This is a reduced-order, quasi-steady model, not CFD or a validated building model.
- The outdoor stadium branch sets the background air temperature to weather dry-bulb;
  solar load enters the team's modified apparent-temperature indicator, not zone air.
- System A fan flow is treated as recirculating air and is NOT counted as outdoor
  dry-air ventilation or moisture dilution.
- AHP pairwise judgments in V2.3.2 are marked DE_XUAT_CHO_NHOM; scores are provisional.
- DNI is not fabricated from GHI. Direct-horizontal radiation is converted to DNI only
  when solar zenith is below the configured 80-degree cutoff; low-sun/night values are
  flagged and excluded from DNI-based calculations. Current comfort calculation uses
  GHI as specified in the workbook's team-modified Qg assumption.
"""
from __future__ import annotations
import argparse, json, math, os, warnings
from pathlib import Path
import numpy as np
import pandas as pd

VENUES = {
    "Ga_Ha_Noi": {"csv": "WeatherData_Ga_Ha_Noi_GiaPhu_08Nov_v1(1).csv", "workbook_col": "Ga_Ha_Noi", "lat": 21.0285, "lon": 105.8542},
    "San_Bay_Honolulu": {"csv": "WeatherData_Honolulu_International_Airport_GiaPhu_08Nov_v1(1).csv", "workbook_col": "San_Bay_Honolulu", "lat": 21.3187, "lon": -157.9225},
    "Sun_Devil_Stadium": {"csv": "WeatherData_Sun_Devil_Stadium_GiaPhu_08Nov_v1(2).csv", "workbook_col": "Sun_Devil_Stadium", "lat": 33.4255, "lon": -111.9400},
}
VENUE_ORDER = list(VENUES)
R_D, R_V, CP_DA, H_FG, P0 = 287.058, 461.495, 1005.0, 2.45e6, 101325.0
G = 9.80665


def saturation_vapor_pressure_pa(Tc):
    # Buck equation over liquid water; Pa
    return 611.21 * np.exp((18.678 - Tc / 234.5) * (Tc / (257.14 + Tc)))


def psychrometrics(Tc, rh_pct, altitude_m):
    Tk = np.asarray(Tc, float) + 273.15
    rh = np.clip(np.asarray(rh_pct, float) / 100.0, 0.0, 1.0)
    P = P0 * (1.0 - 2.25577e-5 * float(altitude_m)) ** 5.25588
    pv = rh * saturation_vapor_pressure_pa(np.asarray(Tc, float))
    pd = np.maximum(P - pv, 100.0)
    w = 0.62198 * pv / pd                         # kg water/kg dry air
    rho_da = pd / (R_D * Tk)                     # kg dry air/m3 mixture
    rho_moist = rho_da * (1.0 + w)               # kg moist air/m3
    cp_moist = CP_DA + 1860.0 * w                # J/(kg dry air K), approx mixture basis
    return P, pv, pd, w, rho_da, rho_moist, cp_moist


def wet_bulb_stull(Tc, rh_pct, pressure_pa=None):
    """Stull (2011) for RH>=5%; psychrometric root fallback below 5% RH."""
    T=np.asarray(Tc,float); RH=np.asarray(rh_pct,float)
    RHc=np.clip(RH,5.0,99.0)
    tw=(T*np.arctan(0.151977*np.sqrt(RHc+8.313659)) + np.arctan(T+RHc)
        - np.arctan(RHc-1.676331) + 0.00391838*RHc**1.5*np.arctan(0.023101*RHc) - 4.686035)
    low=RH<5.0
    if np.any(low):
        # Solve e = es(Tw) - gamma*(T-Tw) with a standard psychrometric coefficient.
        P=float(pressure_pa if pressure_pa is not None else P0)
        pv=(RH[low]/100.0)*saturation_vapor_pressure_pa(T[low])
        tlo=T[low]-80.0; thi=T[low]
        for _ in range(55):
            mid=(tlo+thi)/2
            gamma=0.00066*(1+0.00115*mid)*(P/1000.0)*1000.0 # Pa/K
            f=saturation_vapor_pressure_pa(mid)-gamma*(T[low]-mid)-pv
            above=f>0
            thi=np.where(above,mid,thi); tlo=np.where(above,tlo,mid)
        tw[low]=(tlo+thi)/2
    return tw


def solar_zenith_deg(timestamp_local, lat, lon):
    """NOAA fractional-year approximation; timestamps must already be local standard time."""
    ts = pd.DatetimeIndex(timestamp_local)
    day = ts.dayofyear.to_numpy(); hour = ts.hour.to_numpy() + ts.minute.to_numpy()/60
    gamma = 2*np.pi/365.0*(day-1+(hour-12)/24)
    eqtime = 229.18*(0.000075+0.001868*np.cos(gamma)-0.032077*np.sin(gamma)
                     -0.014615*np.cos(2*gamma)-0.040849*np.sin(2*gamma))
    decl = (0.006918-0.399912*np.cos(gamma)+0.070257*np.sin(gamma)
            -0.006758*np.cos(2*gamma)+0.000907*np.sin(2*gamma)
            -0.002697*np.cos(3*gamma)+0.00148*np.sin(3*gamma))
    # longitude positive east; local standard-time meridian is approximated from timezone offset
    # Use apparent solar time correction relative to UTC offset encoded in timestamp-naive local data.
    # For site-level low-sun flagging, longitude + local standard meridian is supplied separately below.
    # Local standard meridian estimates: Hanoi +105, Honolulu -150, Phoenix -105.
    lsm = 105 if lon > 0 else (-150 if lon < -130 else -105)
    time_offset = eqtime + 4*(lon-lsm)
    tst = hour*60 + time_offset
    ha = np.deg2rad(tst/4 - 180)
    latr = np.deg2rad(lat)
    cosz = np.sin(latr)*np.sin(decl)+np.cos(latr)*np.cos(decl)*np.cos(ha)
    return np.rad2deg(np.arccos(np.clip(cosz,-1,1)))


def load_workbook_parameters(path):
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True, read_only=False)
    def table(sheet_name):
        ws=wb[sheet_name]
        header_row=None
        for r in range(1,ws.max_row+1):
            if ws.cell(r,1).value == "Ky_hieu": header_row=r; break
        if header_row is None: raise ValueError(f"Không tìm thấy header Ky_hieu ở {sheet_name}")
        headers=[ws.cell(header_row,c).value for c in range(1,ws.max_column+1)]
        out={}
        for r in range(header_row+1,ws.max_row+1):
            key=ws.cell(r,1).value
            if not key: continue
            out[str(key)]={str(headers[c-1]):ws.cell(r,c).value for c in range(1,len(headers)+1) if headers[c-1]}
        return out
    t1=table("Cân bằng nhiệt trạm"); t2=table("Hệ làm mát, tiêu chí, phân phối"); t3=table("Chi phí + AHP"); t4=table("Dữ liệu + Mô phỏng")
    return {"geometry":t1,"cooling":t2,"cost":t3,"data":t4}


def val(tables, group, key, venue_col, default=None):
    x=tables[group].get(key,{}).get(venue_col, default)
    if x is None: return default
    return x


def get_venue_params(tables, venue_col):
    p={}
    for key in ["L","W","H","A_mai","A_kinh","A_mo","U_tuong","alpha_tuong","alpha_shelter","U_mai","alpha_abs","h_o","rad_cool_K","U_kinh","SHGC","f_dung","C_d","k_gio","v_min","n_min","n_mode","n_max","q_s","q_lat","P_tb","T_set","RH_set"]:
        p[key]=val(tables,"geometry",key,venue_col)
    for key in ["P_A","P_B","P_C","mix_eff","cluster","A_cover_A","A_cover_B","A_cover_C","v_A_person","T_fan_full","T_fan_zero","V_sys_A","V_sys_B","V_sys_C","eps_evap","theta_comf_AT","theta_danger_AT","dT_rated_A","dT_rated_B","dT_rated_C","op_start","op_end"]:
        p[key]=val(tables,"cooling",key,venue_col)
    for key in ["n_units","days_year","elec_price","water_price","capex_A","capex_B","capex_C","install_A","install_B","install_C","lifetime_yr","maint_rate","discount_rate","co2_factor"]:
        p[key]=val(tables,"cost",key,venue_col)
    for key in ["cao_do_m","n_units"]: p[key]=val(tables,"data",key,venue_col)
    # String occupancy profile is in geometry table row 49.
    raw=val(tables,"geometry","occ_profile_h10_21",venue_col,"")
    p["occ_profile"]=np.array([float(z.strip().replace(",",".")) for z in str(raw).split(";") if z.strip()],float)
    if p["cao_do_m"] is None: p["cao_do_m"]={"Ga_Ha_Noi":10,"San_Bay_Honolulu":0,"Sun_Devil_Stadium":356}.get(venue_col,0)
    p["A_floor"]=float(p["L"])*float(p["W"])
    p["occupancy_scale"]=1.0
    if p.get("A_tuong") is None or pd.isna(p.get("A_tuong")):
        p["A_tuong"]=max(0.0,2*(float(p["L"])+float(p["W"]))*float(p["H"])-2*float(p["A_mo"])-float(p["A_kinh"]))
    p["CRF"]=float(p["discount_rate"])*(1+float(p["discount_rate"]))**int(p["lifetime_yr"])/((1+float(p["discount_rate"]))**int(p["lifetime_yr"])-1)
    return p


def load_weather(path, site):
    df=pd.read_csv(path,parse_dates=["Timestamp"])
    required=["Timestamp","Temp_C","Relative_Humidity_Pct","Solar_Radiation_Wm2","Wind_Speed_ms","Cloud_Cover_Pct","Direct_Horizontal_Wm2","Imputed"]
    missing=[c for c in required if c not in df.columns]
    if missing: raise ValueError(f"{path}: thiếu cột {missing}")
    if df.Timestamp.duplicated().any(): raise ValueError(f"{path}: timestamp trùng")
    if df[required].isna().any().any(): raise ValueError(f"{path}: NaN ở cột bắt buộc")
    # Dataset convention in workbook: local timestamps, June-August, labels 10..21 = intervals 09:00-21:00.
    mask=df.Timestamp.dt.month.isin([6,7,8]) & df.Timestamp.dt.hour.between(10,21)
    out=df.loc[mask].copy().sort_values("Timestamp").reset_index(drop=True)
    if len(out)!=7728: warnings.warn(f"{site}: expected 7728 operating rows, got {len(out)}")
    if (out.Imputed.astype(float)!=0).any(): warnings.warn(f"{site}: operating rows include imputed data")
    if not out.Relative_Humidity_Pct.between(0,100).all(): raise ValueError(f"{site}: RH outside [0,100]")
    if (out.Solar_Radiation_Wm2<0).any() or (out.Wind_Speed_ms<0).any(): raise ValueError(f"{site}: negative solar/wind")
    meta=VENUES[site]
    zen=solar_zenith_deg(out.Timestamp,meta["lat"],meta["lon"])
    cosz=np.cos(np.deg2rad(zen))
    # Direct-horizontal irradiance is not DNI. Convert only above sun elevation 10 degrees.
    valid=(zen<80.0)&(cosz>0)
    out["Solar_Zenith_Deg"]=zen
    out["DNI_Est_Wm2"]=np.where(valid,out.Direct_Horizontal_Wm2.to_numpy(float)/np.maximum(cosz,np.cos(np.deg2rad(80))),np.nan)
    out["DNI_Valid_Flag"]=valid.astype(int)
    z=float({"Ga_Ha_Noi":10,"San_Bay_Honolulu":0,"Sun_Devil_Stadium":356}[meta["workbook_col"]])
    P= P0*(1.0-2.25577e-5*z)**5.25588
    out["WetBulb_C"]=wet_bulb_stull(out.Temp_C.to_numpy(float),out.Relative_Humidity_Pct.to_numpy(float),P)
    return out


def occupancy_series(df,p):
    prof=np.asarray(p["occ_profile"],float)
    if len(prof)!=12: raise ValueError(f"Occupancy profile must have 12 values, got {len(prof)}")
    prof=prof/prof.mean()
    # CSV labels 10..21 map directly to the 12 workbook profile slots.
    slot=df.Timestamp.dt.hour.to_numpy()-10
    base=(float(p["n_min"])+float(p["n_mode"])+float(p["n_max"]))/3.0
    uncapped=base*prof[slot]*float(p.get("occupancy_scale",1.0))
    return np.minimum(float(p["n_max"]),uncapped),uncapped


def ventilation_flow(df,p,Tbase_guess):
    wind=np.maximum(float(p["k_gio"])*float(p["alpha_shelter"])*df.Wind_Speed_ms.to_numpy(float),float(p["v_min"]))
    cd=float(p["C_d"]); area=float(p["A_mo"]); H=float(p["H"])
    vwind=cd*area*wind
    tk=df.Temp_C.to_numpy(float)+273.15
    vstack=cd*area*np.sqrt(G*H*np.abs(np.asarray(Tbase_guess)-df.Temp_C.to_numpy(float))/np.maximum(tk,1))
    return np.sqrt(vwind**2+vstack**2),vwind,vstack


def solve_base_state(df,p,site):
    T=df.Temp_C.to_numpy(float); ghi=df.Solar_Radiation_Wm2.to_numpy(float)
    cloud=np.clip(df.Cloud_Cover_Pct.to_numpy(float)/100,0,1)
    N,_=occupancy_series(df,p)
    # Outdoor stadium: no enclosed zone envelope/air volume; use ambient dry-bulb as background air.
    if site=="Sun_Devil_Stadium" or (float(p["A_floor"])>0 and float(p["A_mai"])+float(p["A_kinh"])+float(p["A_tuong"])==0):
        return T.copy(),np.zeros_like(T),np.zeros_like(T),np.zeros_like(T),np.zeros_like(T)
    P,pv,pd,w,rho_da,rho_m,rho_cp=psychrometrics(T,df.Relative_Humidity_Pct.to_numpy(float),p["cao_do_m"])
    cp_da=CP_DA
    # Sol-air approximation. Vertical irradiance fraction remains a team assumption (f_dung).
    Tsol_roof=T+float(p["alpha_abs"])*ghi/float(p["h_o"])-float(p["rad_cool_K"])*(1-0.8*cloud)
    Tsol_wall=T+float(p["alpha_tuong"])*float(p["f_dung"])*ghi/float(p["h_o"])
    Hroof=float(p["U_mai"])*float(p["A_mai"])
    Hwall=float(p["U_tuong"])*float(p["A_tuong"])
    Hglass=float(p["U_kinh"])*float(p["A_kinh"])
    Henv=Hroof+Hwall+Hglass
    Qsolar_glass=float(p["A_kinh"])*float(p["SHGC"])*float(p["f_dung"])*ghi
    Qinternal=N*float(p["q_s"])+float(p["P_tb"])
    # First estimate and stack ventilation; then resolve once (quasi-steady reduced-order).
    V0=np.maximum(float(p["C_d"])*float(p["A_mo"])*np.maximum(float(p["k_gio"])*float(p["alpha_shelter"])*df.Wind_Speed_ms.to_numpy(float),float(p["v_min"])),1e-6)
    Hvent=rho_da*cp_da*V0
    numerator=Hroof*Tsol_roof+Hwall*Tsol_wall+Hglass*T+Qsolar_glass+Qinternal+Hvent*T
    T0=numerator/np.maximum(Henv+Hvent,1e-9)
    Vtot,_,_=ventilation_flow(df,p,T0)
    Hvent=rho_da*cp_da*Vtot
    Tbase=(Hroof*Tsol_roof+Hwall*Tsol_wall+Hglass*T+Qsolar_glass+Qinternal+Hvent*T)/np.maximum(Henv+Hvent,1e-9)
    return Tbase,Vtot,Henv,Hvent,Qinternal


def apparent_temperature(T,rh,wind,ghi,f_sun):
    # Team-modified Steadman-style index. Not the unmodified Steadman or UTCI.
    pv= (rh/100.0)*saturation_vapor_pressure_pa(T)
    e_hpa=pv/100.0
    Qg=0.1333333333*np.asarray(f_sun)*np.asarray(ghi)
    return T+0.348*e_hpa-0.70*np.asarray(wind)+0.70*Qg/(np.asarray(wind)+10.0)-4.25


def annualize_mean_season(values, timestamps):
    frame=pd.DataFrame({"year":timestamps.dt.year.to_numpy(),"value":np.asarray(values,float)})
    byyear=frame.groupby("year").value.sum()
    return float(byyear.mean())


def simulate_configuration(df,p,site,system,n_units):
    T_out=df.Temp_C.to_numpy(float); RH=df.Relative_Humidity_Pct.to_numpy(float)
    GHI=df.Solar_Radiation_Wm2.to_numpy(float); wind=df.Wind_Speed_ms.to_numpy(float)
    cloud=np.clip(df.Cloud_Cover_Pct.to_numpy(float)/100,0,1)
    N, N_uncapped=occupancy_series(df,p)
    Tbase,Vinf,Henv,Hvent,Qint=solve_base_state(df,p,site)
    P,pv,pd,wout,rho_da,rho_moist,cp_moist=psychrometrics(T_out,RH,p["cao_do_m"])
    Tw=df.WetBulb_C.to_numpy(float) if "WetBulb_C" in df.columns else wet_bulb_stull(T_out,RH,P)
    floor=float(p["A_floor"])
    Acover=float(p.get(f"A_cover_{system}",0) or 0) if system else 0.0
    f_cov=min(1.0,n_units*Acover/max(floor,1e-9)) if system else 0.0
    f_people=min(1.0,float(p["cluster"])*f_cov) if system else 0.0
    # Convert 10-m weather wind to a representative occupied-zone background speed.
    v_background=float(p["k_gio"])*float(p["alpha_shelter"])*wind
    Nserved=N*f_people; Nunserved=N-Nserved
    # Shared background-zone moisture balance for the full occupant load.
    # Recirculating fans do not increase outdoor-air mass flow.
    if site=="Sun_Devil_Stadium":
        w_zone=wout.copy(); RH_zone=RH.copy()
    else:
        mdot_da_background=rho_da*np.maximum(Vinf,0.05)
        w_zone_raw=wout+N*float(p["q_lat"])/H_FG/mdot_da_background
        wsat_zone=0.62198*saturation_vapor_pressure_pa(Tbase)/np.maximum(P-saturation_vapor_pressure_pa(Tbase),100.0)
        w_zone=np.minimum(w_zone_raw,wsat_zone)
        RH_zone=100*np.minimum(1.0,(w_zone*P/(0.62198+w_zone))/saturation_vapor_pressure_pa(Tbase))
    mix=float(p["mix_eff"])
    Vsys=float(p.get(f"V_sys_{system}",0) or 0)*n_units if system else 0.0
    phi=mix*Vsys/np.maximum(Vsys+f_cov*Vinf,1e-9) if system in ("B","C") else 0.0
    Tserved=Tbase.copy(); RHserved=RH.copy(); wserved=wout.copy(); vserved=wind.copy()
    water_rate=np.zeros_like(T_out)
    if system=="A":
        # Recirculating fan: no air-temperature reduction, no outdoor-air dilution.
        vA=float(p["v_A_person"])
        ambient_local=v_background
        fan_factor=np.clip((float(p["T_fan_zero"])-Tbase)/max(float(p["T_fan_zero"])-float(p["T_fan_full"]),1e-9),0,1)
        v_increment=vA*fan_factor
        vserved=np.sqrt(ambient_local**2+v_increment**2)
        # All occupants share the same background-zone humidity; fan flow is not fresh air.
        wserved=w_zone.copy()
        RHserved=RH_zone.copy()
    elif system in ("B","C"):
        rated=float(p[f"dT_rated_{system}"])
        dT_supply=np.maximum(0,np.minimum(rated,float(p["eps_evap"])*(T_out-Tw)))
        Tserved=Tbase-phi*dT_supply
        # Approximate adiabatic humidity rise on dry-air basis; use actual supply mass flow.
        delta_w_supply=CP_DA*dT_supply/H_FG
        mdot_evap=rho_da*Vsys*delta_w_supply
        water_rate=mdot_evap # kg/s; rho_water=1 kg/L
        # Moisture balance includes evaporated water + latent moisture from served people.
        # Device supply flow is included for B/C; System A fan flow is deliberately excluded.
        mdot_da_zone=rho_da*np.maximum(Vsys+f_cov*Vinf,0.05)
        latent_people=Nserved*float(p["q_lat"])/H_FG
        wserved_raw=wout+(mdot_evap+latent_people)/mdot_da_zone
        wsat_local=0.62198*saturation_vapor_pressure_pa(Tserved)/np.maximum(P-saturation_vapor_pressure_pa(Tserved),100.0)
        wserved=np.minimum(wserved_raw,wsat_local)
        RHserved=100*np.minimum(1.0,(wserved*P/(0.62198+wserved))/saturation_vapor_pressure_pa(Tserved))
        vserved=v_background
    # For outdoor stadium, avoid treating occupant moisture as a confined-zone humidity rise.
    if site=="Sun_Devil_Stadium":
        RHserved=RH.copy(); wserved=wout.copy()
    # Team sun-exposure factor: fully roofless stadium exposed; fully roofed indoor zones shaded.
    f_sun=1.0 if float(p["A_mai"])<=0 else max(0.0,1.0-float(p["A_mai"])/max(floor,1e-9))
    AT_unserved=apparent_temperature(Tbase,RH_zone,v_background,GHI,f_sun)
    AT_served=apparent_temperature(Tserved,RHserved,vserved,GHI,f_sun)
    if system=="A":
        # Fan's direct comfort benefit is represented by increased local air speed only.
        # No artificial subtraction of 3 C; the 3 C is a problem-statement upper bound.
        pass
    AT_people=np.where(Nserved>0,AT_served,AT_unserved)
    total_person_hours=float(N.sum())
    comfort_ph=float((Nserved*(AT_served<=float(p["theta_comf_AT"]))+Nunserved*(AT_unserved<=float(p["theta_comf_AT"]))).sum())
    danger_ph=float((Nserved*(AT_served>=float(p["theta_danger_AT"]))+Nunserved*(AT_unserved>=float(p["theta_danger_AT"]))).sum())
    stress_ph=max(0.0,total_person_hours-comfort_ph-danger_ph)
    edh=float((Nserved*np.maximum(0,AT_served-float(p["theta_danger_AT"]))+Nunserved*np.maximum(0,AT_unserved-float(p["theta_danger_AT"]))).sum()/7.0)
    comfort_share=comfort_ph/max(total_person_hours,1e-12)
    danger_share=danger_ph/max(total_person_hours,1e-12)
    # Electrical energy: rated power for every scheduled operating hour; duty-cycle not provided.
    Psys=float(p.get(f"P_{system}",0) or 0) if system else 0.0
    hours_per_season=int(p["days_year"])*12
    E_kwh=Psys*n_units*hours_per_season
    water_L=annualize_mean_season(water_rate*3600.0,timestamps=df.Timestamp)/float(p.get("rho_water",1) or 1) if system in ("B","C") else 0.0
    capex=float(p.get(f"capex_{system}",0) or 0) if system else 0.0
    install=float(p.get(f"install_{system}",0) or 0) if system else 0.0
    CRF=float(p["CRF"]); maint=float(p["maint_rate"])/100.0
    EAC=n_units*(capex+install)*(CRF+maint)+E_kwh*float(p["elec_price"])+water_L*float(p["water_price"])
    co2=E_kwh*float(p["co2_factor"])
    return {"site":site,"system":system or "Baseline","n_units":int(n_units),"coverage":f_cov,"people_served_fraction":f_people,
            "comfort_share":comfort_share,"danger_share":danger_share,"stress_share":stress_ph/max(total_person_hours,1e-12),
            "EDH38_person_C_h_per_season_year":edh,"E_kWh_per_season_year":E_kwh,"water_L_per_season_year":water_L,
            "EAC_USD_per_season_year":EAC,"CO2_kg_per_season_year":co2,"mean_AT_served_C":float(np.average(AT_served,weights=np.maximum(Nserved,1e-9))),
            "mean_AT_unserved_C":float(np.average(AT_unserved,weights=np.maximum(Nunserved,1e-9))) if Nunserved.sum()>0 else float(np.mean(AT_unserved)),
            "mean_Tbase_C":float(np.mean(Tbase)),"mean_rho_moist_kg_m3":float(np.mean(rho_moist)),"mean_rho_dry_kg_m3":float(np.mean(rho_da)),
            "mean_w_out_kgkg":float(np.mean(wout)),"mean_phi":float(np.mean(phi)) if np.ndim(phi) else float(phi),
            "occupancy_clip_fraction":float(np.mean(N_uncapped>float(p["n_max"]))),"occupancy_mean_after_cap":float(np.mean(N)),
            "AHP_status":"PROVISIONAL_PAIRWISE_JUDGMENTS_NOT_APPROVED"}


def pareto_filter(df):
    """Keep nondominated alternatives: comfort higher, danger/cost/energy lower."""
    x=df[["comfort_share","danger_share","EAC_USD_per_season_year","E_kWh_per_season_year"]].to_numpy(float)
    keep=np.ones(len(df),bool)
    for i in range(len(df)):
        for j in range(len(df)):
            if i==j: continue
            # j dominates i if no worse in all objectives and strictly better in >=1.
            no_worse=(x[j,0]>=x[i,0]-1e-12 and x[j,1]<=x[i,1]+1e-12 and x[j,2]<=x[i,2]+1e-12 and x[j,3]<=x[i,3]+1e-12)
            strictly=(x[j,0]>x[i,0]+1e-12 or x[j,1]<x[i,1]-1e-12 or x[j,2]<x[i,2]-1e-12 or x[j,3]<x[i,3]-1e-12)
            if no_worse and strictly: keep[i]=False; break
    out=df.loc[keep].copy(); out["pareto"] = True
    return out


def ahp_weights_provisional():
    # Workbook's proposed judgments: Safety/Comfort=2, Safety/Cost=2, Safety/Energy=3,
    # Comfort/Cost=1, Comfort/Energy=2, Cost/Energy=2. All are DE_XUAT_CHO_NHOM.
    A=np.array([[1,2,2,3],[1/2,1,1,2],[1/2,1,1,2],[1/3,1/2,1/2,1]],float)
    eigvals,eigvecs=np.linalg.eig(A); idx=np.argmax(eigvals.real)
    w=np.abs(eigvecs[:,idx].real); w=w/w.sum(); lam=float(eigvals[idx].real)
    CI=(lam-4)/3; CR=CI/0.90
    return dict(zip(["Safety","Comfort","Cost","Energy"],w)),float(CR)


def add_provisional_ahp(pareto):
    d=pareto.copy()
    def max_benefit(x):
        lo,hi=np.nanmin(x),np.nanmax(x)
        return np.ones(len(x)) if hi-lo<1e-12 else (x-lo)/(hi-lo)
    def min_benefit(x):
        lo,hi=np.nanmin(x),np.nanmax(x)
        return np.ones(len(x)) if hi-lo<1e-12 else (hi-x)/(hi-lo)
    d["b_safety"]=min_benefit(d.danger_share.to_numpy())
    d["b_comfort"]=max_benefit(d.comfort_share.to_numpy())
    d["b_cost"]=min_benefit(d.EAC_USD_per_season_year.to_numpy())
    d["b_energy"]=min_benefit(d.E_kWh_per_season_year.to_numpy())
    w,cr=ahp_weights_provisional()
    d["AHP_score_provisional"]=(w["Safety"]*d.b_safety+w["Comfort"]*d.b_comfort+w["Cost"]*d.b_cost+w["Energy"]*d.b_energy)
    d["AHP_CR_provisional"]=cr
    d["AHP_status"]="PROVISIONAL_ONLY_GROUP_APPROVAL_REQUIRED"
    return d.sort_values("AHP_score_provisional",ascending=False)


def run_site(df,p,site,max_units=40):
    rows=[simulate_configuration(df,p,site,None,0)]
    for system in ["A","B","C"]:
        for n in range(1,max_units+1):
            rows.append(simulate_configuration(df,p,site,system,n))
    all_df=pd.DataFrame(rows)
    frontier=pareto_filter(all_df)
    ranked=add_provisional_ahp(frontier)
    return all_df,frontier,ranked


def run_sensitivity_oat(df,p,site):
    """Deterministic OAT +/-20%; no Monte Carlo. Baseline units from workbook data sheet."""
    cases=[]; base_units=int(p.get("n_units") or 1)
    parameters=["C_d","alpha_shelter","mix_eff","f_dung","U_mai","eps_evap","elec_price","water_price","occupancy_scale"]
    for system in ["A","B","C"]:
        base=simulate_configuration(df,p,site,system,base_units)
        for key in parameters:
            # f_dung is the workbook's vertical irradiance ratio; internally named f_dung in report.
            actual_key="f_dung" if key=="f_dung" else key
            workbook_key="f_dung" if key=="f_dung" else key
            if workbook_key not in p or p[workbook_key] is None: continue
            for factor in [0.8,1.2]:
                pp=dict(p); pp[workbook_key]=float(p[workbook_key])*factor
                result=simulate_configuration(df,pp,site,system,base_units)
                cases.append({"site":site,"system":system,"n_units":base_units,"parameter":key,"factor":factor,
                    "delta_comfort_share":result["comfort_share"]-base["comfort_share"],
                    "delta_danger_share":result["danger_share"]-base["danger_share"],
                    "delta_EAC_USD":result["EAC_USD_per_season_year"]-base["EAC_USD_per_season_year"],
                    "delta_E_kWh":result["E_kWh_per_season_year"]-base["E_kWh_per_season_year"],
                    "delta_water_L":result["water_L_per_season_year"]-base["water_L_per_season_year"]})
    return pd.DataFrame(cases)


def data_audit(data_dir):
    rows=[]
    for site,meta in VENUES.items():
        path=Path(data_dir)/meta["csv"]
        df=pd.read_csv(path,parse_dates=["Timestamp"])
        op=df.loc[df.Timestamp.dt.month.isin([6,7,8]) & df.Timestamp.dt.hour.between(10,21)]
        meta=VENUES[site]
        zen=solar_zenith_deg(df.Timestamp,meta["lat"],meta["lon"])
        cz=np.cos(np.deg2rad(zen)); valid=(zen<80.0)&(cz>0)
        dni=np.where(valid,df.Direct_Horizontal_Wm2.to_numpy(float)/np.maximum(cz,np.cos(np.deg2rad(80))),np.nan)
        rows.append({"site":site,"file":path.name,"rows":len(df),"columns":len(df.columns),"timestamp_start":str(df.Timestamp.min()),"timestamp_end":str(df.Timestamp.max()),
                     "duplicate_timestamps":int(df.Timestamp.duplicated().sum()),"missing_cells":int(df.isna().sum().sum()),"imputed_rows":int((df.Imputed.astype(float)!=0).sum()),
                     "operating_rows":len(op),"T_min_C":float(df.Temp_C.min()),"T_max_C":float(df.Temp_C.max()),"RH_min_pct":float(df.Relative_Humidity_Pct.min()),"RH_max_pct":float(df.Relative_Humidity_Pct.max()),
                     "GHI_max_Wm2":float(df.Solar_Radiation_Wm2.max()),"direct_horizontal_max_Wm2":float(df.Direct_Horizontal_Wm2.max()),"DNI_valid_fraction":float(np.mean(valid)),"DNI_est_p95_valid_Wm2":float(np.nanpercentile(dni,95))})
    return pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",default=".")
    ap.add_argument("--workbook",default="THONG_SO_THONG_NHAT_CHUAN_100_v2.3.2 (1).xlsx")
    ap.add_argument("--max-units",type=int,default=40)
    args=ap.parse_args(); root=Path(args.data_dir); wbpath=root/args.workbook
    if not wbpath.exists(): raise FileNotFoundError(f"Workbook not found: {wbpath}")
    tables=load_workbook_parameters(wbpath)
    audit=data_audit(root); audit.to_csv(root/"DataAudit_LeMinh_13Nov_v2.3.3.csv",index=False,encoding="utf-8-sig")
    all_results=[]; frontiers=[]; rankings=[]; sensitivities=[]
    for site,meta in VENUES.items():
        df=load_weather(root/meta["csv"],site)
        p=get_venue_params(tables,meta["workbook_col"])
        # workbook value for rho_water lives in cooling sheet
        p["rho_water"]=val(tables,"cooling","rho_water",meta["workbook_col"],1.0)
        all_df,front,rank=run_site(df,p,site,args.max_units)
        sens=run_sensitivity_oat(df,p,site)
        all_results.append(all_df); frontiers.append(front); rankings.append(rank); sensitivities.append(sens)
    pd.concat(all_results).to_csv(root/"ConfigurationResults_LeMinh_13Nov_v2.3.3.csv",index=False,encoding="utf-8-sig")
    pd.concat(frontiers).to_csv(root/"ParetoFrontier_LeMinh_13Nov_v2.3.3.csv",index=False,encoding="utf-8-sig")
    pd.concat(rankings).to_csv(root/"AHP_Ranking_Provisional_LeMinh_13Nov_v2.3.3.csv",index=False,encoding="utf-8-sig")
    pd.concat(sensitivities).to_csv(root/"Sensitivity_OAT_LeMinh_13Nov_v2.3.3.csv",index=False,encoding="utf-8-sig")
    print(audit.to_string(index=False))
    print("\nWrote DataAudit, ConfigurationResults, ParetoFrontier, provisional AHP ranking, and OAT sensitivity.")
    print("AHP judgments remain provisional; do not report the provisional ranking as group-approved.")

if __name__=="__main__": main()
