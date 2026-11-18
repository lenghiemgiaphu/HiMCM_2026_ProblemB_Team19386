"""
MonteCarloV24_LeMinh_15Nov_v1.py -- HiMCM 2026 Problem B (Team #19386)
Monte Carlo 10.000 kịch bản CHẠY TRÊN MÔ HÌNH V2.4.0 (simulate_configuration của nhóm), không dùng mô hình giữ chỗ.

Mỗi kịch bản (draw) ngẫu nhiên hóa:
  (a) Thời tiết: year-block bootstrap -- rút có hoàn lại 7 mùa hè (2020..2026) từ dữ liệu thật, giữ nguyên cấu trúc giờ/ngày
  (b) Tham số kém chắc chắn: nhân hệ số ~ Uniform(0.8, 1.2) (khoảng +-20% theo quy ước độ nhạy của nhóm) cho
      C_d, alpha_shelter, mix_eff, f_dung, U_mai, eps_evap, elec_price, water_price, occupancy_scale, capex+install
      (alpha_shelter, mix_eff, f_dung, eps_evap bị chặn trên bằng 1)
Cấu hình đánh giá cho mỗi địa điểm: Baseline + cấu hình tốt nhất của từng hệ A/B/C theo AHP chính thức (giới hạn quét 40).
Đầu ra: phân phối (trung bình, độ lệch, khoảng 2,5-97,5%), CI95 của trung bình, hiệu suất người·°C·giờ/kWh,
        xác suất đứng hạng 1 theo AHP chính thức (chuẩn hóa gộp trên toàn bộ mẫu).
seed = 42 (rng riêng theo [42, địa điểm, khối] nên kết quả tái lập được dù chạy song song).
"""
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

np.random.seed(42)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ModelExt_LeMinh_15Nov_v1 as X  # noqa: E402

OUT = Path(os.environ.get("HIMCM_OUT", "."))
OUT.mkdir(parents=True, exist_ok=True)
TAG = "LeMinh_15Nov_v1"
N_TOTAL = 10_000
CHUNK = 500
T0 = time.time()
PARAMS = ["C_d", "alpha_shelter", "mix_eff", "f_dung", "U_mai", "eps_evap", "elec_price", "water_price",
          "occupancy_scale", "capex_install_total"]


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


PM, TABLES, SITES = X.load_ctx()
SITE_LIST = list(SITES)

# --- cấu hình đánh giá: Baseline + tốt nhất mỗi hệ (AHP chính thức, quét 40)
CONFIGS = {}
for s in SITE_LIST:
    df, p = SITES[s]
    _, _, r = X.sweep_custom(PM, s, df, p, 40)
    cfg = [("Baseline", 0)]
    for sysm in X.SYSTEMS:
        g = r[r.system == sysm]
        if len(g):
            cfg.append((sysm, int(g.sort_values("AHP_score", ascending=False).iloc[0].n_units)))
    CONFIGS[s] = cfg

# --- chuẩn bị khối năm cho bootstrap
YEARS = list(range(2020, 2027))
YB = {}
for s in SITE_LIST:
    df, _ = SITES[s]
    yr = df.Timestamp.dt.year.to_numpy()
    cols = {c: {y: df.loc[yr == y, c].to_numpy(float) for y in YEARS}
            for c in ["Temp_C", "Relative_Humidity_Pct", "Solar_Radiation_Wm2", "Wind_Speed_ms", "Cloud_Cover_Pct", "WetBulb_C"]}
    ts = {}
    for y in YEARS:
        base_ts = pd.DatetimeIndex(df.loc[yr == y, "Timestamp"])
        for k in range(7):
            ts[(y, k)] = base_ts + pd.DateOffset(years=(2000 + k) - y)
    YB[s] = (cols, ts)


def make_df(site, rng, identity=False):
    cols, ts = YB[site]
    ys = YEARS if identity else list(rng.choice(YEARS, size=7, replace=True))
    d = {c: np.concatenate([cols[c][y] for y in ys]) for c in cols}
    d["Timestamp"] = pd.Series(ts[(ys[0], 0)].append([ts[(y, k)] for k, y in enumerate(ys) if k > 0]))
    return pd.DataFrame(d)


def draw_params(site, rng, identity=False):
    p = SITES[site][1]
    for name in PARAMS:
        v = 1.0 if identity else rng.uniform(0.8, 1.2)
        p = X.apply_variant(p, name, v)
    return p


METRICS = ["danger_share", "comfort_share", "EDH38_person_C_h_per_season_year", "E_kWh_per_season_year",
           "water_L_per_season_year", "EAC_USD_per_season_year", "mean_AT_population_C"]


def one_draw(site, rng, identity=False):
    df = make_df(site, rng, identity)
    p = draw_params(site, rng, identity)
    out = []
    for system, n in CONFIGS[site]:
        r = PM.simulate_configuration(df, p, site, None if system == "Baseline" else system, n)
        out.append([r[m] for m in METRICS])
    return np.array(out)    # (n_cfg, n_metrics)


def chunk_task(args):
    si, ci, n = args
    site = SITE_LIST[si]
    rng = np.random.default_rng([42, si, ci])
    return si, ci, np.stack([one_draw(site, rng) for _ in range(n)])


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else N_TOTAL
    # --- kiểm tra tự động: identity draw phải trùng mô hình tất định
    for s in SITE_LIST:
        a = one_draw(s, None, identity=True)
        df, p = SITES[s]
        for (system, n), row in zip(CONFIGS[s], a):
            r = PM.simulate_configuration(df, p, s, None if system == "Baseline" else system, n)
            assert abs(r["danger_share"] - row[0]) < 1e-9 and abs(r["EAC_USD_per_season_year"] - row[5]) < 1e-6, (s, system, n)
    log("kiểm tra identity-draw == mô hình tất định: OK")

    tasks = []
    for si in range(len(SITE_LIST)):
        n_chunks = int(np.ceil(N / CHUNK))
        for ci in range(n_chunks):
            tasks.append((si, ci, min(CHUNK, N - ci * CHUNK)))
    res = {si: {} for si in range(len(SITE_LIST))}
    with Pool(2) as pool:
        for k, (si, ci, arr) in enumerate(pool.imap_unordered(chunk_task, tasks), 1):
            res[si][ci] = arr
            if k % 8 == 0 or k == len(tasks):
                log(f"{k}/{len(tasks)} khối xong")

    W, _ = PM.ahp_weights_official()
    long_rows, summ = [], []
    for si, s in enumerate(SITE_LIST):
        A = np.concatenate([res[si][c] for c in sorted(res[si])])       # (N, n_cfg, n_metrics)
        cfgs = CONFIGS[s]
        base = A[:, 0, :]
        det = {}
        df, p = SITES[s]
        for (system, n) in cfgs:
            r = PM.simulate_configuration(df, p, s, None if system == "Baseline" else system, n)
            det[(system, n)] = r
        # chuẩn hóa gộp để tính P(hạng 1) theo AHP chính thức
        crit = np.stack([A[:, :, 0], A[:, :, 1], A[:, :, 5], A[:, :, 3]], axis=2)   # danger, comfort, EAC, kWh
        lo, hi = crit.min(axis=(0, 1)), crit.max(axis=(0, 1))
        span = np.where(hi - lo < 1e-12, 1.0, hi - lo)
        b = np.empty_like(crit)
        b[:, :, 0] = (hi[0] - crit[:, :, 0]) / span[0]      # Safety: danger thấp tốt
        b[:, :, 1] = (crit[:, :, 1] - lo[1]) / span[1]      # Comfort cao tốt
        b[:, :, 2] = (hi[2] - crit[:, :, 2]) / span[2]      # Cost thấp tốt
        b[:, :, 3] = (hi[3] - crit[:, :, 3]) / span[3]      # Energy thấp tốt
        score = b @ np.array([W["Safety"], W["Comfort"], W["Cost"], W["Energy"]])
        win = (score == score.max(axis=1, keepdims=True)).astype(float)
        win = win / win.sum(axis=1, keepdims=True)
        for j, (system, n) in enumerate(cfgs):
            x = A[:, j, :]
            dEDH = base[:, 2] - x[:, 2]
            eff = np.where(x[:, 3] > 0, dEDH / np.where(x[:, 3] > 0, x[:, 3], 1), np.nan)
            dpp = (base[:, 0] - x[:, 0]) * 100
            series = {"danger_share": x[:, 0], "comfort_share": x[:, 1], "EDH38": x[:, 2], "kWh": x[:, 3],
                      "water_L": x[:, 4], "EAC_USD": x[:, 5], "mean_AT_pop_C": x[:, 6],
                      "danger_pp_avoided": dpp, "dEDH38": dEDH, "Eff_per_kWh": eff, "AHP_score_pooled": score[:, j]}
            detm = det[(system, n)]
            detmap = {"danger_share": detm["danger_share"], "comfort_share": detm["comfort_share"],
                      "EDH38": detm["EDH38_person_C_h_per_season_year"], "kWh": detm["E_kWh_per_season_year"],
                      "water_L": detm["water_L_per_season_year"], "EAC_USD": detm["EAC_USD_per_season_year"],
                      "mean_AT_pop_C": detm["mean_AT_population_C"]}
            for m, v in series.items():
                v = v[~np.isnan(v)]
                if len(v) == 0:
                    continue
                summ.append(dict(site=s, system=system, n_units=n, metric=m, N=len(v), mean=v.mean(), sd=v.std(ddof=1),
                                 CI95_mean_half=1.96 * v.std(ddof=1) / np.sqrt(len(v)),
                                 p2_5=np.percentile(v, 2.5), p50=np.percentile(v, 50), p97_5=np.percentile(v, 97.5),
                                 deterministic=detmap.get(m, np.nan)))
            summ.append(dict(site=s, system=system, n_units=n, metric="P_rank1_AHP", N=len(A), mean=win[:, j].mean(),
                             sd=np.nan, CI95_mean_half=np.nan, p2_5=np.nan, p50=np.nan, p97_5=np.nan, deterministic=np.nan))
            long_rows.append(pd.DataFrame({"site": s, "system": system, "n_units": n, "draw": np.arange(len(A)),
                                           **{k: v for k, v in series.items()}}))
    pd.DataFrame(summ).to_csv(OUT / f"MonteCarloSummary_{TAG}.csv", index=False, encoding="utf-8-sig")
    pd.concat(long_rows, ignore_index=True).to_csv(OUT / f"MonteCarloSamples_{TAG}.csv.gz", index=False, compression="gzip")
    log(f"DONE  N={N} draws/site  -> MonteCarloSummary_{TAG}.csv , MonteCarloSamples_{TAG}.csv.gz")
