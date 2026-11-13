#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3D (Bulk) 材料彈性常數分析與繪圖腳本 (Cubic Projection 3 參數計算 VRH 版)

功能特點：
1. 自動檢測並略過缺失 (.sout 不存在) 或未完成 (無 'JOB DONE.') 的計算檔案，避免未收斂應力污染擬合。
2. 若預設擬合區間內有效點數不足，自動擴展應變區間進行降級擬合。
3. 輸出 9 組獨立彈性常數 (C11, C22, C33, C44, C55, C66, C12/C21, C13/C31, C23/C32) 與對稱化數值。
4. 以立方對稱投影公式 (Cubic Projection) 計算等效立方彈性常數：
   - C11_cubic = (C11 + C22 + C33) / 3
   - C12_cubic = (C12_sym + C13_sym + C23_sym) / 3
   - C44_cubic = (C44 + C55 + C66) / 3
   - C' (Tetragonal Shear Modulus) = (C11_cubic - C12_cubic) / 2
   - Zener 各向異性因子 A = 2 * C44_cubic / (C11_cubic - C12_cubic)
5. 僅由上述 3 個立方投影平均參數 (C11_cubic, C12_cubic, C44_cubic) 計算 Voigt-Reuss-Hill (VRH) 性質：
   - 體模量 (Bulk Modulus): B_VRH (= B_V = B_R)
   - 剪切模量 (Shear Modulus): G_VRH (含 G_V, G_R)
   - 楊氏模量 (Young's Modulus): E_VRH = 9*B*G / (3*B + G)
   - 蒲松比 (Poisson's Ratio): nu_VRH = (3*B - 2*G) / (2*(3*B + G))
   - Pugh 延脆性比 (B/G)、柯西壓力 (Cauchy Pressure)、普適各向異性指數 (A^U) 與 Born 穩定性判據。
"""

import os
import glob
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")  # 支援無 GUI 的遠端伺服器環境直接出圖
import matplotlib.pyplot as plt

# ================= 設定區域 =================
WORK_DIR = "/home/engfox/elastic_4revise_3D/Al05NbTiVZr05/relax_inputs"
SIGN_CORRECTION = -1.0           # QE 輸出應力正值代表壓力(壓縮)，乘以 -1 轉為標準拉伸正應力

MAX_TENSION = 0.05
NUM_TENSION = 5
MAX_COMPRESSION = 0.05
NUM_COMPRESSION = 5

# --- [擬合與容錯判定設定] ---
FIT_RANGE = [-0.03, 0.03]        # 預設線性擬合區間
MIN_FIT_POINTS = 3               # 擬合區間內最少需要的有效點數 (若不足則自動放寬區間，最低 2 點即可擬合)
REQUIRE_JOB_DONE = True          # 是否強制要求 .sout 內必須包含 "JOB DONE." 才納入計算
YIELD_THRESHOLD = 0.01           # 降伏點判定相對誤差閾值 (1%)


def generate_strain_mapping(max_t, num_t, max_c, num_c):
    """建立資料夾後綴 (如 m05, p01, orig) 與應變數值的對應字典"""
    mapping = {"orig": 0.0}
    if num_c > 0:
        c_vals = np.linspace(-max_c, 0, num_c + 1)[:-1]
        for i, val in enumerate(c_vals):
            mapping[f"m{len(c_vals)-i:02d}"] = round(float(val), 6)
    if num_t > 0:
        t_vals = np.linspace(0, max_t, num_t + 1)[1:]
        for i, val in enumerate(t_vals):
            mapping[f"p{i+1:02d}"] = round(float(val), 6)
    return mapping


def linear_fit(x, y):
    """
    穩健的一階線性回歸 (y = slope * x + intercept)
    不依賴 scipy，僅用 numpy 即可回傳 slope, intercept, r_squared
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-12 else 1.0
    return float(slope), float(intercept), float(r_squared)


def parse_sout(file_path, require_job_done=True):
    """
    解析 QE .sout 檔案：
    1. 檢查檔案是否存在且非空
    2. 檢查是否包含 'JOB DONE.' (避免讀取到中途崩潰或未收斂的應力)
    3. 從檔案末尾倒序擷取最後一組收斂的 3x3 應力張量 (kbar)
    回傳: (status_ok, reason_msg, (s_xx, s_yy, s_zz, s_yz, s_xz, s_xy))
    """
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return False, "檔案不存在或為空檔案", (None,) * 6

    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
    except Exception as e:
        return False, f"讀取錯誤: {e}", (None,) * 6

    if require_job_done:
        has_job_done = any("JOB DONE." in line for line in lines[-500:])
        if not has_job_done:
            has_job_done = any("JOB DONE." in line for line in lines)
        if not has_job_done:
            return False, "計算未完成 (缺少 JOB DONE.)", (None,) * 6

    stress_xx = stress_yy = stress_zz = None
    stress_xy = stress_yz = stress_xz = None

    for i in range(len(lines) - 1, -1, -1):
        line = lines[i]
        if "total" in line and "stress" in line and "kbar" in line:
            try:
                line_xx = lines[i + 1].split()
                line_yy = lines[i + 2].split()
                line_zz = lines[i + 3].split()
                stress_xx, stress_xy, stress_xz = float(line_xx[3]), float(line_xx[4]), float(line_xx[5])
                stress_yy, stress_yz = float(line_yy[4]), float(line_yy[5])
                stress_zz = float(line_zz[5])
                break
            except Exception:
                continue

    if stress_xx is None:
        return False, "找不到應力張量輸出 (total stress)", (None,) * 6

    return True, "OK", (stress_xx, stress_yy, stress_zz, stress_yz, stress_xz, stress_xy)


def get_info_from_name(folder_name, strain_mapping):
    """從資料夾名稱解析 (應變量, 應變方向類別, 材料名稱)"""
    if folder_name.endswith("-orig"):
        return 0.0, "orig", folder_name[:-5]

    for tag, stype in [
        ("-scaleX", "X"),
        ("-scaleY", "Y"),
        ("-scaleZ", "Z"),
        ("-shearYZ", "ShearYZ"),
        ("-shearXZ", "ShearXZ"),
        ("-shearXY", "ShearXY"),
    ]:
        if tag in folder_name:
            parts = folder_name.split(tag)
            mat_id = parts[0]
            suffix = parts[1]
            return strain_mapping.get(suffix), stype, mat_id

    return None, None, None


def select_fit_data(data_array, fit_range=FIT_RANGE, min_points=MIN_FIT_POINTS):
    """
    篩選線性擬合區間內的數據點。
    若預設區間內因缺檔導致點數少於 min_points，自動逐步放寬應變區間。
    """
    mask = (data_array[:, 0] >= fit_range[0] - 1e-6) & (data_array[:, 0] <= fit_range[1] + 1e-6)
    filtered = data_array[mask]
    if len(filtered) >= min_points:
        return filtered, fit_range

    for bound in [0.03, 0.04, 0.05]:
        mask_wide = (data_array[:, 0] >= -bound - 1e-6) & (data_array[:, 0] <= bound + 1e-6)
        filtered_wide = data_array[mask_wide]
        if len(filtered_wide) >= min_points:
            print(f"  [提示] 預設區間 {fit_range} 點數不足 ({len(filtered)} 點)，自動放寬擬合區間至 [{-bound:.2f}, {bound:.2f}] ({len(filtered_wide)} 點)")
            return filtered_wide, [-bound, bound]

    return data_array, [float(np.min(data_array[:, 0])), float(np.max(data_array[:, 0]))]


def find_yield_point(strains, stresses, slope, intercept):
    """尋找正向拉伸超越線性擬合區間後，偏離理想線性超過閾值 (1%) 的降伏點"""
    valid_indices = np.where(strains > FIT_RANGE[1] + 1e-6)[0]
    for idx in valid_indices:
        e = strains[idx]
        s_actual = stresses[idx]
        s_ideal = slope * e + intercept
        if abs(s_ideal) > 0.5:
            error = abs((s_actual - s_ideal) / s_ideal)
            if error > YIELD_THRESHOLD:
                return e, s_actual
    return None, None


def compute_polycrystalline_properties(c_dict):
    """
    由擬合出的 C_ij 分量計算：
    1. 對稱化 C12_sym, C13_sym, C23_sym
    2. 立方投影等效常數 (Cubic Projection):
       - C11_cubic = (C11 + C22 + C33) / 3
       - C12_cubic = (C12_sym + C13_sym + C23_sym) / 3
       - C44_cubic = (C44 + C55 + C66) / 3
       - C_prime = (C11_cubic - C12_cubic) / 2
       - Zener A = 2 * C44_cubic / (C11_cubic - C12_cubic)
    3. 僅以 C11_cubic, C12_cubic, C44_cubic 3 個參數計算 Voigt-Reuss-Hill (VRH) 多晶宏觀性質
    """
    def avg_valid(*vals):
        valid = [v for v in vals if v is not None and not np.isnan(v)]
        return float(np.mean(valid)) if valid else np.nan

    C11 = c_dict.get("C11", np.nan)
    C22 = c_dict.get("C22", np.nan)
    C33 = c_dict.get("C33", np.nan)
    C44 = c_dict.get("C44", np.nan)
    C55 = c_dict.get("C55", np.nan)
    C66 = c_dict.get("C66", np.nan)

    # 跨軸偶合項對稱化 (若單側缺失自動退回單側有效值)
    C12_sym = avg_valid(c_dict.get("C12"), c_dict.get("C21"))
    C13_sym = avg_valid(c_dict.get("C13"), c_dict.get("C31"))
    C23_sym = avg_valid(c_dict.get("C23"), c_dict.get("C32"))

    # === [1] 立方對稱投影平均 (Cubic Projection 3 參數) ===
    C11_cubic = avg_valid(C11, C22, C33)
    C12_cubic = avg_valid(C12_sym, C13_sym, C23_sym)
    C44_cubic = avg_valid(C44, C55, C66)

    C_prime = (C11_cubic - C12_cubic) / 2.0 if not (np.isnan(C11_cubic) or np.isnan(C12_cubic)) else np.nan
    Zener_A = (2.0 * C44_cubic) / (C11_cubic - C12_cubic) if (not np.isnan(C_prime) and C_prime != 0) else np.nan
    B_cubic = (C11_cubic + 2.0 * C12_cubic) / 3.0 if not (np.isnan(C11_cubic) or np.isnan(C12_cubic)) else np.nan

    props = {
        "C12_sym": C12_sym,
        "C13_sym": C13_sym,
        "C23_sym": C23_sym,
        "C11_cubic": C11_cubic,
        "C12_cubic": C12_cubic,
        "C44_cubic": C44_cubic,
        "C_prime": C_prime,
        "Zener_A": Zener_A,
        "B_cubic": B_cubic,
    }

    if any(np.isnan(x) for x in [C11_cubic, C12_cubic, C44_cubic]):
        return props

    # 僅以立方投影平均 3 個參數建立等效立方 6x6 勁度矩陣 C_mat
    C_mat = np.array([
        [C11_cubic, C12_cubic, C12_cubic, 0.0,       0.0,       0.0      ],
        [C12_cubic, C11_cubic, C12_cubic, 0.0,       0.0,       0.0      ],
        [C12_cubic, C12_cubic, C11_cubic, 0.0,       0.0,       0.0      ],
        [0.0,       0.0,       0.0,       C44_cubic, 0.0,       0.0      ],
        [0.0,       0.0,       0.0,       0.0,       C44_cubic, 0.0      ],
        [0.0,       0.0,       0.0,       0.0,       0.0,       C44_cubic],
    ], dtype=float)

    # Born 機械穩定性判據：等價於 C11 - C12 > 0, C11 + 2*C12 > 0, C44 > 0
    eigvals = np.linalg.eigvalsh(C_mat)
    born_stable = bool(np.all(eigvals > 0))

    # === [2] 立方 3 參數 Voigt 近似 (上限) ===
    B_V = (C11_cubic + 2.0 * C12_cubic) / 3.0
    G_V = (C11_cubic - C12_cubic + 3.0 * C44_cubic) / 5.0

    # === [3] 立方 3 參數 Reuss 近似 (下限) ===
    B_R = B_V
    denom_gr = 4.0 * C44_cubic + 3.0 * (C11_cubic - C12_cubic)
    G_R = (5.0 * (C11_cubic - C12_cubic) * C44_cubic) / denom_gr if denom_gr != 0 else np.nan

    # === [4] Hill 平均 (VRH) 模量與衍生性質 ===
    B_VRH = 0.5 * (B_V + B_R)
    G_VRH = 0.5 * (G_V + G_R) if not np.isnan(G_R) else G_V

    E_VRH = (9.0 * B_VRH * G_VRH) / (3.0 * B_VRH + G_VRH) if (3.0 * B_VRH + G_VRH) != 0 else np.nan
    nu_VRH = (3.0 * B_VRH - 2.0 * G_VRH) / (2.0 * (3.0 * B_VRH + G_VRH)) if (3.0 * B_VRH + G_VRH) != 0 else np.nan
    pugh_ratio = B_VRH / G_VRH if G_VRH != 0 else np.nan
    cauchy_pressure = C12_cubic - C44_cubic
    anisotropy_AU = 5.0 * (G_V / G_R) + (B_V / B_R) - 6.0 if (G_R > 0 and B_R > 0) else np.nan

    props.update({
        "C_matrix": C_mat,
        "B_V": B_V,
        "B_R": B_R,
        "B_VRH": B_VRH,
        "G_V": G_V,
        "G_R": G_R,
        "G_VRH": G_VRH,
        "E_VRH": E_VRH,
        "nu_VRH": nu_VRH,
        "Pugh_B_G": pugh_ratio,
        "Cauchy_Pressure": cauchy_pressure,
        "Anisotropy_AU": anisotropy_AU,
        "Born_Stable": born_stable,
    })
    return props


def plot_stress_strain(mat_id, strains_dict, fit_results, poly_props, work_dir):
    """繪製 6 個應變模式的應力-應變圖、擬合線與 Cubic / VRH 統計文字"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    c11_c = poly_props.get("C11_cubic", np.nan)
    c12_c = poly_props.get("C12_cubic", np.nan)
    c44_c = poly_props.get("C44_cubic", np.nan)
    b_vrh = poly_props.get("B_VRH", np.nan)
    g_vrh = poly_props.get("G_VRH", np.nan)
    e_vrh = poly_props.get("E_VRH", np.nan)

    title_str = (
        f"3D Bulk Mechanics: {mat_id} (GPa)\n"
        f"Cubic: C11={c11_c:.1f}, C12={c12_c:.1f}, C44={c44_c:.1f} GPa  |  "
        f"VRH (from Cubic): B={b_vrh:.1f}, G={g_vrh:.1f}, E={e_vrh:.1f} GPa"
    )
    fig.suptitle(title_str, fontsize=14, fontweight='bold')

    plot_configs = [
        ('X', 'X-Axis Tension', ['C11', 'C12', 'C13'], [1, 2, 3], ['Stress_XX', 'Stress_YY', 'Stress_ZZ'], ['#1f77b4', '#ff7f0e', '#2ca02c']),
        ('Y', 'Y-Axis Tension', ['C22', 'C21', 'C23'], [2, 1, 3], ['Stress_YY', 'Stress_XX', 'Stress_ZZ'], ['#ff7f0e', '#1f77b4', '#2ca02c']),
        ('Z', 'Z-Axis Tension', ['C33', 'C31', 'C32'], [3, 1, 2], ['Stress_ZZ', 'Stress_XX', 'Stress_YY'], ['#2ca02c', '#1f77b4', '#ff7f0e']),
        ('ShearYZ', 'YZ Shear', ['C44'], [4], ['Stress_YZ'], ['purple']),
        ('ShearXZ', 'XZ Shear', ['C55'], [5], ['Stress_XZ'], ['brown']),
        ('ShearXY', 'XY Shear', ['C66'], [6], ['Stress_XY'], ['red']),
    ]

    for ax, (key, title, c_labels, data_cols, y_labels, colors) in zip(axes, plot_configs):
        if not strains_dict.get(key):
            ax.set_title(f"{title} (No Valid Data)")
            ax.grid(True, linestyle=':', alpha=0.7)
            continue

        data = np.array(strains_dict[key])
        data = data[data[:, 0].argsort()]
        strains = data[:, 0]

        for col, y_label, color, c_label in zip(data_cols, y_labels, colors, c_labels):
            ax.plot(strains, data[:, col], marker='o', label=y_label, markersize=4, color=color)
            max_idx = np.argmax(data[:, col])
            ax.plot(strains[max_idx], data[max_idx, col], marker='*', color=color, markersize=10)

            if key in fit_results and c_label in fit_results[key]:
                slope = fit_results[key][c_label]
                intercept = fit_results[key][f"{c_label}_int"]
                ideal_line = slope * strains + intercept
                ax.plot(strains, ideal_line, '--', color=color, alpha=0.4, linewidth=1.5)

                if (key == 'X' and c_label == 'C11') or (key == 'Y' and c_label == 'C22') or \
                   (key == 'Z' and c_label == 'C33') or ('Shear' in key):
                    y_strain, y_stress = find_yield_point(strains, data[:, col], slope, intercept)
                    if y_strain is not None:
                        ax.axvline(y_strain, color='red', linestyle='-.', alpha=0.7)
                        ax.annotate(
                            f"Yield ({YIELD_THRESHOLD*100:.0f}% dev)\nStrain: {y_strain:.2f}",
                            xy=(y_strain, y_stress), xytext=(10, -30), textcoords="offset points",
                            ha='left', fontsize=8, color='red',
                            arrowprops=dict(arrowstyle="->", color='red')
                        )

        if key in fit_results:
            res = fit_results[key]
            n_fit = res.get("n_fit", 0)
            n_tot = len(strains)
            text_lines = [f"Points: {n_fit} fit / {n_tot} valid"]
            for c_label in c_labels:
                if c_label in res:
                    slope = res[c_label]
                    r2 = res.get(f"{c_label}_r2", np.nan)
                    text_lines.append(f"{c_label} = {slope:.1f} GPa (R²={r2:.3f})")
            if text_lines:
                props = dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray')
                ax.text(0.05, 0.95, "\n".join(text_lines), transform=ax.transAxes, fontsize=9.5,
                        verticalalignment='top', bbox=props)

        used_range = fit_results.get(key, {}).get("used_range", FIT_RANGE)
        ax.axvspan(used_range[0], used_range[1], color='gray', alpha=0.15, label='Fit Range')
        ax.set_title(title)
        ax.set_xlabel("Strain")
        ax.set_ylabel("3D Stress (GPa)")
        ax.grid(True, linestyle=':', alpha=0.7)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.axvline(0, color='black', linewidth=0.8)
        ax.legend(loc='lower right', fontsize=8)

    plt.tight_layout()
    plot_path = os.path.join(work_dir, f"{mat_id}_Analysis_3D_GPa.png")
    plt.savefig(plot_path, dpi=600)
    plt.close()
    print(f" 📊 分析圖表已儲存: {plot_path}")


def main():
    strain_mapping = generate_strain_mapping(MAX_TENSION, NUM_TENSION, MAX_COMPRESSION, NUM_COMPRESSION)
    data_groups = {}
    skipped_records = {}

    subdirs = sorted(glob.glob(os.path.join(WORK_DIR, "*")))
    for subdir in subdirs:
        if not os.path.isdir(subdir):
            continue
        folder_name = os.path.basename(subdir)
        strain_val, strain_type, mat_id = get_info_from_name(folder_name, strain_mapping)

        if mat_id is None or strain_val is None:
            continue

        if mat_id not in skipped_records:
            skipped_records[mat_id] = []

        sout_files = glob.glob(os.path.join(subdir, "*.sout"))
        if not sout_files:
            skipped_records[mat_id].append((folder_name, "找不到 .sout 檔案"))
            continue

        ok, msg, stresses = parse_sout(sout_files[0], require_job_done=REQUIRE_JOB_DONE)
        if not ok:
            skipped_records[mat_id].append((folder_name, msg))
            continue

        s_xx, s_yy, s_zz, s_yz, s_xz, s_xy = stresses

        # 轉換 kbar 到 GPa (1 kbar = 0.1 GPa) 並修正正負號
        sx_val = s_xx * 0.1 * SIGN_CORRECTION
        sy_val = s_yy * 0.1 * SIGN_CORRECTION
        sz_val = s_zz * 0.1 * SIGN_CORRECTION
        syz_val = s_yz * 0.1 * SIGN_CORRECTION
        sxz_val = s_xz * 0.1 * SIGN_CORRECTION
        sxy_val = s_xy * 0.1 * SIGN_CORRECTION

        if mat_id not in data_groups:
            data_groups[mat_id] = {'X': [], 'Y': [], 'Z': [], 'ShearYZ': [], 'ShearXZ': [], 'ShearXY': []}

        point_tuple = (strain_val, sx_val, sy_val, sz_val, syz_val, sxz_val, sxy_val)
        if strain_type == "orig":
            for k in data_groups[mat_id].keys():
                data_groups[mat_id][k].append(point_tuple)
        else:
            data_groups[mat_id][strain_type].append(point_tuple)

    # 輸出略過清單摘要
    for mat_id, skips in skipped_records.items():
        if mat_id not in data_groups:
            print(f"\n⚠️ 材料 [{mat_id}] 無任何完成的計算檔案 (共略過 {len(skips)} 個資料夾)，跳過分析。")
        elif skips:
            print(f"\nℹ️ 材料 [{mat_id}] 自動略過 {len(skips)} 個未完成/缺失的案例：")
            for fname, reason in skips:
                print(f"   - {fname}: {reason}")

    # CSV 欄位清單
    fieldnames = [
        "Material",
        # --- [1] 立方等效彈性常數 (Cubic Projection) ---
        "C11_cubic_GPa",
        "C12_cubic_GPa",
        "C44_cubic_GPa",
        "C_prime_cubic_GPa",
        "Zener_A_cubic",
        # --- [2] Voigt-Reuss-Hill (VRH) 模量與力學性質 (由 Cubic 3 參數計算) ---
        "B_VRH_GPa",
        "G_VRH_GPa",
        "E_VRH_GPa",
        "Poisson_nu_VRH",
        "Pugh_B_G",
        "B_V_GPa",
        "B_R_GPa",
        "G_V_GPa",
        "G_R_GPa",
        "Cauchy_Pressure_GPa",
        "Anisotropy_AU",
        "Born_Stable",
        # --- [3] 原始 9 組獨立分量與對稱化分量 ---
        "C11_GPa", "C22_GPa", "C33_GPa",
        "C44_GPa", "C55_GPa", "C66_GPa",
        "C12_sym_GPa", "C13_sym_GPa", "C23_sym_GPa",
        "C12_GPa", "C21_GPa", "C13_GPa", "C31_GPa", "C23_GPa", "C32_GPa"
    ]

    csv_data_list = []

    for mat_id, strains_dict in data_groups.items():
        print(f"\n=================================================================")
        print(f" 材料: {mat_id} (Cubic 3 參數彈性常數與 VRH 分析中... 單位: GPa)")
        print(f"=================================================================")
        row_data = {f: "N/A" for f in fieldnames}
        row_data["Material"] = mat_id

        fit_results = {'X': {}, 'Y': {}, 'Z': {}, 'ShearYZ': {}, 'ShearXZ': {}, 'ShearXY': {}}
        raw_c = {}

        # 1. X 方向拉伸/壓縮 -> C11, C12, C13
        if len(strains_dict['X']) >= 2:
            dx, u_range = select_fit_data(np.array(strains_dict['X']))
            if len(dx) >= 2:
                s11, i11, r11 = linear_fit(dx[:, 0], dx[:, 1])
                s12, i12, r12 = linear_fit(dx[:, 0], dx[:, 2])
                s13, i13, r13 = linear_fit(dx[:, 0], dx[:, 3])
                fit_results['X'] = {
                    'C11': s11, 'C11_int': i11, 'C11_r2': r11,
                    'C12': s12, 'C12_int': i12, 'C12_r2': r12,
                    'C13': s13, 'C13_int': i13, 'C13_r2': r13,
                    'n_fit': len(dx), 'used_range': u_range
                }
                raw_c.update({'C11': s11, 'C12': s12, 'C13': s13})
                row_data["C11_GPa"], row_data["C12_GPa"], row_data["C13_GPa"] = f"{s11:.2f}", f"{s12:.2f}", f"{s13:.2f}"

        # 2. Y 方向拉伸/壓縮 -> C22, C21, C23
        if len(strains_dict['Y']) >= 2:
            dy, u_range = select_fit_data(np.array(strains_dict['Y']))
            if len(dy) >= 2:
                s22, i22, r22 = linear_fit(dy[:, 0], dy[:, 2])
                s21, i21, r21 = linear_fit(dy[:, 0], dy[:, 1])
                s23, i23, r23 = linear_fit(dy[:, 0], dy[:, 3])
                fit_results['Y'] = {
                    'C22': s22, 'C22_int': i22, 'C22_r2': r22,
                    'C21': s21, 'C21_int': i21, 'C21_r2': r21,
                    'C23': s23, 'C23_int': i23, 'C23_r2': r23,
                    'n_fit': len(dy), 'used_range': u_range
                }
                raw_c.update({'C22': s22, 'C21': s21, 'C23': s23})
                row_data["C22_GPa"], row_data["C21_GPa"], row_data["C23_GPa"] = f"{s22:.2f}", f"{s21:.2f}", f"{s23:.2f}"

        # 3. Z 方向拉伸/壓縮 -> C33, C31, C32
        if len(strains_dict['Z']) >= 2:
            dz, u_range = select_fit_data(np.array(strains_dict['Z']))
            if len(dz) >= 2:
                s33, i33, r33 = linear_fit(dz[:, 0], dz[:, 3])
                s31, i31, r31 = linear_fit(dz[:, 0], dz[:, 1])
                s32, i32, r32 = linear_fit(dz[:, 0], dz[:, 2])
                fit_results['Z'] = {
                    'C33': s33, 'C33_int': i33, 'C33_r2': r33,
                    'C31': s31, 'C31_int': i31, 'C31_r2': r31,
                    'C32': s32, 'C32_int': i32, 'C32_r2': r32,
                    'n_fit': len(dz), 'used_range': u_range
                }
                raw_c.update({'C33': s33, 'C31': s31, 'C32': s32})
                row_data["C33_GPa"], row_data["C31_GPa"], row_data["C32_GPa"] = f"{s33:.2f}", f"{s31:.2f}", f"{s32:.2f}"

        # 4. YZ 剪切 -> C44
        if len(strains_dict['ShearYZ']) >= 2:
            ds, u_range = select_fit_data(np.array(strains_dict['ShearYZ']))
            if len(ds) >= 2:
                s44, i44, r44 = linear_fit(ds[:, 0], ds[:, 4])
                fit_results['ShearYZ'] = {'C44': s44, 'C44_int': i44, 'C44_r2': r44, 'n_fit': len(ds), 'used_range': u_range}
                raw_c['C44'] = s44
                row_data["C44_GPa"] = f"{s44:.2f}"

        # 5. XZ 剪切 -> C55
        if len(strains_dict['ShearXZ']) >= 2:
            ds, u_range = select_fit_data(np.array(strains_dict['ShearXZ']))
            if len(ds) >= 2:
                s55, i55, r55 = linear_fit(ds[:, 0], ds[:, 5])
                fit_results['ShearXZ'] = {'C55': s55, 'C55_int': i55, 'C55_r2': r55, 'n_fit': len(ds), 'used_range': u_range}
                raw_c['C55'] = s55
                row_data["C55_GPa"] = f"{s55:.2f}"

        # 6. XY 剪切 -> C66
        if len(strains_dict['ShearXY']) >= 2:
            ds, u_range = select_fit_data(np.array(strains_dict['ShearXY']))
            if len(ds) >= 2:
                s66, i66, r66 = linear_fit(ds[:, 0], ds[:, 6])
                fit_results['ShearXY'] = {'C66': s66, 'C66_int': i66, 'C66_r2': r66, 'n_fit': len(ds), 'used_range': u_range}
                raw_c['C66'] = s66
                row_data["C66_GPa"] = f"{s66:.2f}"

        # 計算對稱化矩陣、Cubic 投影等效常數與多晶宏觀彈性模量 (VRH)
        poly = compute_polycrystalline_properties(raw_c)

        # 填入 CSV 數值
        for k_csv, k_poly in [
            ("C11_cubic_GPa", "C11_cubic"),
            ("C12_cubic_GPa", "C12_cubic"),
            ("C44_cubic_GPa", "C44_cubic"),
            ("C_prime_cubic_GPa", "C_prime"),
            ("C12_sym_GPa", "C12_sym"),
            ("C13_sym_GPa", "C13_sym"),
            ("C23_sym_GPa", "C23_sym"),
            ("B_V_GPa", "B_V"),
            ("B_R_GPa", "B_R"),
            ("B_VRH_GPa", "B_VRH"),
            ("G_V_GPa", "G_V"),
            ("G_R_GPa", "G_R"),
            ("G_VRH_GPa", "G_VRH"),
            ("E_VRH_GPa", "E_VRH"),
            ("Cauchy_Pressure_GPa", "Cauchy_Pressure"),
        ]:
            val = poly.get(k_poly, np.nan)
            if not np.isnan(val):
                row_data[k_csv] = f"{val:.2f}"

        for k_csv, k_poly in [
            ("Zener_A_cubic", "Zener_A"),
            ("Poisson_nu_VRH", "nu_VRH"),
            ("Pugh_B_G", "Pugh_B_G"),
            ("Anisotropy_AU", "Anisotropy_AU"),
        ]:
            val = poly.get(k_poly, np.nan)
            if not np.isnan(val):
                row_data[k_csv] = f"{val:.4f}"

        if "Born_Stable" in poly:
            row_data["Born_Stable"] = "Yes" if poly["Born_Stable"] else "No"

        # 輸出終端機完整結構化報告
        print(f"\n  [1. 等效立方彈性常數 (Cubic Projection)]")
        print(f"    C11_cubic = {row_data['C11_cubic_GPa']:>7} GPa  (三向正向拉壓勁度平均)")
        print(f"    C12_cubic = {row_data['C12_cubic_GPa']:>7} GPa  (跨軸泊松耦合勁度平均)")
        print(f"    C44_cubic = {row_data['C44_cubic_GPa']:>7} GPa  (三向純剪切勁度平均)")
        print(f"    C' (Shear)= {row_data['C_prime_cubic_GPa']:>7} GPa  ((C11 - C12) / 2)")
        print(f"    Zener A   = {row_data['Zener_A_cubic']:>7}      (2*C44 / (C11 - C12))")

        print(f"\n  [2. Voigt-Reuss-Hill (VRH) 多晶宏觀力學性質 (由 Cubic 3 參數計算)]")
        if "B_VRH" in poly:
            print(f"    體模量 (Bulk Modulus, B)   : B_VRH = {poly['B_VRH']:>7.2f} GPa  (B_V = {poly['B_V']:.2f}, B_R = {poly['B_R']:.2f})")
            print(f"    剪切模量 (Shear Mod., G)   : G_VRH = {poly['G_VRH']:>7.2f} GPa  (G_V = {poly['G_V']:.2f}, G_R = {poly['G_R']:.2f})")
            print(f"    楊氏模量 (Young's Mod., E) : E_VRH = {poly['E_VRH']:>7.2f} GPa")
            print(f"    蒲松比 (Poisson Ratio, ν)  : ν_VRH = {poly['nu_VRH']:>7.4f}")
            print(f"    Pugh 比 (B/G)              : {poly['Pugh_B_G']:>7.4f}  ({'延展性 Ductile (>1.75)' if poly['Pugh_B_G'] > 1.75 else '脆性 Brittle (<1.75)'})")
            print(f"    柯西壓力 (Cauchy Pressure) : {poly['Cauchy_Pressure']:>7.2f} GPa  (C12_cubic - C44_cubic)")
            print(f"    普適各向異性 (Anisotropy)  : A^U   = {poly['Anisotropy_AU']:>7.4f}")
            print(f"    Born 機械穩定性判據        : {'穩定 (Stable - 立方條件滿足)' if poly['Born_Stable'] else '不穩定 (Unstable)'}")

        print(f"\n  [3. 原始獨立彈性常數 (9 Components & Symmetry)]")
        print(f"    正向分量: C11 = {row_data['C11_GPa']:>7} | C22 = {row_data['C22_GPa']:>7} | C33 = {row_data['C33_GPa']:>7}")
        print(f"    剪切分量: C44 = {row_data['C44_GPa']:>7} | C55 = {row_data['C55_GPa']:>7} | C66 = {row_data['C66_GPa']:>7}")
        print(f"    耦合分量: C12 = {row_data['C12_GPa']:>7} (C21 = {row_data['C21_GPa']:>7}) -> 對稱化 C12_sym = {row_data['C12_sym_GPa']:>7}")
        print(f"              C13 = {row_data['C13_GPa']:>7} (C31 = {row_data['C31_GPa']:>7}) -> 對稱化 C13_sym = {row_data['C13_sym_GPa']:>7}")
        print(f"              C23 = {row_data['C23_GPa']:>7} (C32 = {row_data['C32_GPa']:>7}) -> 對稱化 C23_sym = {row_data['C23_sym_GPa']:>7}")

        plot_stress_strain(mat_id, strains_dict, fit_results, poly, WORK_DIR)
        csv_data_list.append(row_data)

    if csv_data_list:
        csv_filepath = os.path.join(WORK_DIR, "elastic_summary_3D_GPa.csv")
        with open(csv_filepath, mode='w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_data_list)
        print(f"\n✅ 分析完成！總表已匯出至: {csv_filepath}")


if __name__ == "__main__":
    main()