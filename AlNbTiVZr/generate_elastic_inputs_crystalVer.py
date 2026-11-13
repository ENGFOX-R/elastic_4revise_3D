import os
import glob
import numpy as np
import re

# ================= 設定區域 =================
INPUT_ROOT = "/home/engfox/elastic_4revise_3D/vc-relax_output/"
OUTPUT_ROOT = "/home/engfox/elastic_4revise_3D/relax_inputs/"

MAX_TENSION = 0.05       # 最大拉伸 5%
NUM_TENSION = 5          # 拉伸切 5 分

MAX_COMPRESSION = 0.05   # 最大壓縮 5%
NUM_COMPRESSION = 5      # 壓縮切 5 分

def generate_strain_list(max_t, num_t, max_c, num_c):
    strains = []
    if num_c > 0:
        c_vals = np.linspace(-max_c, 0, num_c + 1)[:-1]
        for i, val in enumerate(c_vals):
            strains.append((val, f"m{len(c_vals)-i:02d}"))
    strains.append((0.0, "orig"))
    if num_t > 0:
        t_vals = np.linspace(0, max_t, num_t + 1)[1:]
        for i, val in enumerate(t_vals):
            strains.append((val, f"p{i+1:02d}"))
    return strains

def parse_sout_final_structure(sout_path):
    with open(sout_path, 'r') as f:
        lines = f.readlines()
    cell_params, atom_names, atom_coords = [], [], []
    cell_idx, pos_idx = -1, -1
    pos_type = "angstrom"
    
    # 從檔案末尾倒找最後一組更新的結構
    for i in range(len(lines)-1, -1, -1):
        if "CELL_PARAMETERS" in lines[i] and cell_idx == -1: cell_idx = i
        if "ATOMIC_POSITIONS" in lines[i] and pos_idx == -1: 
            pos_idx = i
            if "crystal" in lines[i].lower(): pos_type = "crystal"
        if cell_idx != -1 and pos_idx != -1: break
    if cell_idx == -1 or pos_idx == -1: return None, None, None, None
    
    try:
        # 1. 讀取晶格向量 (3 行)
        for i in range(1, 4):
            parts = lines[cell_idx + i].split()
            cell_params.append([float(parts[0]), float(parts[1]), float(parts[2])])
        cell_params = np.array(cell_params)
        
        # 2. 讀取原子座標 (使用 try-except 判斷是否為浮點數，避免遇到 Writing config-only... 崩潰)
        current_idx = pos_idx + 1
        while current_idx < len(lines):
            line = lines[current_idx].strip()
            parts = line.split()
            if len(parts) >= 4:
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    atom_names.append(parts[0])
                    atom_coords.append([x, y, z])
                except ValueError:
                    # 遇到字串無法轉浮點數（例如 Writing config-only...），代表原子區塊結束
                    break
            else:
                break 
            current_idx += 1
        return cell_params, atom_names, np.array(atom_coords), pos_type
    except Exception as e:
        print(f"  [解析失敗] {sout_path} 錯誤訊息: {e}")
        return None, None, None, None

def get_fractional_coords(cell, coords_angstrom):
    return np.dot(coords_angstrom, np.linalg.inv(cell))

def apply_strain_tensor(cell, strain_val, strain_type):
    strain_tensor = np.zeros((3, 3))
    if strain_type == 'X': strain_tensor[0, 0] = strain_val
    elif strain_type == 'Y': strain_tensor[1, 1] = strain_val
    elif strain_type == 'Z': strain_tensor[2, 2] = strain_val
    elif strain_type == 'ShearYZ':
        strain_tensor[1, 2] = strain_val / 2.0
        strain_tensor[2, 1] = strain_val / 2.0
    elif strain_type == 'ShearXZ':
        strain_tensor[0, 2] = strain_val / 2.0
        strain_tensor[2, 0] = strain_val / 2.0
    elif strain_type == 'ShearXY':
        strain_tensor[0, 1] = strain_val / 2.0
        strain_tensor[1, 0] = strain_val / 2.0
    return np.dot(cell, np.eye(3) + strain_tensor)

def create_relax_input(original_in_path, output_path, new_cell, atom_names, new_atoms_angstrom, strain_info):
    with open(original_in_path, 'r') as f:
        content = f.read()
    
    # 修正為 relax 計算並開啟應力輸出
    content = re.sub(r'calculation\s*=\s*[\'"]vc-relax[\'"]', 'calculation = "relax"', content, flags=re.IGNORECASE)
    if "tstress" not in content.lower():
        content = re.sub(r'(calculation\s*=\s*"relax")', r'\1\n tstress = .true.', content, flags=re.IGNORECASE)
    
    lines = content.split('\n')
    new_lines = []
    skip_mode = None  
    
    for line in lines:
        u_line = line.upper().strip()
        
        # 遇到 &CELL 區塊，跳過直到 '/'
        if u_line.startswith("&CELL"):
            skip_mode = "NAMELIST_CELL"
            continue
        if skip_mode == "NAMELIST_CELL":
            if u_line == "/": skip_mode = None
            continue
            
        # 遇到晶格或原子卡片開頭，切換為跳過模式
        if u_line.startswith("CELL_PARAMETERS"):
            skip_mode = "CELL_PARAM"
            continue
        if u_line.startswith("ATOMIC_POSITIONS"):
            skip_mode = "ATOMIC_POS"
            continue
            
        # 在跳過模式下，只要遇到其他合法關鍵字卡片/區塊即解除跳過
        if skip_mode in ["CELL_PARAM", "ATOMIC_POS"]:
            if any(u_line.startswith(kw) for kw in ["&", "ATOMIC_SPECIES", "K_POINTS", "CONSTRAINTS", "OCCUPATIONS", "ATOMIC_FORCES"]):
                skip_mode = None
            else:
                continue
                
        if skip_mode is None:
            new_lines.append(line)

    with open(output_path, 'w') as f:
        for line in new_lines: f.write(line + "\n")
        f.write(f"\n! --- 3D Bulk Elastic Data: {strain_info} ---\n")
        f.write("CELL_PARAMETERS (angstrom)\n")
        for v in new_cell: f.write(f"{v[0]:14.9f} {v[1]:14.9f} {v[2]:14.9f}\n")
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for name, coord in zip(atom_names, new_atoms_angstrom):
            f.write(f"{name:3s} {coord[0]:14.9f} {coord[1]:14.9f} {coord[2]:14.9f}\n")

def main():
    print("========== 3D Material Elastic Constants Input Generator ==========")
    STRAIN_LIST = generate_strain_list(MAX_TENSION, NUM_TENSION, MAX_COMPRESSION, NUM_COMPRESSION)
    print(f"[系統] 已規劃 {len(STRAIN_LIST)} 組應變梯度 (拉伸 {NUM_TENSION}, 壓縮 {NUM_COMPRESSION})")

    if not os.path.exists(OUTPUT_ROOT): os.makedirs(OUTPUT_ROOT)
    
    # 修正搜尋邏輯：同時檢查直接放在目錄下的 .sout 檔，以及子目錄中的 .sout 檔
    sout_files = glob.glob(os.path.join(INPUT_ROOT, "*.sout")) + glob.glob(os.path.join(INPUT_ROOT, "*/*.sout"))
    if not sout_files:
        print(f"[錯誤] 在 {INPUT_ROOT} 下找不到任何 .sout 檔案！")
        return
        
    for sout_path in sout_files:
        folder_name = os.path.splitext(os.path.basename(sout_path))[0]
        dir_name = os.path.dirname(sout_path)
        in_path = os.path.join(dir_name, f"{folder_name}.in")
        
        if not os.path.exists(in_path):
            print(f"  [跳過] 找不到 {sout_path} 對應的 {in_path}")
            continue
        
        base_cell, atom_names, base_coords, pos_type = parse_sout_final_structure(sout_path)
        if base_cell is None:
            print(f"  [跳過] 無法解析結構: {sout_path}")
            continue
            
        frac_coords = base_coords if pos_type == "crystal" else get_fractional_coords(base_cell, base_coords)

        for strain, suffix in STRAIN_LIST:
            for def_type, folder_ext in [('X', '-scaleX'), ('Y', '-scaleY'), ('Z', '-scaleZ'), 
                                         ('ShearYZ', '-shearYZ'), ('ShearXZ', '-shearXZ'), ('ShearXY', '-shearXY')]:
                if strain == 0.0:
                    if def_type != 'X': continue
                    new_folder_name = f"{folder_name}-orig"
                else:
                    new_folder_name = f"{folder_name}{folder_ext}{suffix}"
                
                new_cell = apply_strain_tensor(base_cell, strain, def_type)
                new_atoms_ang = np.dot(frac_coords, new_cell)
                
                new_folder_path = os.path.join(OUTPUT_ROOT, new_folder_name)
                if not os.path.exists(new_folder_path): os.makedirs(new_folder_path)
                
                create_relax_input(in_path, os.path.join(new_folder_path, f"{new_folder_name}.in"), 
                                   new_cell, atom_names, new_atoms_ang, f"{def_type} {strain*100:.2f}%")
        print(f"  [成功生成] {folder_name} -> 相關應變目錄")

if __name__ == "__main__":
    main()