# 3D (Bulk) 材料彈性常數計算 標準作業流程 (SOP)

本流程包含 5 個步驟，從結構最佳化到生成應變結構、排程計算與最終資料分析。
以下詳細列出各步驟的輸入、輸出資料夾與檔案關聯。

---

### Step 1. vc-relax (尋找零應變參考點)
*   **動作**：對材料進行完全放鬆晶格與原子位置 (vc-relax) 計算。
*   **所在目錄 / 讀取路徑**：`/home/engfox/elastic_4revise_3D/vc-relax_output/`
*   **讀取檔案**：使用者自行準備的原始 QE 輸入檔 (`*.in`)。
*   **生成檔案**：計算完成後的輸出檔 (`*.sout`)，此檔案將作為後續加載應變的無應力平衡態結構。

### Step 2. generate_elastic_inputs_crystalVer.py (產生 3D 應變結構)
*   **動作**：讀取 vc-relax 結果，並套用 6 個方向 (X, Y, Z, ShearYZ, ShearXZ, ShearXY) 的應變矩陣，設定預設應變 -5% 到 5%，每次間隔 1%。
*   **讀取路徑**：`/home/engfox/elastic_4revise_3D/vc-relax_output/` (`INPUT_ROOT`)
*   **讀取檔案**：各個材料子資料夾內的 `*.sout` (提取最佳化後的晶格參數) 與 `*.in` (做為新輸入檔的模板)。
*   **生成路徑**：`/home/engfox/elastic_4revise_3D/relax_inputs/` (`OUTPUT_ROOT`)
*   **生成檔案**：會在 `relax_inputs/` 下為每個材料與每個應變條件產生**專屬的子資料夾** (例如 `材料名-scaleXp01`)，並在該資料夾內生成加上應變的新 QE 輸入檔 `*.in` (設定為 `relax` 計算並開啟應力計算 `tstress = .true.`)。

### Step 3. generate_sh.py (產生 Slurm 排程腳本)
*   **動作**：掃描所有剛生成的應變資料夾，為每個 QE 輸入檔產生對應的 Slurm `.sh` 腳本，並建立一鍵提交程式。
*   **讀取/生成路徑**：`/home/engfox/elastic_4revise_3D/relax_inputs/` (`TARGET_ROOT`)
*   **讀取檔案**：掃描目錄下所有包含 `scale`, `shear` 或 `orig` 的子資料夾內的 `*.in` 檔。
*   **生成檔案**：
    1.  於各子資料夾內生成對應的 Slurm 腳本 `*.sh` (例如 `材料名-scaleXp01.sh`)。
    2.  於 `relax_inputs/` 母目錄下生成總提交腳本 `submit_all_elastic.sh`。

### Step 4. sbatch submit_all_elastic.sh (提交運算)
*   **動作**：執行總提交腳本，透過 slurm 系統將所有應變結構的 QE 運算排入佇列。
*   **執行路徑**：`/home/engfox/elastic_4revise_3D/relax_inputs/`
*   **讀取檔案**：讀取 `submit_all_elastic.sh`，該腳本會自動切換到各個子資料夾執行 `sbatch *.sh`。
*   **生成檔案**：QE 計算執行完畢後，會在各個子資料夾中生成計算結果檔 `*.sout`。

### Step 5. analyze_elastic_plotver.py (計算 3D 彈性常數矩陣與繪圖)
*   **動作**：讀取各應變條件下的計算結果，提取 3x3 應力張量，進行線性擬合並換算為 GPa 單位。
*   **讀取/生成路徑**：`/home/engfox/elastic_4revise_3D/relax_inputs/` (`WORK_DIR`)
*   **讀取檔案**：掃描目錄下所有子資料夾內剛計算完成的 `*.sout` 檔，並從中擷取各方向的 stress (kbar)。
*   **生成檔案**：
    1.  **資料彙整**：在母目錄生成彈性常數總表 `elastic_summary_3D_GPa.csv`，內含 C11, C12, C13, C22, C23, C33, C44, C55, C66 數值。
    2.  **應力應變圖**：在母目錄針對每種材料生成 `材料名_Analysis_3D_GPa.png` 圖檔，顯示 6 個方向的應力-應變擬合線與降伏點預測。
