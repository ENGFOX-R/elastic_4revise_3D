import os
import glob
import stat

# ================= 設定區域 =================
TARGET_ROOT = "/home/engfox/elastic_4revise_3D/relax_inputs/"

# 你的 .sh 模板內容
TEMPLATE_CONTENT = """#!/bin/sh
#SBATCH --output={new_name}.sout
#SBATCH --job-name=MS_QE_Elastic_{new_name}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=All
#source /opt/intel/oneapi/setvars.sh
rm -rf pwscf*
node=1
threads=1
processors=$(nproc)
np=$(($node*$processors/$threads))

export OMP_NUM_THREADS=$threads
export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=AUTO
export LD_LIBRARY_PATH=/opt/mpich-4.0.3/lib:/opt/intel/oneapi/mkl/latest/lib:$LD_LIBRARY_PATH
export PATH=/opt/mpich-4.0.3/bin:$PATH
mpiexec -np $np /opt/thermoPW-7-2_intel/bin/pw.x -in {new_name}.in
rm -rf pwscf*
perl /opt/qe_perl/QEout_analysis.pl
perl /opt/qe_perl/QEout2data.pl
"""

def main():
    if not os.path.exists(TARGET_ROOT):
        print(f"錯誤：找不到路徑 {TARGET_ROOT}")
        return

    generated_scripts = []
    subdirs = glob.glob(os.path.join(TARGET_ROOT, "*"))
    
    for subdir in subdirs:
        if not os.path.isdir(subdir):
            continue
            
        folder_name = os.path.basename(subdir)
        
        valid_keywords = ["scale", "shear", "orig"]
        if not any(key in folder_name for key in valid_keywords):
            continue

        target_in_file = os.path.join(subdir, f"{folder_name}.in")
        
        if not os.path.exists(target_in_file):
            print(f"[Skip] 在 {folder_name} 中找不到對應的 .in 檔 ({folder_name}.in)")
            continue

        sh_filename = f"{folder_name}.sh"
        sh_path = os.path.join(subdir, sh_filename)
        
        new_content = TEMPLATE_CONTENT.format(new_name=folder_name)
        
        with open(sh_path, 'w') as f:
            f.write(new_content)
        
        st = os.stat(sh_path)
        os.chmod(sh_path, st.st_mode | stat.S_IEXEC)
        
        print(f"生成腳本: {sh_path}")
        generated_scripts.append(sh_path)

    if generated_scripts:
        submit_all_path = os.path.join(TARGET_ROOT, "submit_all_elastic.sh")
        with open(submit_all_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 自動提交所有 3D Elastic Constant 相關任務 (ScaleX/Y/Z, ShearYZ/XZ/XY, Orig)\n\n")
            
            for script in generated_scripts:
                folder = os.path.abspath(os.path.dirname(script))
                script_name = os.path.basename(script)
                
                f.write(f"echo \"Submitting job in {folder}...\"\n")
                f.write(f"cd {folder}\n")
                f.write(f"sbatch {script_name}\n")
            
            f.write(f"\ncd {os.path.abspath(TARGET_ROOT)}\n")
            f.write("echo 'All jobs submitted.'\n")
        
        st = os.stat(submit_all_path)
        os.chmod(submit_all_path, st.st_mode | stat.S_IEXEC)
        
        print("-" * 30)
        print(f"完成！共處理 {len(generated_scripts)} 個任務。")
        print(f"請執行以下指令開始計算：")
        print(f"  {submit_all_path}")
    else:
        print("未發現符合條件的資料夾（需包含 scale, shear 或 orig）。")

if __name__ == "__main__":
    main()
