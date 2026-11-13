# elastic_4revise_3D

**3D Bulk Elastic Constants and Mechanical Properties of High-Entropy Alloys (HEAs)**  
Calculated using Density Functional Theory (DFT) via **Quantum ESPRESSO**.

---

## 1. Overview & Studied Materials

This repository contains the complete calculation inputs, execution scripts, structural models, and post-processing analysis for the 3D elastic stiffness tensors and polycrystalline mechanical properties of two body-centered cubic (BCC) High-Entropy Alloys (HEAs):

1. **Al₀.₅NbTiVZr₀.₅** (Nominal: Al₅Nb₉Ti₉V₉Zr₄, 36-atom $3 \times 3 \times 2$ supercell)
   - Relaxed equilibrium lattice constant: $a_{\text{eq}} = 3.2401\text{ \AA}$ ($V = 612.254\text{ \AA}^3$)
2. **AlNbTiVZr** (Equiatomic: Al₇Nb₈Ti₇V₇Zr₇, 36-atom $3 \times 3 \times 2$ supercell)
   - Relaxed equilibrium lattice constant: $a_{\text{eq}} = 3.2858\text{ \AA}$ ($V = 638.576\text{ \AA}^3$)

---

## 2. Computational Methodology

- **DFT Code**: Quantum ESPRESSO (`pw.x`)
- **Exchange-Correlation Functional**: GGA-PBE (SSSP efficiency pseudopotential library)
- **Cutoffs**: $E_{\text{cut,wfc}} = 40\text{ Ry}$, $E_{\text{cut},\rho} = 320\text{ Ry}$
- **Brillouin Zone Sampling**: $3 \times 3 \times 4$ Monkhorst-Pack k-point mesh
- **Step 1 - vc-relax**: Variable-cell relaxation (`press = 0 kbar`) to obtain stress-free baseline cells.
- **Step 2 - Strain Application**: 6 deformation modes (Scale X, Scale Y, Scale Z, Shear YZ, Shear XZ, Shear XY) applied from $-5\%$ to $+5\%$ with step size $1\%$ (total 61 strained supercells per material).
- **Step 3 - Elastic Stiffness Matrix ($C_{ij}$)**: Extracted via linear regression of the 3D Cauchy stress tensor versus applied strain ($\sigma_{i} = C_{ij} \varepsilon_{j}$).
- **Step 4 - Cubic Projection & VRH Homogenization**:
  - $C_{11}^{\text{cubic}} = \frac{1}{3}(C_{11} + C_{22} + C_{33})$
  - $C_{12}^{\text{cubic}} = \frac{1}{3}(C_{12}^{\text{sym}} + C_{13}^{\text{sym}} + C_{23}^{\text{sym}})$
  - $C_{44}^{\text{cubic}} = \frac{1}{3}(C_{44} + C_{55} + C_{66})$
  - $C' = \frac{1}{2}(C_{11} - C_{12})$, $\quad A = \frac{2 C_{44}}{C_{11} - C_{12}}$ (Zener anisotropy)
  - Voigt-Reuss-Hill (VRH) isotropic polycrystalline averages: Bulk Modulus ($B$), Shear Modulus ($G$), Young's Modulus ($E$), Poisson's Ratio ($\nu$), Pugh ratio ($B/G$), Cauchy Pressure ($C_{12} - C_{44}$), and Born mechanical stability criteria.

---

## 3. Summary of Mechanical Properties (Across 1%–5% Fitting Windows)

### (A) Al₀.₅NbTiVZr₀.₅
| Fit Range | $C_{11}^{\text{cub}}$ (GPa) | $C_{12}^{\text{cub}}$ (GPa) | $C_{44}^{\text{cub}}$ (GPa) | $B_{\text{VRH}}$ (GPa) | $G_{\text{VRH}}$ (GPa) | $E_{\text{VRH}}$ (GPa) | $\nu_{\text{VRH}}$ | Pugh ($B/G$) | Born Stable |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **±1%** | 190.10 | 129.70 | 33.45 | 149.15 | 24.90 | 70.77 | 0.4209 | 5.99 | Yes |
| **±2%** | 164.05 | 119.38 | 37.39 | 134.19 | 28.93 | 80.98 | 0.3994 | 4.64 | Yes |
| **±3%** | 158.54 | 117.16 | 48.05 | 130.79 | 32.81 | 90.84 | 0.3842 | 3.99 | Yes |
| **±4%** | 161.48 | 117.81 | 48.73 | 132.35 | 34.28 | 94.68 | 0.3808 | 3.86 | Yes |
| **±5%** | 159.61 | 117.41 | 47.35 | 131.47 | 33.42 | 92.44 | 0.3828 | 3.93 | Yes |

### (B) AlNbTiVZr
| Fit Range | $C_{11}^{\text{cub}}$ (GPa) | $C_{12}^{\text{cub}}$ (GPa) | $C_{44}^{\text{cub}}$ (GPa) | $B_{\text{VRH}}$ (GPa) | $G_{\text{VRH}}$ (GPa) | $E_{\text{VRH}}$ (GPa) | $\nu_{\text{VRH}}$ | Pugh ($B/G$) | Born Stable |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **±1%** | 143.98 | 106.92 | 54.80 | 119.20 | 34.95 | 95.52 | 0.3664 | 3.41 | Yes |
| **±2%** | 137.48 | 108.34 | 52.69 | 116.84 | 29.74 | 82.24 | 0.3827 | 3.93 | Yes |
| **±3%** | 135.98 | 108.44 | 48.95 | 116.32 | 27.57 | 76.66 | 0.3902 | 4.22 | Yes |
| **±4%** | 137.83 | 108.96 | 51.51 | 118.31 | 29.96 | 82.89 | 0.3832 | 3.95 | Yes |
| **±5%** | 138.51 | 109.31 | 51.34 | 118.95 | 30.33 | 83.85 | 0.3825 | 3.92 | Yes |

---

## 4. Repository Structure

```
elastic_4revise_3D/
├── README.md                                  # Repository documentation
├── .gitignore                                 # Excludes large raw .sout outputs (~5 GB)
├── Al05NbTiVZr05/
│   ├── analyze_elastic_plotver_robust.py      # Regression analysis and plotting tool
│   ├── result/                                # Final 3D figures and summary tables
│   │   ├── Al05NbTiVZr05_Analysis_3D_GPa.png
│   │   └── Al05NbTiVZr05_elastic_summary_3D_GPa.csv
│   └── relax_inputs/                          # 61 strained input subfolders & 1%–5% intervals
├── AlNbTiVZr/
│   ├── SOP.md                                 # Standard Operating Procedure
│   ├── generate_elastic_inputs_crystalVer.py  # Strain generation script
│   ├── generate_sh.py                         # Slurm job generator
│   ├── analyze_elastic_plotver_robust.py      # Elastic constant evaluation script
│   ├── result/                                # Final 3D figures and summary tables
│   │   ├── AlNbTiVZr_Analysis_3D_GPa.png
│   │   └── elastic_summary_3D_GPa.csv
│   └── relax_inputs/                          # 61 strained input subfolders & 1%–5% intervals
```

---

## 5. License & Attribution

Research dataset created by ENGFOX-R. For questions or citation inquiries, please open an issue in this repository.
