# Sn–vacancy pair in fcc Al: five-frequency model with KineCluE

Transport coefficients of a Sn solute–vacancy pair in fcc aluminium, computed with KineCluE using the
Lidiard–LeClaire five-frequency model, 500–900 K, with DFT barriers and attempt frequencies.

All five barriers and attempt frequencies are the Al–Sn values from the Wu/Morgan DFT dilute-diffusion
database (`../DFT_solute_diffusivity_all.xlsx`, sheet `Al-X`).

## Results (kinetic range 4, thermodynamic range 1)

| T (K) | L_vs/L_ss | Solute correlation f₂ |
|---|---|---|
| 500 | 0.843 | 0.0018 |
| 550 | 0.760 | 0.0035 |
| 600 | 0.662 | 0.0060 |
| 650 | 0.555 | 0.0097 |
| 700 | 0.442 | 0.0147 |
| 750 | 0.328 | 0.0212 |
| 800 | 0.217 | 0.0291 |
| 850 | 0.109 | 0.0387 |
| 900 | 0.007 | 0.0498 |

- **Vacancy drag:** L_vs/L_ss is positive over the whole range, so vacancies drag Sn along with them.
  But it falls steeply with temperature and is close to zero at 900 K. Extrapolating, drag reverses
  at about 905–910 K (≈ 0.97 T_m, with T_m = 933 K), so in practice it holds over the whole
  solid-state range.
- **Correlation:** f₂ is very small (0.002–0.05) because Sn–vacancy exchange (w2 = 0.30 eV) is far
  faster than rotation (w1 = 0.80 eV) and dissociation (w3 = 0.60 eV): w2/w1 ≈ 2000 at 700 K. The
  vacancy swaps back and forth with Sn many times before it leaves, and most of those swaps cancel.
- **Absolute Lij** (m²/s) are in `CALC/AlSn_5FM_0.dat`; ratios are in `CALC/AlSn_5FM_0DR.dat`.
  Absolute values are per-cluster quantities and change with kinetic range, so use the ratios for
  physics.

## Checks

**Thermodynamic range 1 vs full range:** a run with ThRa = KiRa = 4 (no reduction of variables)
gives drag ratios identical to every printed digit. The reduction merges the 2NN, 3NN and 4NN
dissociation jumps into a single w3/w4 frequency, which is the five-frequency model anyway.

**Kinetic range convergence** (drag ratio L_vs/L_ss; kira 5 and 6 from a separate kira-6 run):

| T (K) | kira 2 | kira 3 | kira 4 | kira 5 | kira 6 |
|---|---|---|---|---|---|
| 500 | 0.8157 | 0.8375 | 0.8427 | 0.8447 | 0.8455 |
| 700 | 0.3649 | 0.4269 | 0.4423 | 0.4482 | 0.4506 |
| 800 | 0.1213 | 0.1973 | 0.2165 | 0.2238 | 0.2269 |
| 900 | −0.0984 | −0.0148 | 0.0068 | 0.0150 | 0.0184 |

Kinetic range 4 is within 0.012 of kira 6 at every temperature. Below 800 K that is under 3%.
Near 900 K, where the ratio is close to zero, the relative error is large and the
value still rises with kira. If the exact reversal temperature matters, use kira ≥ 6.

**Against Manning's analytical f₂** for the five-frequency model
(f₂ = (2w1 + 7F w3)/(2w1 + 2w2 + 7F w3), with Manning's F(w4/w0)), using KineCluE's L_ss/L_ss⁰:

| T (K) | Manning | kira 2 | kira 3 | kira 4 |
|---|---|---|---|---|
| 500 | 0.00179 | 0.00181 | 0.00180 | 0.00179 |
| 700 | 0.01468 | 0.01506 | 0.01478 | 0.01472 |
| 900 | 0.04957 | 0.05159 | 0.05016 | 0.04981 |

KineCluE converges to the analytical result, which confirms the jump assignments and prefactors.

**500 K convergence file:** KineCluE's `convergence_analysis` used to skip any temperature where
every |Lij| < 1e-12 m²/s (an absolute tolerance), which left the 500 K file empty. `kinepy.py` now
makes that check relative to the largest Lij. The fixed output matches a rerun in `& UNITS debug`
(Lij ×1e8) to machine precision.

## Parameters

| Jump | Barrier (eV) | Attempt frequency (THz) |
|---|---|---|
| w0 bulk Al | 0.58143 | 4.2242 |
| w1 rotation 1NN → 1NN | 0.80046 | 4.1025 |
| w2 Sn–V exchange | 0.29848 | 2.0330 |
| w3 dissociation 1NN → 2NN/3NN/4NN | 0.60063 | 3.8920 |
| w4 association | 0.36456 | 4.3066 |

Lattice parameter 4.05 Å (experimental Al). Sn–V binding at 1NN is Eb = E3 − E4 = 0.236 eV, zero
beyond 1NN. Every jump not touching the 1NN shell is w0, as in the strict five-frequency model.

The attempt frequencies differ, so the 1NN configuration carries an entropy prefactor
S = ν4/ν3 = 1.1065. Detailed balance between w3 and w4 needs it. KineCluE divides a jump's prefactor
by the starting configuration's S, so `make_5fm.py` writes ν·S for w1 and w2 and ν4 for w3/w4. All
five rates then equal ν·exp(−E/kT) with the dataset values.

The database's own Q for Sn (1.070 eV) is about 0.2 eV above what these frequencies give with the
sheet's vacancy formation energy (0.865 eV). Pure Al shows the same offset (1.266 vs 1.066 eV), so
the database seems to use a different formation energy for Q. It doesn't affect Lij ratios or f₂.

## Files and how to run

| File | Purpose |
|---|---|
| `main_AlSn_input.txt` | Analytical input: fcc lattice, vacancy + Sn, kinetic range 4, thermodynamic range 1 |
| `num_AlSn_input.txt` | Numerical input: 500–900 K, lattice parameter, drag-ratio output |
| `make_5fm.py` | Assigns each of KineCluE's jump classes to w0–w4 by its start and end configurations, and writes the energies into `CALC/configurations.txt` and `CALC/jump_frequencies.txt` |
| `run.sh` | Runs analytical stage → `make_5fm.py` → numerical stage |

To use other barriers or attempt frequencies, edit the parameter block at the top of `make_5fm.py` and run `./run.sh`
(about 10 s).

Requires the `sympy-numpy-compat` branch of KineCluE (fixes for current SymPy/SciPy).

Energy convention (KineCluE): configuration energy = −Eb (Eb > 0 is attraction), and saddle-point
energies are absolute, with the dissociated pair at 0 eV. So a barrier from the 1NN configuration
is E_saddle + Eb. Jump prefactors are in THz.

## Sources

- [Wu, Mayeshiba & Morgan, High-throughput ab-initio dilute solute diffusion database (Sci. Data 2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4950524/)
- [Energetics of point defects in aluminum via orbital-free DFT](https://arxiv.org/pdf/1611.07631)
- [First-principles calculations of solute–vacancy interactions in aluminum](https://cpb.iphy.ac.cn/article/2020/2032\cpb_29_6_066103.html)
- [Uni Halle open-data record (vacancy–In/Sn complexes in Al)](https://opendata.uni-halle.de//handle/1981185920/81764)
- [Simonovic & Sluiter, Impurity diffusion activation energies in Al from first principles (PRB 2009)](https://link.aps.org/doi/10.1103/PhysRevB.79.054304)
- [Mantina et al., First principles impurity diffusion coefficients](http://www.mmm.psu.edu/MMantina2009_Acta_DiffusionCoefficient.pdf)
- T. Schuler, L. Messina, M. Nastar, KineCluE, Comput. Mater. Sci. (2019), doi:10.1016/j.commatsci.2019.109191
