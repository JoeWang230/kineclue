# Sn–vacancy pair in fcc Al: five-frequency model with KineCluE

Transport coefficients of a Sn solute–vacancy pair in fcc aluminium, computed with KineCluE using the
Lidiard–LeClaire five-frequency model, 500–900 K.

> **Status: provisional.** Three of the five barriers (w1, w2, w3) are estimates, not published Sn
> values. The model and code are validated (see *Checks*), but the physical numbers will change once
> real DFT barriers are used.

## Results (kinetic range 4)

| T (K) | L_vs/L_ss | Solute correlation f₂ |
|---|---|---|
| 500 | 0.948 | ~0.04 |
| 550 | 0.915 | 0.064 |
| 600 | 0.874 | — |
| 650 | 0.825 | — |
| 700 | 0.771 | 0.118 |
| 750 | 0.712 | — |
| 800 | 0.651 | — |
| 850 | 0.589 | — |
| 900 | 0.527 | 0.195 |

- **Vacancy drag:** L_vs/L_ss is positive across the whole range, so vacancies drag Sn along with
  them. Under irradiation or quenching this would push Sn toward vacancy sinks. The drag weakens as
  temperature rises but doesn't reverse below 900 K.
- **Correlation:** f₂ is small because Sn–vacancy exchange (w2) is faster than the jumps that carry
  the vacancy around the solute. The vacancy keeps swapping back and forth with Sn, and most of those
  swaps cancel out.
- **Absolute Lij** (m²/s) are in `CALC/AlSn_5FM_0.dat`; ratios are in `CALC/AlSn_5FM_0DR.dat`.
  Absolute values are per-cluster quantities and change with kinetic range, so use the ratios for
  physics.

## Checks

**Kinetic range convergence** (drag ratio L_vs/L_ss):

| T (K) | kira 2 | kira 3 | kira 4 | kira 5 |
|---|---|---|---|---|
| 500 | 0.9404 | 0.9467 | 0.9482 | 0.9488 |
| 700 | 0.7440 | 0.7653 | 0.7707 | 0.7728 |
| 900 | 0.4842 | 0.5179 | 0.5268 | 0.5302 |

Going from kinetic range 4 to 5 changes the ratio by less than 0.7%.

**Against Manning's analytical f₂** for the five-frequency model
(f₂ = (2w1 + 7F w3)/(2w1 + 2w2 + 7F w3), with Manning's F(w4/w0)), using KineCluE's L_ss/L_ss⁰:

| T (K) | Manning | kira 2 | kira 3 | kira 4 | kira 5 |
|---|---|---|---|---|---|
| 550 | 0.0639 | 0.0641 | 0.0639 | 0.0639 | 0.0639 |
| 700 | 0.1175 | 0.1184 | 0.1178 | 0.1176 | 0.1175 |
| 900 | 0.1944 | 0.1969 | 0.1952 | 0.1948 | 0.1946 |

KineCluE converges to the analytical result, which confirms the jump assignments are correct.

## Parameters

| Quantity | Value | Source |
|---|---|---|
| Lattice parameter | 4.05 Å | experimental Al |
| Attempt frequency | 10 THz, all jumps | assumed |
| w0 (bulk Al jump) | 0.581 eV | PBE value, Wu et al. 2016 |
| Sn–V binding at 1NN | 0.26 eV | ab initio value quoted in a 2021 study; its own experiment gave 0.32 ± 0.10 eV |
| w1 (rotation 1NN → 1NN) | 0.581 eV | **estimate** (set equal to w0) |
| w2 (Sn–V exchange) | 0.45 eV | **estimate** (Sn assumed a moderately fast diffuser) |
| w3 / w4 (dissociation / association) | 0.711 / 0.451 eV | **estimate**: w0 ± half the binding energy |

Binding is zero beyond 1NN, and every jump not touching the 1NN shell is w0, as in the strict
five-frequency model. f₂ and the drag ratio depend strongly on the barriers, especially w2/w1.
Published Sn-specific w1–w3 values were not found; the Wu/Morgan diffusion database likely has them
in its downloadable dataset.

## Files and how to run

| File | Purpose |
|---|---|
| `main_AlSn_input.txt` | Analytical input: fcc lattice, vacancy + Sn, kinetic range 4 |
| `num_AlSn_input.txt` | Numerical input: 500–900 K, lattice parameter, drag-ratio output |
| `make_5fm.py` | Assigns each of KineCluE's jump classes to w0–w4 by its start and end configurations, and writes the energies into `CALC/configurations.txt` and `CALC/jump_frequencies.txt` |
| `run.sh` | Runs analytical stage → `make_5fm.py` → numerical stage |

To use your own barriers, edit the parameter block at the top of `make_5fm.py` and run `./run.sh`
(about 10 s).

Requires the `sympy-numpy-compat` branch of KineCluE (fixes for current SymPy/SciPy).

Energy convention (KineCluE): configuration energy = −Eb (Eb > 0 is attraction), and saddle-point
energies are absolute, with the dissociated pair at 0 eV. So a barrier from the 1NN configuration
is E_saddle + Eb.

## Sources

- [Wu, Mayeshiba & Morgan, High-throughput ab-initio dilute solute diffusion database (Sci. Data 2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4950524/)
- [Energetics of point defects in aluminum via orbital-free DFT](https://arxiv.org/pdf/1611.07631)
- [First-principles calculations of solute–vacancy interactions in aluminum](https://cpb.iphy.ac.cn/article/2020/2032\cpb_29_6_066103.html)
- [Uni Halle open-data record (vacancy–In/Sn complexes in Al)](https://opendata.uni-halle.de//handle/1981185920/81764)
- [Simonovic & Sluiter, Impurity diffusion activation energies in Al from first principles (PRB 2009)](https://link.aps.org/doi/10.1103/PhysRevB.79.054304)
- [Mantina et al., First principles impurity diffusion coefficients](http://www.mmm.psu.edu/MMantina2009_Acta_DiffusionCoefficient.pdf)
- T. Schuler, L. Messina, M. Nastar, KineCluE, Comput. Mater. Sci. (2019), doi:10.1016/j.commatsci.2019.109191
