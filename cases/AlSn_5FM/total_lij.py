"""Total Onsager coefficients and Sn diffusivity in dilute Al-Sn from KineCluE's pair coefficients.

Usage:  python3 total_lij.py [CALC_DIR] [--input FILE] [--name NAME] [--csn X] [--ef E] [--sf S] [--shift]
  CALC_DIR  KineCluE output directory (default CALC/)
  --input   numerical input the run used, for the lattice parameter (default num_AlSn_input.txt; the copy
            KineCluE leaves in CALC_DIR is not refreshed on reruns)
  --name    KineCluE output name (& OUTPUT) to read (default AlSn_5FM)
  --csn     Sn atomic fraction (default 1e-3)
  --ef      vacancy formation energy, eV (default 0.4847, DFT value from the database's Host Information sheet)
  --sf      vacancy formation entropy, kB (default 0)
  --shift   apply the database's self-diffusion correction (A_shift = 12, E_shift = 0.2 eV) to every coefficient
Reads CALC_DIR/NAME_0.dat and writes CALC_DIR/NAME_total.dat.

KineCluE's cluster coefficients l_ij (m2/s) are for one pair, averaged over the Z configurations within
the kinetic range, so the pair contributes C_pair * l_ij with C_pair = C_v * C_Sn * Z (local
equilibrium, C_v = exp(Sf/kB - Ef/kT)). The Z0 vacancy sites that the pair counts are removed from the
free vacancies. To first order in C_Sn, in units of m2/s (L*kT per site):
  L_vv   = C_v * [a^2 w0 + C_Sn * (Z l_vv - Z0 a^2 w0)]
  L_vSn  = C_v * C_Sn * Z * l_vSn
  L_SnSn = C_v * C_Sn * Z * l_SnSn
  D_Sn   = L_SnSn / C_Sn = C_v * Z * l_SnSn
  D_v    = L_vv / C_v_tot, with C_v_tot = C_v * [1 + C_Sn * (Z - Z0)] (free vacancies plus those bound to Sn)
The free-vacancy term a^2 w0 is the uncorrelated vacancy coefficient in fcc ((1/6) * 12 * (a^2/2) * w0); it is
also D_v in pure Al. D_v does not depend on Ef and is never shifted: E_shift corrects self-diffusion as a
whole and cannot be split between vacancy formation and migration.
Species index 0 is the vacancy and 1 is Sn, as in main_AlSn_input.txt.
"""
import argparse
import numpy as np

KB = 1 / 11604.5221  # eV/K, as in kineclue_num.py
A_SHIFT, E_SHIFT = 12, 0.2  # Al row of the database's Host Information sheet

p = argparse.ArgumentParser()
p.add_argument('calc', nargs='?', default='CALC/')
p.add_argument('--input', default='num_AlSn_input.txt')
p.add_argument('--name', default='AlSn_5FM')
p.add_argument('--csn', type=float, default=1e-3)
p.add_argument('--ef', type=float, default=0.4847)
p.add_argument('--sf', type=float, default=0.0)
p.add_argument('--shift', action='store_true')
args = p.parse_args()
calc = args.calc.rstrip('/') + '/'

# Lattice parameter (Angstrom) from the numerical input, w0 prefactor (THz) and barrier (eV) from the jump file
alat = next(float(l.split()[2]) for l in open(args.input) if l.startswith('& LATTPARAM')) * 1e-10
nu0, e0 = next((float(l.split()[1]), float(l.split()[2])) for l in open(calc + 'jump_frequencies.txt')
               if l.rstrip().endswith('# w0'))

d = np.loadtxt(calc + args.name + '_0.dat', ndmin=2)
T, Z, Z0, l_vv, l_vs, l_ss = d[:, 0], d[:, 3], d[:, 4], d[:, 7], d[:, 8], d[:, 9]

cv = np.exp(args.sf - args.ef / (KB * T))
lv_free = alat**2 * nu0 * 1e12 * np.exp(-e0 / (KB * T))
L_vv = cv * (lv_free + args.csn * (Z * l_vv - Z0 * lv_free))
L_vs = cv * args.csn * Z * l_vs
L_ss = cv * args.csn * Z * l_ss
D_sn = cv * Z * l_ss
D_al = cv * 0.7815 * lv_free  # pure-Al self-diffusion, fcc vacancy correlation factor f0 = 0.7815
D_v0 = lv_free  # vacancy diffusivity in pure Al
D_v = L_vv / (cv * (1 + args.csn * (Z - Z0)))  # vacancy diffusivity in the alloy, slowed by trapping on Sn
if args.shift:
    corr = A_SHIFT * np.exp(-E_SHIFT / (KB * T))
    L_vv, L_vs, L_ss, D_sn, D_al = (x * corr for x in (L_vv, L_vs, L_ss, D_sn, D_al))

out = calc + args.name + '_total.dat'
cols = ['1) T [K]', '2) C_v_pure', '3) Ltot_vv', '4) Ltot_vSn', '5) Ltot_SnSn', '6) Ltot_vSn/Ltot_SnSn',
        '7) D_Sn [m2/s]', '8) D_Al [m2/s]', '9) D_v_pure [m2/s]', '10) D_v_alloy [m2/s]']
hdr = ('Dilute Al-Sn, C_Sn = {:g}, Ef = {:g} eV, Sf = {:g} kB, shift = {} (applies to columns 3-8 only)\n'
       'C_v_pure: equilibrium vacancy fraction in pure Al, exp(Sf - Ef/kT)\n'
       'Ltot_*: total (alloy) Onsager coefficients, not KineCluE cluster coefficients; L*kT per lattice site, m2/s\n'
       'D_v_pure: vacancy diffusivity in pure Al; D_v_alloy: vacancy diffusivity including trapping on Sn\n'
       ).format(args.csn, args.ef, args.sf, args.shift) + '{:>9s}'.format(cols[0]) + ''.join(' {:>21s}'.format(c) for c in cols[1:])
np.savetxt(out, np.column_stack([T, cv, L_vv, L_vs, L_ss, L_vs / L_ss, D_sn, D_al, D_v0, D_v]), header=hdr,
           fmt=['%11.1f'] + ['%21.4e'] * 9)
print(open(out).read().rstrip())

# Arrhenius fits, to compare with the database's D0 and Q
for name, D in (('Sn', D_sn), ('Al', D_al), ('v(pure)', D_v0), ('v(alloy)', D_v)):
    slope, icpt = np.polyfit(1 / (KB * T), np.log(D), 1)
    print('D_{}: D0 = {:.4g} cm2/s, Q = {:.4f} eV'.format(name, np.exp(icpt) * 1e4, -slope))
print('Wrote', out)
