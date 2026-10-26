"""Plot the Sn-vacancy five-frequency-model results: drag ratio, f2 vs Manning, kinetic-range convergence.

Usage:  python3 plot_results.py [CONV_CALC_DIR]
Reads CALC/ (the kinetic-range-4 run) and, for the convergence data in panels (a) and (c), CALC_kira6/
(the same model run to kinetic range 6). Pass another CALC directory to use it instead; if
CALC_kira6/ is missing, CALC/ is used and the plots stop at kinetic range 4.
Writes AlSn_5FM_results.png.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker

# Five-frequency values, as in make_5fm.py (eV, THz)
E = {'w0': 0.58143, 'w1': 0.80046, 'w2': 0.29848, 'w3': 0.60063, 'w4': 0.36456}
NU = {'w0': 4.2242429813276337, 'w1': 4.1025418572786059, 'w2': 2.0330434324753228, 'w3': 3.8920239816921365,
      'w4': 4.3065647103498712}
KB = 8.617333e-5
TEMPS = np.arange(500, 901, 50)

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
S1, S2, S3 = '#2a78d6', '#eb6834', '#1baf7a'  # categorical slots 1-3


def manning_f2(T):
    w = {k: NU[k] * np.exp(-E[k] / (KB * T)) for k in E}
    x = w['w4'] / w['w0']
    F = 1 - (10*x**4 + 180.5*x**3 + 927*x**2 + 1341*x) / (7*(2*x**4 + 40.2*x**3 + 254*x**2 + 597*x + 436))
    return (2*w['w1'] + 7*F*w['w3']) / (2*w['w1'] + 2*w['w2'] + 7*F*w['w3'])


def load(calc):
    """Return {T: array of rows [kira, drag ratio L_vs/L_ss, f2]} from the convergence study."""
    out = {}
    for T in TEMPS:
        a = np.loadtxt(calc + 'CONVERGENCE_STUDY/S0.0000_T{:d}_D0.dat'.format(T), ndmin=2)
        out[T] = np.column_stack([a[:, 0], a[:, 6] / a[:, 7], a[:, 7] / a[:, 10]])
    return out


def at_kira(data, k, col):
    return np.array([data[T][np.argmin(abs(data[T][:, 0] - k)), col] for T in TEMPS])


main = load('CALC/')
conv_dir = sys.argv[1] if len(sys.argv) > 1 else 'CALC_kira6'
conv = load(conv_dir.rstrip('/') + '/') if os.path.isdir(conv_dir) else main
kmax = conv[TEMPS[0]][:, 0].max()

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'text.color': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'legend.frameon': False, 'lines.linewidth': 2})
fig, ax = plt.subplots(1, 3, figsize=(11, 3.6), facecolor=SURFACE)
for a in ax:
    a.set_facecolor(SURFACE)
    a.set_axisbelow(True)

# (a) drag ratio vs temperature
a = ax[0]
a.axhline(0, color=INK2, lw=0.8)
a.plot(TEMPS, at_kira(main, 4, 1), 'o-', color=S1, ms=5, label='kinetic range 4')
if kmax > 4:
    a.plot(TEMPS, at_kira(conv, kmax, 1), 's--', color=S2, ms=5, label='kinetic range {:g}'.format(kmax))
a.set(xlabel='T (K)', ylabel=r'$L_{vs}/L_{ss}$', title='(a) Vacancy drag ratio')
a.legend(loc='upper right')

# (b) solute correlation factor vs Manning
a = ax[1]
Tf = np.linspace(500, 900, 200)
a.semilogy(Tf, manning_f2(Tf), color=INK2, lw=1.5, label='Manning (analytical)')
a.semilogy(TEMPS, at_kira(main, 4, 2), 'o', color=S1, ms=6, label='KineCluE, kinetic range 4')
a.set(xlabel='T (K)', ylabel=r'$f_2 = L_{ss}/L^0_{ss}$', title=r'(b) Solute correlation factor $f_2$')
a.set_yticks([0.002, 0.005, 0.01, 0.02, 0.05])
a.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g'))
a.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.legend(loc='lower right')

# (c) convergence of the drag ratio with kinetic range
a = ax[2]
for T, c in zip((500, 700, 900), (S1, S2, S3)):
    d = conv[T][np.argsort(conv[T][:, 0])]
    d = d[d[:, 0] >= 1.9]
    a.plot(d[:, 0], d[:, 1], '-', color=c, label='{:d} K'.format(T))
    a.annotate('{:d} K'.format(T), (d[-1, 0], d[-1, 1]), xytext=(4, 0), textcoords='offset points',
               va='center', color=INK)
a.axhline(0, color=INK2, lw=0.8)
a.axvline(4, color=INK2, lw=0.8, ls=':')
a.set(xlabel='kinetic range (a)', ylabel=r'$L_{vs}/L_{ss}$', title='(c) Convergence with kinetic range')
a.set_xlim(right=kmax + 0.8)
a.legend(loc='lower right', bbox_to_anchor=(0.88, 0.12))

fig.suptitle('Sn–vacancy pair in fcc Al, five-frequency model (DFT barriers, Wu et al. 2016)', color=INK)
fig.tight_layout()
fig.savefig('AlSn_5FM_results.png', dpi=200, facecolor=SURFACE)
print('Wrote AlSn_5FM_results.png')
