"""Plot the vacancy diffusivity in pure Al and in Al-0.04 at.% Sn, 300-800 K.

Usage:  python3 plot_dv.py
Reads CALC/AlSn_5FM_lowT_total.dat (from total_lij.py --name AlSn_5FM_lowT --input num_AlSn_lowT_input.txt
--csn 4e-4). Writes AlSn_5FM_Dv.png: (a) Arrhenius plot of D_v, (b) D_v(alloy)/D_v(pure Al).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
S1, S2 = '#2a78d6', '#eb6834'  # categorical slots 1-2, as in plot_results.py

d = np.loadtxt('CALC/AlSn_5FM_lowT_total.dat')
T, dv_pure, dv_alloy = d[:, 0], d[:, 8], d[:, 9]

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'text.color': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'legend.frameon': False, 'lines.linewidth': 2})
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), facecolor=SURFACE)
for a in ax:
    a.set_facecolor(SURFACE)
    a.set_axisbelow(True)

# (a) Arrhenius plot, temperature on a top axis
a = ax[0]
x = 1000 / T
a.semilogy(x, dv_pure, 'o-', color=S1, ms=3, label='pure Al')
a.semilogy(x, dv_alloy, 's-', color=S2, ms=3, label='Al–0.04 at.% Sn')
i = np.argmin(abs(T - 375))
a.annotate('pure Al', (x[i], dv_pure[i]), xytext=(0, 10), textcoords='offset points', ha='center', color=INK)
a.annotate('Al–0.04 at.% Sn', (x[i], dv_alloy[i]), xytext=(-8, -10), textcoords='offset points', ha='right',
           color=INK)
a.set(xlabel='1000/T (1/K)', ylabel=r'$D_v$ (m$^2$/s)', title='(a) Vacancy diffusivity')
a.legend(loc='upper right')
top = a.secondary_xaxis('top', functions=(lambda v: 1000 / np.maximum(v, 1e-9), lambda t: 1000 / np.maximum(t, 1e-9)))
top.set_xticks([300, 400, 500, 600, 800])
top.set_xlabel('T (K)', color=INK2)
top.spines['top'].set_visible(True)
top.spines['top'].set_color(INK2)

# (b) slowdown from trapping on Sn
a = ax[1]
a.plot(T, dv_alloy / dv_pure, 'o-', color=S2, ms=4)
for t in (300, 500, 800):
    i = np.argmin(abs(T - t))
    a.annotate('{:.2f}'.format(dv_alloy[i] / dv_pure[i]), (T[i], dv_alloy[i] / dv_pure[i]), xytext=(0, 8),
               textcoords='offset points', ha='center', color=INK)
a.set(xlabel='T (K)', ylabel=r'$D_v^{\rm alloy} / D_v^{\rm pure}$', title='(b) Slowdown from vacancy trapping on Sn')
a.set_ylim(0, 1.05)

fig.suptitle('Vacancy diffusivity in Al and Al–0.04 at.% Sn (five-frequency model, DFT values)', color=INK)
fig.tight_layout()
fig.savefig('AlSn_5FM_Dv.png', dpi=200, facecolor=SURFACE)
print('Wrote AlSn_5FM_Dv.png')
