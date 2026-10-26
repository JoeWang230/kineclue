"""Fill KineCluE's configurations.txt and jump_frequencies.txt with five-frequency-model values.

Run after kineclue_main.py and before kineclue_num.py:  python3 make_5fm.py CALC/
Energies follow KineCluE's convention: configuration energy = -Eb (Eb > 0 attractive), and the
saddle-point energy is absolute, with the dissociated pair at 0 eV.
"""
import sys
import numpy as np

# ---- Model parameters (eV, THz): Al-Sn row of DFT_solute_diffusivity_all.xlsx (Wu et al. 2016) ----
EM = {         # migration barriers measured from the initial configuration
    'w0': 0.58143,  # bulk Al vacancy jump
    'w1': 0.80046,  # vacancy rotation 1NN -> 1NN around Sn
    'w2': 0.29848,  # Sn-vacancy exchange
    'w3': 0.60063,  # dissociation 1NN -> 2NN/3NN/4NN
    'w4': 0.36456,  # association, reverse of w3
}
NU = {         # attempt frequencies
    'w0': 4.22424,
    'w1': 4.10254,
    'w2': 2.03304,
    'w3': 3.89202,
    'w4': 4.30656,
}
EB = EM['w3'] - EM['w4']  # Sn-vacancy 1NN binding energy (0.236 eV); zero beyond 1NN
# KineCluE divides a jump prefactor by the entropy prefactor of the starting configuration.
# Detailed balance between w3 and w4 needs the 1NN configuration to carry S1NN = nu4/nu3, so jumps
# leaving 1NN are given nu*S1NN to recover the dataset rates (w3/w4 jumps then use nu4).
S1NN = NU['w4'] / NU['w3']
JUMP_NU = {'w0': NU['w0'], 'w1': NU['w1'] * S1NN, 'w2': NU['w2'] * S1NN, 'w3/w4': NU['w4']}
# ------------------------------------------------------------------------

calc = sys.argv[1] if len(sys.argv) > 1 else 'CALC/'
nn1 = np.sqrt(0.5)  # 1NN distance in units of a

# Configurations: find the 1NN class from the vacancy-solute separation
lines = open(calc + 'configurations.txt').read().splitlines()
header, rows, nn1_class = lines[0], [], None
for line in lines[1:]:
    if not line.strip():
        continue
    f = line.split()
    pos = [float(x) for x in line.replace('[', ' ').replace(']', ' ').split()[-6:]]
    dist = np.linalg.norm(np.array(pos[3:]) - np.array(pos[:3]))
    if abs(dist - nn1) < 1e-6 and f[3] == '-':
        nn1_class = f[0]
        line = '{:6s} {:7.4f} {:6.4f}'.format(f[0], S1NN, EB) + line.split(None, 3)[3].join(['  ', ''])
    rows.append(line)
assert nn1_class is not None, '1NN configuration not found'
open(calc + 'configurations.txt', 'w').write('\n'.join([header] + rows) + '\n')

# Jump frequencies: classify each by its end configurations and mechanism
lines = open(calc + 'jump_frequencies.txt').read().splitlines()
header, rows, count = lines[0], [], {}
for line in lines[1:]:
    if not line.strip():
        continue
    f = line.split()
    num, c1, c2, mech, mult = f[0], f[3], f[4], f[5], f[6]
    if mech == 'solute_jump':
        w, esp = 'w2', -EB + EM['w2']
    elif c1 == nn1_class and c2 == nn1_class:
        w, esp = 'w1', -EB + EM['w1']
    elif nn1_class in (c1, c2):
        w, esp = 'w3/w4', -EB + EM['w3']
    else:
        w, esp = 'w0', EM['w0']
    count[w] = count.get(w, 0) + 1
    rows.append('{:6s} {:7.4f} {:7.4f}   {:5s} {:5s} {:14s} {:>4s}   # {}'.format(num, JUMP_NU[w], esp, c1, c2, mech, mult, w))
open(calc + 'jump_frequencies.txt', 'w').write('\n'.join([header] + rows) + '\n')

print('1NN configuration class:', nn1_class)
print('Jump classes assigned:', count)
print('Barriers (eV): ' + ', '.join('{}={:.4f}'.format(k, v) for k, v in EM.items()) + ', Eb={:.4f}, S1NN={:.4f}'.format(EB, S1NN))
