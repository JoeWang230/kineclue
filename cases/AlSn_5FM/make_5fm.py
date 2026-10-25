"""Fill KineCluE's configurations.txt and jump_frequencies.txt with five-frequency-model values.

Run after kineclue_main.py and before kineclue_num.py:  python3 make_5fm.py CALC/
Energies follow KineCluE's convention: configuration energy = -Eb (Eb > 0 attractive), and the
saddle-point energy is absolute, with the dissociated pair at 0 eV.
"""
import sys
import numpy as np

# ---- Model parameters (eV, THz) - replace with your own DFT values ----
NU = 10.0      # attempt frequency, THz (same for all jumps)
EB = 0.26      # Sn-vacancy 1NN binding energy (positive = attraction); zero beyond 1NN
EM = {         # migration barriers measured from the initial configuration
    'w0': 0.581,  # bulk Al vacancy jump
    'w1': 0.581,  # vacancy rotation 1NN -> 1NN around Sn
    'w2': 0.45,   # Sn-vacancy exchange
    'w3': 0.711,  # dissociation 1NN -> 2NN/3NN/4NN  (w4 = w3 - EB follows from detailed balance)
}
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
        line = '{:6s} {:7.3f} {:6.3f}'.format(f[0], 1.0, EB) + line.split(None, 3)[3].join(['  ', ''])
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
    rows.append('{:6s} {:7.3f} {:7.4f}   {:5s} {:5s} {:14s} {:>4s}   # {}'.format(num, NU, esp, c1, c2, mech, mult, w))
open(calc + 'jump_frequencies.txt', 'w').write('\n'.join([header] + rows) + '\n')

print('1NN configuration class:', nn1_class)
print('Jump classes assigned:', count)
print('Barriers (eV): ' + ', '.join('{}={:.3f}'.format(k, v) for k, v in EM.items()) + ', w4={:.3f}'.format(EM['w3'] - EB))
