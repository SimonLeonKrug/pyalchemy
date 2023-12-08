import numpy as np
from scipy.integrate import romb

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from pyalchemy.kernels import kernel_1D
from pyalchemy.potentials import Morse

# ----------------------------------Parameters----------------------------------
# Use Hartree atomic units throughout!!!

D_A = 22
a_A = 1.0
r_A = 0
D_B_list = [22.1, 22.5, 23, 25, 44]
a_B_list = [1.0, 1.2]
r_B_list = [0.0]
n_list = [0,1,2,3]


# ------------------------------------------------------------------------------

result_list = []

Morse_A = Morse(D_A, a_A, r_A)

for D_B in D_B_list:
    for a_B in a_B_list:
        for r_B in r_B_list:
            Morse_B = Morse(D_B, a_B, r_B)

            for n in n_list:

                # The true value
                Delta_E_analyt = Morse_B.E(n) - Morse_A.E(n)

                def v_A(k, x):
                    return Morse_A.v(k, x)

                def v_B(k, x):
                    return Morse_B.v(k, x)

                def integrand(x):
                    return Morse_A.rho(n, x)*kernel_1D(v_A, v_B , x)

                # The value computed via AIT, romb
                steps = 2**13 + 1
                limit_low = -30
                limit_high = 70
                dx = (limit_high - limit_low)/steps
                romb_list = [integrand(x) for x in np.linspace(limit_low, limit_high, steps)]
                Delta_E_AIT = romb(romb_list, dx=dx)


                # save in list
                result_list.append([n, D_A, a_A, r_A, D_B, a_B, r_B, Delta_E_AIT, Delta_E_analyt])

# print(result_list)

fs = 16

fig, (ax, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(5,8.5), sharex=True, sharey=False, gridspec_kw={'height_ratios': [5, 3.2]})
fig.subplots_adjust(hspace=0.05)
ax.plot(np.linspace(0.01,1000,100), np.linspace(0.01,1000,100), color='black')

for n in n_list:
    x = [result_list[j][8] for j in range(len(result_list)) if result_list[j][0] == n]
    y = [result_list[j][7] for j in range(len(result_list)) if result_list[j][0] == n]
    deviation = [abs(result_list[j][7]-result_list[j][8]) for j in range(len(result_list)) if result_list[j][0] == n]
    ax.scatter(x,y, marker='+', s = fs*25, label=r'$n = '+str(n)+'$')
    ax2.scatter(x, deviation, marker='x',s = fs*15)

ax.set_ylabel(r'$\Delta E_{\textrm{AIT}} \,\, \textrm{ }[\textrm{Ha}]$', fontsize=fs)
ax.set_xlim([0.01,15])
ax.set_ylim([0.01,15])
ax.set_xscale('log')
ax.set_yscale('log')
ax.tick_params(axis='y', labelsize=fs)
ax.set_aspect('equal', adjustable='box')
ax.legend(loc='upper left', frameon=False, fontsize=fs)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r'$\Delta E_{\textrm{exact}} \,\, \textrm{ }[\textrm{Ha}]$', fontsize=fs)
ax2.set_ylabel(r'$\Delta \Delta E \,\, \textrm{ }[\textrm{Ha}]$', fontsize=fs)
ax2.tick_params(axis='x', labelsize=fs)
ax2.tick_params(axis='y', labelsize=fs)

#plt.savefig('Morse_proof_errors_romb.png',dpi = 300, bbox_inches='tight')
plt.show()
