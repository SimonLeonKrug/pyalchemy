import numpy as np
from scipy.integrate import romb

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from pyalchemy.kernels import kernel_1D
from pyalchemy.potentials import QHO

# ----------------------------------Parameters----------------------------------
# Use Hartree atomic units throughout!!!

omega_A = 10.0
omega_B_list = [10.2,10.4,11.0,12.0,15]
n_list = [0,1,2,3]






# ------------------------------------------------------------------------------

result_list = []

QHO_A = QHO(omega_A)

for omega_B in omega_B_list:
    QHO_B = QHO(omega_B)

    for n in n_list:

        # The true value
        Delta_E_analyt = QHO_B.E(n) - QHO_A.E(n)

        def v_A(k, x):
            return QHO_A.v(k, x)

        def v_B(k, x):
            return QHO_B.v(k, x)

        def integrand(x):
            return QHO_A.rho(n, x)*kernel_1D(v_A, v_B , x, max_order=10)

        # The value computed via AIT, romb; accurate, but slow
        steps = 2**8 + 1
        limit_low = -30
        limit_high = 30
        dx = (limit_high - limit_low)/steps
        romb_list = [integrand(x) for x in np.linspace(limit_low, limit_high, steps)]
        Delta_E_AIT = romb(romb_list, dx=dx)

        deviation = abs(Delta_E_AIT-Delta_E_analyt)

        # save in list
        result_list.append([n, omega_A, omega_B, Delta_E_AIT, Delta_E_analyt])

# print(result_list)

fs = 16

fig, (ax, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(5,8.5), sharex=True, sharey=False, gridspec_kw={'height_ratios': [5, 3.2]})
fig.subplots_adjust(hspace=0.05)
ax.plot(np.linspace(0.01,1000,100), np.linspace(0.01,1000,100), color='black')

for n in n_list:
    x = [result_list[j][4] for j in range(len(result_list)) if result_list[j][0] == n]
    y = [result_list[j][3] for j in range(len(result_list)) if result_list[j][0] == n]
    deviation = [abs(result_list[j][3]-result_list[j][4]) for j in range(len(result_list)) if result_list[j][0] == n]
    ax.scatter(x,y, marker='+', s = fs*25, label=r'$n = '+str(n)+'$')
    ax2.scatter(x, deviation, marker='x',s = fs*15)

ax.set_ylabel(r'$\Delta E_{\textrm{AIT}} \,\, \textrm{ }[\textrm{Ha}]$', fontsize=fs)
ax.set_xlim([0.05,50])
ax.set_ylim([0.05,50])
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

#plt.savefig('QHO_proof_errors_romb.png',dpi = 300, bbox_inches='tight')
plt.show()
