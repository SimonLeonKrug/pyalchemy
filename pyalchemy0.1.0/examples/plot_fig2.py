import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
plt.rcParams['text.usetex'] = True
from scipy.optimize import curve_fit

#fontsize
fs= 16
anno_fs = 8

# Dict for the elements
elements = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10}
inv_elements = {v: k for k, v in elements.items()}


# A function to get the right exponent for charged atoms
def exponent(DeltaZ):
    if DeltaZ == 0:
        return ''
    elif DeltaZ == 1:
        return '+'
    elif DeltaZ == -1:
        return '-'
    else:
        return str(int(abs(DeltaZ))) + '{0:+}'.format(DeltaZ)[0]

# ------------------------------------------------------------------------------
# The data setup

# List of number of electrons to consider
N_e_list = [1,2,3,4,5,6,7,8,9]

# Open the right file
# Note that the nuclear charge of the initial atom equals the number of electrons!
# The data can be found under https://zenodo.org/records/6779769
df = pd.read_csv('multielectron_atom_H_to_Ne.txt', usecols=['initial_charge', 'final_charge', 'DeltaE_SCF'])

# All the initial Z needed
Z_A = []
# All the c's that are fitted below
c_list = []

# ------------------------------------------------------------------------------
# The first plot: DeltaE vs. the nuclear charge of the final atom

fig, ax = plt.subplots()
for N_e in N_e_list:

    # Get the data:
    if N_e == 1:
        # If number of electrons N_e is 1, use analytical expression of the
        # hydrogen-like atom
        Z_B = [1,2,3,4,5,6,7,8,9]
        DeltaE = [-0.5*(Z**2 - 1**2) for Z in Z_B]
    else:
        # Remember: Z_A = N_e!
        df_DeltaE = df[df['initial_charge'] == N_e]
        # Remove all lines where the final charge is no integer
        Z_B = [df_DeltaE['final_charge'].tolist()[i] for i in range(len(df_DeltaE['DeltaE_SCF'].tolist())) if df_DeltaE['final_charge'].tolist()[i] == int(df_DeltaE['final_charge'].tolist()[i])]
        DeltaE = [df_DeltaE['DeltaE_SCF'].tolist()[i] for i in range(len(df_DeltaE['DeltaE_SCF'].tolist())) if df_DeltaE['final_charge'].tolist()[i] == int(df_DeltaE['final_charge'].tolist()[i])]
        # Add the zeros, i.e. Z_A = Z_B
        Z_B.append(N_e)
        DeltaE.append(0)

        # Annotate a few points as examples
        if N_e >= 2:
            for i in range(len(Z_B)):
                if Z_B[i]-N_e > 1:
                    continue
                else:
                    ax.annotate(r'$\textrm{' +str(inv_elements[int(Z_B[i])]) +  '}^{' + exponent(Z_B[i]-N_e) + '}$', xy=(Z_B[i]+0.1,DeltaE[i]), va='bottom', ha='left', fontsize=anno_fs)

    # Plot the data:
    ax.scatter(Z_B, DeltaE, label=r'$N_e = ' + str(N_e) + '$', marker='2', s=100)

    # Fit the data:
    def fit_DeltaE_to_Z(Z,a,c):
        return a*Z**2 + c

    popt, pcov = curve_fit(fit_DeltaE_to_Z, Z_B, DeltaE)
    perr = np.sqrt(np.diag(pcov))
    residuals = [DeltaE[i] - fit_DeltaE_to_Z(Z_B[i], *popt) for i in range(len(DeltaE))]
    ss_res = np.sum([res**2 for res in residuals])
    ss_tot = np.sum((np.array(DeltaE)-np.mean(DeltaE))**2)
    r_squared = 1 - (ss_res / ss_tot)

    print('N_e = '+str(N_e))
    print('a,c = ', *popt)
    # print(pcov)
    print('Error = ', *perr)
    print('R_squared = ', r_squared)
    print()

    # Annotate the data
    if N_e == 1:
        va = 'top'
    else:
        va = 'center'
    ax.annotate(r'$N_e = ' + str(N_e) + '$', xy=(0.7, fit_DeltaE_to_Z(0.8, *popt)), va=va, ha='right')


    # Plot the fit:
    Z_B_fit = np.linspace(0.8, 11, 1000)
    DeltaE_fit = [fit_DeltaE_to_Z(Z, *popt) for Z in Z_B_fit]
    plt.plot(Z_B_fit, DeltaE_fit)

    # Annotate the figure
    ax.annotate(r'$\Delta E = a Z_B^2 + c$', xy=(0.5,-20), fontsize=fs)

    # Add and annotate a zero line
    plt.plot(np.linspace(1,12,100),np.zeros(100), ls='dashed', lw=0.7, color='black')
    plt.annotate(r'$\Delta Z = 0$', xy=(10,0), va='bottom', ha='center', fontsize=anno_fs)

    # Gather data for the fit of the parameter c vs. N_e
    c_list.append(popt[1])


# Parameters for the plot
# ax.legend(loc='upper right', frameon=False, fontsize=fs)
ax.set_xlim([-0.3,10.5])
ax.set_ylim([-30,100])
ax.set_xlabel(r'$Z_B$', fontsize=fs)
ax.set_ylabel(r'$\Delta E = E(Z_B, N_e) - E(Z_{B^\prime} = N_e, N_e) \,\, \textrm{[Ha]}$', fontsize=fs)
ax.set_xticks([1,2,3,4,5,6,7,8,9,10],[r'$1$',r'$2$',r'$3$',r'$4$',r'$5$',r'$6$',r'$7$',r'$8$',r'$9$',r'$10$'])
ax.tick_params(axis='both', labelsize=fs)
ax_twin = ax.twiny()
ax_twin.set_xlim([-0.3,10.5])
ax_twin.set_xticks([1,2,3,4,5,6,7,8,9,10],[r'$\textrm{H}$',r'$\textrm{He}$',r'$\textrm{Li}$',r'$\textrm{Be}$',r'$\textrm{B}$',r'$\textrm{C}$',r'$\textrm{N}$',r'$\textrm{O}$',r'$\textrm{F}$',r'$\textrm{Ne}$'], fontsize=fs)

fig.tight_layout()
fig.savefig('DeltaE_vs_Z_quadratic.png', dpi=300)


# ------------------------------------------------------------------------------
# The second plot: parameter c vs the electron number N_e

# Plot the data
fig2, ax2 = plt.subplots()
ax2.scatter(N_e_list, c_list, marker='+', s=300, color='black')

# Fit the data
def fit_c_to_num_elec(N_e, gamma):
    return gamma*N_e**(7/3)

popt_c, pcov_c = curve_fit(fit_c_to_num_elec, N_e_list, c_list)
perr_c = np.sqrt(np.diag(pcov_c))
residuals_c = [c_list[i] - fit_c_to_num_elec(N_e_list[i], *popt_c) for i in range(len(c_list))]
ss_res_c = np.sum([res**2 for res in residuals_c])
ss_tot_c = np.sum((np.array(c_list)-np.mean(c_list))**2)
r_squared_c = 1 - (ss_res_c / ss_tot_c)

print('--------------')
print('gamma = ', popt_c[0])
# print(pcov_c)
print('Error = ', perr_c[0])
print('R_squared = ', r_squared_c)
print()

# Plot the fit
N_e_fit = np.linspace(1,9,100)
c_fit = [fit_c_to_num_elec(num, popt_c[0]) for num in N_e_fit]
plt.plot(N_e_fit,c_fit, ls='solid', lw=3, color='black')

# Annotate the figure
ax2.annotate(r'$c = \gamma N_e^{7/3}$', xy=(1.3,60), fontsize=60)

ax2.set_xlabel(r'$N_e$', fontsize=3*fs)
ax2.set_ylabel(r'$c \,\, \textrm{[Ha]}$', fontsize=2*fs)
ax2.set_xticks([1,2,3,4,5,6,7,8,9],[r'$1$',r'$2$',r'$3$',r'$4$',r'$5$',r'$6$',r'$7$',r'$8$',r'$9$'])
ax2.tick_params(axis='both', labelsize=3*fs)
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(3)
fig2.tight_layout()
fig2.savefig('c_fit.png', dpi=300, transparent=True)
