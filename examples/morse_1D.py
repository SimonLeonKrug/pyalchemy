import numpy as np
from pyalchemy import kernel_1D
from scipy.integrate import quad
from scipy.special import gamma


# ---------------------------------Parameters-----------------------------------
regulator = 100

# excitation state
n = 0

# equilibrium bond distance
D_initial = 15
D_final = 14

# force constants
a_initial = 1.0
a_final = 1.1

# minimum
r_initial = 0
r_final = -0.05


# ----------------------------setup of Morse potential--------------------------

# maximum excitation guard
max_n = int(np.sqrt(2*D_initial)/a_initial-0.5)
if n > max_n:
    print('Warning: Excitation n = '+str(n)+' exceeds maximum excitation n_max = '+str(max_n))

# General Laguerre polynomials
def L(n, alpha, x):
    if n == 0:
        return 1
    elif n == 1:
        return 1 + alpha - x
    else:
        return ((2*(n - 1) + 1 + alpha - x)*L(n - 1, alpha, x) - (n - 1 + alpha)*L(n - 2, alpha, x))/n

# External potential and derivatives of the Morse potential
def partial_v(D, r_e, a, r, num_der):
    result = D*((-2*a)**(num_der)*np.exp(-2*a*(r-r_e)) - 2*(-a)**(num_der)*np.exp(-a*(r-r_e))) + regulator
    return result

# External potential of the initial system
def partial_v_A(mu_r,r):
    return partial_v(D_initial,r_initial,a_initial,r,mu_r)

# External potential of the final system
def partial_v_B(mu_r,r):
    return partial_v(D_final,r_final,a_final,r,mu_r)

# Analytical solution of the Morse potential's energy
def Delta_E_analyt():
    def E(D,a,n):
        return D*(-a**2/(2*D))*(np.sqrt(2*D)/a - n - 0.5)**2 + 1
    return E(D_final, a_final, n) - E(D_initial, a_initial, n)

# Electron density of the initial system
def rho_initial(r):
    lambda_initial = np.sqrt(2*D_initial)/a_initial
    z_initial = 2*lambda_initial*np.exp(-a_initial*(r - r_initial))
    N_squared = gamma(n+1)*(2*lambda_initial - 2*n - 1)/(gamma(2*lambda_initial - n))
    return N_squared*z_initial**(2*lambda_initial - 2*n - 1)*np.exp(-z_initial)*(L(n,2*lambda_initial-2*n-1,z_initial))**2*(a_initial)


# -------------------------Alchemical Integral Transform------------------------

# Define an auxiliary integrand function
def integrand(r):
    return rho_initial(r)*kernel_1D(partial_v_A, partial_v_B, r)

# the true value from the analytical solution
print(Delta_E_analyt())

# the predicted value from AIT
print(quad(integrand, -20,20))
