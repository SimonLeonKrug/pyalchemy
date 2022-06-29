import numpy as np
from pyalchemy import kernel_1D
from scipy.integrate import quad


# ---------------------------------Parameters-----------------------------------
# excitation state
n = 5

# nuclear charges
Z_A = 9
Z_B = 6


# -------------------------setup of the hydrogen-like atom----------------------

# Factorial function
def fc(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n*fc(n-1)

# Generalized Laguerre polynomials
def L(n, alpha, x):
    if n == 0:
        return 1
    elif n == 1:
        return 1 + alpha - x
    else:
        return ((2*(n - 1) + 1 + alpha - x)*L(n - 1, alpha, x) - (n - 1 + alpha)*L(n - 2, alpha, x))/n

# Potential and derivatives of the hydrogen-like atom
def partial_v(r, Z, num_der):
    return Z*(-1)**(num_der+1)*fc(num_der)/(r**(num_der+1))

# External potential of the initial system
def partial_v_A(num_der, r):
    return partial_v(r, Z_A, num_der)

# External potential of the final system
def partial_v_B(num_der, r):
    return partial_v(r, Z_B, num_der)

# Analytical solution of the hydrogen-like atom
def Delta_E_analyt():
    return -0.5*(Z_B**2 - Z_A**2)/(n**2)

# Spherically averaged electron density of the initial system
def rho_initial(r):
    def av_rho_0_l(Z_0, n, l, r):
        return ((2*l+1)/(4*np.pi))*(2*Z_0/n)**3*((fc(n-l-1))/(2*n*fc(n+l)))*np.exp(-2*Z_0*r/n)*(2*Z_0*r/n)**(2*l)*(L(n-l-1, 2*l+1, 2*Z_0*r/n))**2
    #Sum over all angular momentum states
    return sum([av_rho_0_l(Z_A, n, h, r) for h in range(0,n)])/(n**2)


# -------------------------Alchemical Integral Transform------------------------

# Define an auxiliary integrand function
def integrand(r):
    # Do not forget the Jacobian and the angular integration
    return 4*np.pi*r**2*rho_initial(r)*kernel_1D(partial_v_A, partial_v_B, r)

# the true value from the analytical solution
print(Delta_E_analyt())

# the predicted value from AIT
print(quad(integrand, 1e-6,30))
