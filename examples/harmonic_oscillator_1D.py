import numpy as np
from pyalchemy import kernel_1D
from scipy.integrate import quad


# ---------------------------------Parameters-----------------------------------
regulator = 25

# excitation state
n = 0

# angular frequencies
omega_A = 13
omega_B = 12


# --------------------setup of Quantum harmonic oscillator----------------------

# Factorial function
def fc(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n*fc(n-1)

# Physicist's Hermite polynomials
def H(n,x):
    if n == 0:
        return 1
    elif n == 1:
        return 2*x
    else:
        return 2*x*H(n-1,x) - 2*(n-1)*H(n-2,x)

# Potential and derivatives of the quantum harmonic oscillator
def partial_v(x, omega, num_der):
    result = 0
    if num_der == 0:
        result = 0.5*(omega**2)*(x**2) + regulator
    elif num_der == 1:
        result = (omega**2)*x
    elif num_der == 2:
        result = omega**2
    else:
        pass
    return result

# External potential of the initial system
def partial_v_A(num_der,x):
    return partial_v(x, omega_A, num_der)

# External potential of the final system
def partial_v_B(num_der,x):
    return partial_v(x, omega_B, num_der)

# Analytical solution of the quantum harmonic oscillator's energy
def Delta_E_analyt():
    return (omega_B - omega_A)*(n + 0.5)

# Electron density of the initial system
def rho_initial(x):
    return (H(n,np.sqrt(omega_A)*x))**2*np.exp(-omega_A*x**2)*np.sqrt(omega_A/np.pi)/(2**n * fc(n))


# -------------------------Alchemical Integral Transform------------------------

# Define an auxiliary integrand function
def integrand(r):
    return rho_initial(r)*kernel_1D(partial_v_A, partial_v_B, r)

# the true value from the analytical solution
print(Delta_E_analyt())

# the predicted value from AIT
print(quad(integrand, -30,30))
