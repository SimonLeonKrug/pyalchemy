"""
A module which provides implementations of the some potentials
​
Throughout this code, Hartree atomic units are used.
​
"""

from scipy.special import gamma
from numpy import sqrt, exp, pi

# Regulator for numerically instable fractions
_reg = 1e-15
float_prec = 18 # guaranteed floating point precision in ciritical steps

# factorial function
def _fc(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n*_fc(n-1)

# Hermite polynomials
def _H(n,x):
    if n == 0:
        return 1
    elif n == 1:
        return 2*x
    else:
        return 2*x*_H(n-1,x) - 2*(n-1)*_H(n-2,x)

# Associated Laguerre polynomials
def _L(n, alpha, x):
    if n == 0:
        return 1
    elif n == 1:
        return 1 + alpha - x
    else:
        return ((2*(n - 1) + 1 + alpha - x)*_L(n - 1, alpha, x) - (n - 1 + alpha)*_L(n - 2, alpha, x))/n


# Built-in class for the 1D quantum harmonic oscillator
class QHO:
    def __init__(self, omega):
        self.omega = omega

    # Return the energy of the QHO
    def E(self, n):
        return (n + 0.5)*self.omega

    # Return the 1D potential of the QHO
    def v(self, x):
        return 0.5*(self.omega*x)**2

    def rho(self, n, x):
        return (_H(n,sqrt(self.omega)*x))**2*exp(-self.omega*x**2)*sqrt(self.omega/pi)/(2**n * _fc(n))


# Built-in class for the 1D Morse potential
class Morse:
    def __init__(self, D, a, r_e):
        self.D = D
        self.a = a
        self.r_e = r_e

    def E(self, n):
        l = sqrt(2*self.D)/self.a
        if n > int(l-0.5):
            print("Eigenenergy with n = " + str(n) + " does not exist!")
            return 0
        else:
            nu = self.a*sqrt(2*self.D)
            return ((n+0.5) - ((n+0.5)**2)/(2*l))*nu

    def v(self, x):
        result = self.D*((-2*self.a)*exp(-2*self.a*(x-self.r_e)) - 2*(-self.a)*exp(-self.a*(x-self.r_e))) + self.D
        return result

    def rho(self, n, x):
        l = sqrt(2*self.D)/self.a
        z = 2*l*exp(-self.a*(x - self.r_e))
        N_squared = _fc(n)*(2*l - 2*n - 1)/(gamma(2*l - n))
        return self.a*N_squared*z**(2*l - 2*n - 1)*exp(-z)*(_L(n,2*l-2*n-1,z))**2


# Built-in function for nD potentials of molecules
class hydlike:
    def __init__(self, Z):
        self.Z = Z

    def E(self, n):
        if n == 0:
            print("Eigenenergy with n = 0 does not exist in the hydrogen-like atom!")
            return 0
        else:
            return -self.Z/(2*n**2)

    def v(self, r):
        if r <= 0:
            print("Only positive radii are allowed in the hydrogen-like atom")
            return 0
        else:
            return -self.Z/r

    def rho(self, n, r):
        xi = 2*self.Z/n
        return sum([(2*l+1)*((xi*r)**(2*l))*(xi**3)*(exp(-xi*r))*(_L(n-l-1,2*l+1,xi*r))**2*_fc(n-l-1)/(2*n*_fc(n+l)) for l in range(0,n)])/(4*pi*n**2)


class Coulomb_3D:
    """
    A class for the external potential in 3D of the given molecule
​
    Parameters:
            mol : array of shape (..., 4)
                A list of lists of the 4D coordinates (nuclear charge :math:`Z_i`, coordinates :math:`x_i, y_i, z_i` of all atoms,
                i.e. ``mole = [[Z_1, x_1, y_1, z_1], [Z_2, x_2, y_2, z_2], ...]``
    """

    def __init__(self, mol):
        self.mol = mol

    def v(self, r):
        """
        A function for the external potential in 3D of the given molecule.
​
        Parameters:
                mol : array of shape (..., 4)
                    A list of lists of the 4D coordinates (nuclear charge :math:`Z_i`, coordinates :math:`x_i, y_i, z_i` of all atoms,
                    i.e. ``mole = [[Z_1, x_1, y_1, z_1], [Z_2, x_2, y_2, z_2], ...]``
                r : list of three floats x, y, z
                    coordinates

        Returns:
                float
                    the external potential at ``\bm{r} = [x,y,z]``,
​
        """
        for i in range(len(self.mol)):
            sum += -self.mol[i][0]/sqrt((x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)
        return sum
