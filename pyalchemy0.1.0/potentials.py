"""
A module which provides implementations of the some potentials

Throughout this code, Hartree atomic units are used.

"""

from scipy.special import gamma

# Regulator for numerically instable fractions
_reg = 1e-8

pi = 3.14159265
e = 2.71828182


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

    # Return the 1D potential of the QHO and its k-th derivative
    def v(self, k, x):
        if k == 0:
            return 0.5*(self.omega*x)**2
        elif k == 1:
            return x*self.omega**2
        elif k == 2:
            return self.omega**2
        else:
            return 0

    def rho(self, n, x):
        return (_H(n,(self.omega**0.5)*x))**2*e**(-self.omega*x**2)*(self.omega/pi)**0.5/(2**n * _fc(n))


# Built-in class for the 1D Morse potential
class Morse:
    def __init__(self, D, a, r_e):
        self.D = D
        self.a = a
        self.r_e = r_e

    def E(self, n):
        l = (2*self.D)**0.5/self.a
        if n > int(l-0.5):
            print("Eigenenergy with n = " + str(n) + " does not exist!")
            return 0
        else:
            nu = self.a*(2*self.D)**0.5
            return ((n+0.5) - ((n+0.5)**2)/(2*l))*nu

    def v(self, k, x):
        result = self.D*((-2*self.a)**k*e**(-2*self.a*(x-self.r_e)) - 2*(-self.a)**k*e**(-self.a*(x-self.r_e)))
        if k == 0:
            result += self.D
        return result

    def rho(self, n, x):
        l = (2*self.D)**0.5/self.a
        z = 2*l*e**(-self.a*(x - self.r_e))
        N_squared = _fc(n)*(2*l - 2*n - 1)/(gamma(2*l - n))
        return self.a*N_squared*z**(2*l - 2*n - 1)*e**(-z)*(_L(n,2*l-2*n-1,z))**2


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

    def v(self, k, r):
        if r <= 0:
            print("Only positive radii are allowed in the hydrogen-like atom")
            return 0
        else:
            return -self.Z*(-1)**k*_fc(k)/(r**(k+1))

    def rho(self, n, r):
        xi = 2*self.Z/n
        return sum([(2*l+1)*((xi*r)**(2*l))*(xi**3)*(e**(-xi*r))*(_L(n-l-1,2*l+1,xi*r))**2*_fc(n-l-1)/(2*n*_fc(n+l)) for l in range(0,n)])/(4*pi*n**2)


class Coulomb_3D:
    """
    A class for the external potential in 3D of the given molecule

    Parameters:
            mol : array of shape (..., 4)
                A list of lists of the 4D coordinates (nuclear charge :math:`Z_i`, coordinates :math:`x_i, y_i, z_i` of all atoms,
                i.e. ``mole = [[Z_1, x_1, y_1, z_1], [Z_2, x_2, y_2, z_2], ...]``
    """

    def __init__(self, mol):
        self.mol = mol

    def v(self, n, r):
        """
        A function for the external potential in 3D of the given molecule and its spatial derivatives.
        These derivatives are analytical up to and including third order :math:`n_x+n_y+n_z \leq 3`,
        and defined recursively via finite differences for higher orders.

        Parameters:
                mol : array of shape (..., 4)
                    A list of lists of the 4D coordinates (nuclear charge :math:`Z_i`, coordinates :math:`x_i, y_i, z_i` of all atoms,
                    i.e. ``mole = [[Z_1, x_1, y_1, z_1], [Z_2, x_2, y_2, z_2], ...]``
                n : list of three ints n_x, n_y, n_z
                    Order of the derivative
                r : list of three floats x, y, z
                    coordinates
                nuc_rad : float, optional
                    An optional nuclear radius :math:`\eta` such that the Coulomb potential is rendered finite everywhere:
                    :math:`\\frac{-Z_i}{\\sqrt{(x - x_i)^2 + (y - y_i)^2 + (z - z_i)^2}} \\rightarrow \\frac{-Z_i}{\sqrt{(x - x_i)^2 + (y - y_i)^2 + (z - z_i)^2 + \eta^2}}`

        Returns:
                float
                    the :math:`\bm{n} = [n_x+n_y+n_z]`-th derivative of the external potential
                    with nuclear radius ``nuc_rad`` at ``\bm{r} = [x,y,z]``,
                    i.e. :math:`\\frac{\\partial^{n_x + n_y + n_z} }{\\partial x^{n_x} \\partial y^{n_y} \\partial z^{n_z} } v_{\\text{mol}}(x,y,z)`

        """
        n[0] = n_x
        n[1] = n_y
        n[2] = n_z
        r[0] = x
        r[1] = y
        r[2] = z
        sum = 0
        # ----------------------------0-th derivatives------------------------------
        if n_x == 0 and n_y == 0 and n_z == 0:
            for i in range(len(self.mol)):
                sum += -self.mol[i][0]/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(1/2) + _reg)
        # ----------------------------1-st derivatives------------------------------
        elif n_x == 1 and n_y == 0 and n_z == 0:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*(x-self.mol[i][1])/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(3/2) + _reg)
        elif n_x == 0 and n_y == 1 and n_z == 0:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*(y-self.mol[i][2])/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(3/2) + _reg)
        elif n_x == 0 and n_y == 0 and n_z == 1:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*(z-self.mol[i][3])/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(3/2) + _reg)
        # ----------------------------2-nd derivatives------------------------------
        elif n_x == 2 and n_y == 0 and n_z == 0:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*(-2*(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(5/2) + _reg)
        elif n_x == 0 and n_y == 2 and n_z == 0:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*(-2*(y-self.mol[i][2])**2+(x-self.mol[i][1])**2+(z-self.mol[i][3])**2)/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(5/2) + _reg)
        elif n_x == 0 and n_y == 0 and n_z == 2:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*(-2*(z-self.mol[i][3])**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2)/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(5/2) + _reg)
        elif n_x == 1 and n_y == 1 and n_z == 0:
            for i in range(len(self.mol)):
                sum += -3*self.mol[i][0]*(x-self.mol[i][1])*(y-self.mol[i][2])/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(5/2) + _reg)
        elif n_x == 1 and n_y == 0 and n_z == 1:
            for i in range(len(self.mol)):
                sum += -3*self.mol[i][0]*(x-self.mol[i][1])*(z-self.mol[i][3])/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(5/2) + _reg)
        elif n_x == 0 and n_y == 1 and n_z == 1:
            for i in range(len(self.mol)):
                sum += -3*self.mol[i][0]*(y-self.mol[i][2])*(z-self.mol[i][3])/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(5/2) + _reg)
        # ----------------------------3-rd derivatives------------------------------
        elif n_x == 3 and n_y == 0 and n_z == 0:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*6*((x-self.mol[i][1])**3 - (x-self.mol[i][1])*(y-self.mol[i][2])**2 - (x-self.mol[i][1])*(z-self.mol[i][3])**2)/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 0 and n_y == 3 and n_z == 0:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*6*((y-self.mol[i][2])**3 - (y-self.mol[i][2])*(x-self.mol[i][1])**2 - (y-self.mol[i][2])*(z-self.mol[i][3])**2)/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 0 and n_y == 0 and n_z == 3:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*6*((z-self.mol[i][3])**3 - (z-self.mol[i][3])*(x-self.mol[i][1])**2 - (z-self.mol[i][3])*(y-self.mol[i][2])**2)/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 2 and n_y == 1 and n_z == 0:
            for i in range(len(self.mol)):
                sum += -self.mol[i][0]*3*(y - self.mol[i][2])*( (-4)*(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2 )/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 2 and n_y == 0 and n_z == 1:
            for i in range(len(self.mol)):
                sum += -self.mol[i][0]*3*(z - self.mol[i][3])*( (-4)*(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2 )/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 1 and n_y == 2 and n_z == 0:
            for i in range(len(self.mol)):
                sum += -self.mol[i][0]*3*(x - self.mol[i][1])*( (x-self.mol[i][1])**2+(-4)*(y-self.mol[i][2])**2+(z-self.mol[i][3])**2 )/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 0 and n_y == 2 and n_z == 1:
            for i in range(len(self.mol)):
                sum += -self.mol[i][0]*3*(z - self.mol[i][3])*( (x-self.mol[i][1])**2+(-4)*(y-self.mol[i][2])**2+(z-self.mol[i][3])**2 )/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 1 and n_y == 0 and n_z == 2:
            for i in range(len(self.mol)):
                sum += -self.mol[i][0]*3*(x - self.mol[i][1])*( (x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(-4)*(z-self.mol[i][3])**2 )/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 0 and n_y == 1 and n_z == 2:
            for i in range(len(self.mol)):
                sum += -self.mol[i][0]*3*(y - self.mol[i][2])*( (x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(-4)*(z-self.mol[i][3])**2 )/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        elif n_x == 1 and n_y == 1 and n_z == 1:
            for i in range(len(self.mol)):
                sum += self.mol[i][0]*15*( (x-self.mol[i][1])*(y-self.mol[i][2])*(z-self.mol[i][3]) )/((nuc_rad**2+(x-self.mol[i][1])**2+(y-self.mol[i][2])**2+(z-self.mol[i][3])**2)**(7/2) + _reg)
        # --------------------------higher derivatives------------------------------
        elif n_x + n_y + n_z > 3:
            h = 0.01
            if n_x > 3:
                sum += (self.v(self.mol,n_x - 1, n_y, n_z, x+h, y, z) - self.v(self.mol,n_x - 1, n_y, n_z, x-h, y, z))/(2*h)
            if n_y > 3:
                sum += (self.v(self.mol,n_x, n_y - 1, n_z, x, y+h, z) - self.v(self.mol,n_x, n_y - 1, n_z, x, y-h, z))/(2*h)
            if n_z > 3:
                sum += (self.v(self.mol,n_x, n_y, n_z - 1, x, y, z+h) - self.v(self.mol,n_x, n_y, n_z - 1, x, y, z-h))/(2*h)
        # --------------------------------------------------------------------------
        else:
            raise ValueError("The "+str(n_x + n_y + n_z)+"-th derivative is not supported!")
        return sum
