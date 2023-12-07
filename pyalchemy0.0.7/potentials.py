"""
A module which provides implementations of the Coulombic potential in 3D

Throughout this code, Hartree atomic units are used.

"""

# Regulator for numerically instable fractions like v_B/v_A
_reg = 1e-8


# Built-in function for 3D potentials of molecules
def partial_v_mol_3D(mole, n_x, n_y, n_z, x, y, z, nuc_rad = 0):
    """
    A function for the external potential in 3D of a given molecule and its spatial derivatives.
    These derivatives are analytical up to and including third order :math:`n_x+n_y+n_z \leq 3`,
    and defined recursively via finite differences for higher orders.

    Parameters:
            mole : array of shape (..., 4)
                A list of lists of the 4D coordinates (nuclear charge :math:`Z_i`, coordinates :math:`x_i, y_i, z_i` of all atoms,
                i.e. ``mole = [[Z_1, x_1, y_1, z_1], [Z_2, x_2, y_2, z_2], ...]``
            n_x, n_y, n_z : int
                Order of the derivative
            x, y, z : float
                coordinates
            nuc_rad : float, optional
                An optional nuclear radius :math:`\eta` such that the Coulomb potential is rendered finite everywhere:
                :math:`\\frac{-Z_i}{\\sqrt{(x - x_i)^2 + (y - y_i)^2 + (z - z_i)^2}} \\rightarrow \\frac{-Z_i}{\sqrt{(x - x_i)^2 + (y - y_i)^2 + (z - z_i)^2 + \eta^2}}`

    Returns:
            float
                the :math:`n_x+n_y+n_z`-th derivative of the external potential of ``mole``
                with nuclear radius ``nuc_rad`` at ``x,y,z``,
                i.e. :math:`\\frac{\\partial^{n_x + n_y + n_z} }{\\partial x^{n_x} \\partial y^{n_y} \\partial z^{n_z} } v_{\\text{mole}}(x,y,z)`

    """
    sum = 0
    # ----------------------------0-th derivatives------------------------------
    if n_x == 0 and n_y == 0 and n_z == 0:
        for i in range(len(mole)):
            sum += -mole[i][0]/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(1/2) + _reg)
    # ----------------------------1-st derivatives------------------------------
    elif n_x == 1 and n_y == 0 and n_z == 0:
        for i in range(len(mole)):
            sum += mole[i][0]*(x-mole[i][1])/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(3/2) + _reg)
    elif n_x == 0 and n_y == 1 and n_z == 0:
        for i in range(len(mole)):
            sum += mole[i][0]*(y-mole[i][2])/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(3/2) + _reg)
    elif n_x == 0 and n_y == 0 and n_z == 1:
        for i in range(len(mole)):
            sum += mole[i][0]*(z-mole[i][3])/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(3/2) + _reg)
    # ----------------------------2-nd derivatives------------------------------
    elif n_x == 2 and n_y == 0 and n_z == 0:
        for i in range(len(mole)):
            sum += mole[i][0]*(-2*(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(5/2) + _reg)
    elif n_x == 0 and n_y == 2 and n_z == 0:
        for i in range(len(mole)):
            sum += mole[i][0]*(-2*(y-mole[i][2])**2+(x-mole[i][1])**2+(z-mole[i][3])**2)/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(5/2) + _reg)
    elif n_x == 0 and n_y == 0 and n_z == 2:
        for i in range(len(mole)):
            sum += mole[i][0]*(-2*(z-mole[i][3])**2+(x-mole[i][1])**2+(y-mole[i][2])**2)/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(5/2) + _reg)
    elif n_x == 1 and n_y == 1 and n_z == 0:
        for i in range(len(mole)):
            sum += -3*mole[i][0]*(x-mole[i][1])*(y-mole[i][2])/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(5/2) + _reg)
    elif n_x == 1 and n_y == 0 and n_z == 1:
        for i in range(len(mole)):
            sum += -3*mole[i][0]*(x-mole[i][1])*(z-mole[i][3])/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(5/2) + _reg)
    elif n_x == 0 and n_y == 1 and n_z == 1:
        for i in range(len(mole)):
            sum += -3*mole[i][0]*(y-mole[i][2])*(z-mole[i][3])/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(5/2) + _reg)
    # ----------------------------3-rd derivatives------------------------------
    elif n_x == 3 and n_y == 0 and n_z == 0:
        for i in range(len(mole)):
            sum += mole[i][0]*6*((x-mole[i][1])**3 - (x-mole[i][1])*(y-mole[i][2])**2 - (x-mole[i][1])*(z-mole[i][3])**2)/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 0 and n_y == 3 and n_z == 0:
        for i in range(len(mole)):
            sum += mole[i][0]*6*((y-mole[i][2])**3 - (y-mole[i][2])*(x-mole[i][1])**2 - (y-mole[i][2])*(z-mole[i][3])**2)/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 0 and n_y == 0 and n_z == 3:
        for i in range(len(mole)):
            sum += mole[i][0]*6*((z-mole[i][3])**3 - (z-mole[i][3])*(x-mole[i][1])**2 - (z-mole[i][3])*(y-mole[i][2])**2)/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 2 and n_y == 1 and n_z == 0:
        for i in range(len(mole)):
            sum += -mole[i][0]*3*(y - mole[i][2])*( (-4)*(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2 )/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 2 and n_y == 0 and n_z == 1:
        for i in range(len(mole)):
            sum += -mole[i][0]*3*(z - mole[i][3])*( (-4)*(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2 )/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 1 and n_y == 2 and n_z == 0:
        for i in range(len(mole)):
            sum += -mole[i][0]*3*(x - mole[i][1])*( (x-mole[i][1])**2+(-4)*(y-mole[i][2])**2+(z-mole[i][3])**2 )/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 0 and n_y == 2 and n_z == 1:
        for i in range(len(mole)):
            sum += -mole[i][0]*3*(z - mole[i][3])*( (x-mole[i][1])**2+(-4)*(y-mole[i][2])**2+(z-mole[i][3])**2 )/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 1 and n_y == 0 and n_z == 2:
        for i in range(len(mole)):
            sum += -mole[i][0]*3*(x - mole[i][1])*( (x-mole[i][1])**2+(y-mole[i][2])**2+(-4)*(z-mole[i][3])**2 )/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 0 and n_y == 1 and n_z == 2:
        for i in range(len(mole)):
            sum += -mole[i][0]*3*(y - mole[i][2])*( (x-mole[i][1])**2+(y-mole[i][2])**2+(-4)*(z-mole[i][3])**2 )/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    elif n_x == 1 and n_y == 1 and n_z == 1:
        for i in range(len(mole)):
            sum += mole[i][0]*15*( (x-mole[i][1])*(y-mole[i][2])*(z-mole[i][3]) )/((nuc_rad**2+(x-mole[i][1])**2+(y-mole[i][2])**2+(z-mole[i][3])**2)**(7/2) + _reg)
    # --------------------------higher derivatives------------------------------
    elif n_x + n_y + n_z > 3:
        h = 0.01
        if n_x > 3:
            sum += (partial_v_mol_3D(mole,n_x - 1, n_y, n_z, x+h, y, z) - partial_v_mol_3D(mole,n_x - 1, n_y, n_z, x-h, y, z))/(2*h)
        if n_y > 3:
            sum += (partial_v_mol_3D(mole,n_x, n_y - 1, n_z, x, y+h, z) - partial_v_mol_3D(mole,n_x, n_y - 1, n_z, x, y-h, z))/(2*h)
        if n_z > 3:
            sum += (partial_v_mol_3D(mole,n_x, n_y, n_z - 1, x, y, z+h) - partial_v_mol_3D(mole,n_x, n_y, n_z - 1, x, y, z-h))/(2*h)
    # --------------------------------------------------------------------------
    else:
        raise ValueError("The "+str(n_x + n_y + n_z)+"-th derivative is not supported!")
    return sum
