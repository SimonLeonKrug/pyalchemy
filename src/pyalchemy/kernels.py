"""
A module which provides implementations of the kernel of the
Alchemical Integral Transform (AIT) in arbitrary dimensions.

Throughout this code, Hartree atomic units are used.

"""

import numpy as np


def kernel_nD(Delta_v, x, A=None, b=None, rtol=1e-6):
    """
    The kernel of AIT in n dimensions.
​
    Parameters:
            Delta_v : callable
                The difference in external potentials, i.e. :math:`v_B(x) - v_A(x)` which takes the nD position as
                argument
            x : array of size n
                The nD position
            A : callable
                
            b : callable

            rtol : float, optional
                The relative tolerance of the kernel. It determines the number of steps used
                in the midpoint rule of the :math:`\lambda`-integration
​
    Returns:
            float
                the kernel in nD at position :math:`x`
​
    """
    # number of steps in the integration, use the midpoint rule
    # the error of the midpoint rule is at most
    # Error \leq max(f'')*(b-a)^3/24*steps**2, i.e. we can find steps
    # via the rtol demanded above, i.e. rtol := Error/max(f'')
    steps = int(1/np.sqrt(24*rtol))+1
    h = 1/steps
    integral = 0

    # invert the matrix A
    for i in range(0, steps):
        A_inv = np.linalg.inv(A(steps*h))
        new_vec = A_inv @ (x - b(steps*h))
        integral += Delta_v(new_vec)
    return integral
