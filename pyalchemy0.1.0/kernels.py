"""
A module which provides implementations of the kernel of the
Alchemical Integral Transform (AIT) in arbitrary dimensions.

Throughout this code, Hartree atomic units are used.

"""

_reg = 1e-8


# factorial function
def _fc(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n*_fc(n-1)


# binomial coefficient
def _bn(n,k):
    if n == k:
        return 1
    elif k > n:
        return 0
    elif n < 0:
        print('Warning')
    else:
        return _fc(n)/(_fc(k)*_fc(n-k))


# Bell polynomials
def _Bell(n, k, x):
    """
    Return the Bell polynomials of orders n,k where x is a list of length n-k+1.
    x must be at least of length n-k+1, but can be longer, although these elements
    are ignored.
    """
    if n == 0 and k == 0:
        return 1
    elif n > 0 and k == 0:
        return 0
    elif n == 0 and k > 0:
        return 0
    else:
        summe = 0
        for i in range(0, n-k+1):
            # the new x must be of length n-i-k+2
            summe += _bn(n-1, i)*x[i]*_Bell(n-i-1, k-1, x)
        return summe


def _g(i, f, x):
    """
    Return the coefficients of the inverse of the function f whose coefficients
    in return can be computed as its derivatives. f needs to be of form
    f(k, x); k stores the order of the derivative w.r.t. the variable x
    """
    f1 = f(1, x)
    if abs(f1) < _reg:
        # first derivative too small for Lagrange inversion to hold; just return zero,
        # as the number of points where f1 = 0 should be a non-measurable subset
        return 0
    if i == 1:
        return 1/(f1)
    else:
        # initialize all necessary f_i's
        f_hat = [f(k+1, x)/((k+1)*(f1 + _reg)) for k in range(1,i)]
        sum = 0
        for k in range(1,i):
            sum += (-1)**k * (_fc(i+k-1)/_fc(k-1)) * _Bell(i-1,k,f_hat)
        return sum/((f1)**i)


def kernel_1D(v_A, v_B, x, max_order = 4):
    """

    One-dimensional kernel of the Alchemical Integral Transform

    **Parameters:**

    - `v_A` **: callable**
      A scalar function of the initial system's external potential in 1D. It expects two arguments, `k` and `x` such that `v_A(k, x)` $=\frac{\partial^k}{\partial x^k} v_A(x)$
    - `v_B` **: callable**
      A scalar function of the final system's external potential in 1D. It expects two arguments, `k` and `x` such that `v_B(k, x)` $=\frac{\partial^k}{\partial x^k} v_B(x)$
    - `x`**: float**
      coordinate $x$
    - `max_order` **: int, optional**
      Maximum order $p_{max}$ in the kernel to be summed over. Default is 4.
      Caution! Increasing $p_{max}$ also increases the diverging contributions (if there are any) of the integral; higher $p_{max}$ does not necessarily mean higher accuracy!

    **Returns:**

    - **float**
      The 1D kernel of AIT between systems $A$ and $B$ at $x$ up to order $p_{max}$.

    """
    # first order is just Delta v
    summe = v_B(0, x) - v_A(0, x)
    for p in range(2, int(max_order)+1):
        for k in range(1,p):
            def f(y):
                g_vec = [-_g(i, v_A, y) for i in range(1,p-k+1)]
                Delta_v = (v_B(0, y) - v_A(0, y))
                return _Bell(p-1,k,g_vec)*Delta_v**p
            # Approximate derivative with central finite differences
            h = 0.01
            derivative = sum([(-1)**i * _bn(k,i) * f(x+(k/2 - i)*h) for i in range(0,k+1)])/(h**k)
        summe += derivative/_fc(p)
    return summe


def param(v_A, v_B, x, Lambda, max_order=4):
    # zeroth order is just the coordinate
    summe = [x][0]
    Delta_v = v_B(0,x) - v_A(0,x)
    for i in range(1, max_order+1):
        summe += _g(i, v_A, x)*(Delta_v * Lambda)**k/_fc(k)
    return summe
