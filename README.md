# PyAlchemy
A library which provides implementations of the kernel $\space\mathcal{K} $ of the Alchemical Integral Transform (AIT) in 1D, 2D, 3D.

Throughout this README and the code, [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are used.

## Introduction
Instead of calculating electronic energies of systems one at a time, this kernel provides a shortcut. By using an initial system's ($A$) electron density $\rho_A(\vec{r})$, one can calculate the energy difference to any other system ($B$) within the radius of convergence of AIT. A complete explanation and introduction of the concept can be found under https://arxiv.org/abs/2203.13794 .

Consider the two system's $A$ and $B$ with their external potentials $v_A$ and $v_B$. Then their electronic energy difference is given by

$$ E_B - E_A = \int_{\mathbb{R}^3} d\vec{r} \enspace \rho_A \left( \vec{r} \right) \mathcal{K} \left( \vec{r}, v_A, v_B \right) $$

with 3D-kernel

$$ \mathcal{K} \left( \vec{r}, v_A, v_B \right) = \sum_{p = 1}^{\infty} \frac{1}{p} \left(1 - \frac{v_B(\vec{r})}{v_A(\vec{r})} \right)^{p-1} \sum_{S_p} \frac{\partial^{\mu_x + \mu_y + \mu_z} v_B(\vec{r}) - v_A(\vec{r})}{\partial x^{\mu_x} \partial y^{\mu_y} \partial z^{\mu_z}}
    \left[\prod_{i = 1}^{p-1} \frac{ \left( x + y + z \right)^{k_i}}{k_i!} \right]$$

$$S_p := \left\lbrace \mu_x, \mu_y, \mu_z, k_1, \dots, k_{p-1} \in \mathbb{N}\_{0}  \enspace \Bigg\vert \enspace p-1 = \sum_{i=1}^{p-1} i \cdot k_i , \enspace \mu_x + \mu_y + \mu_z = \sum_{i=1}^{p-1} k_i\right\rbrace$$

Analytical expressions for 2D and 1D can be achieved by dropping $\enspace z $- or $\enspace y,z $-dependencies and forcing $\enspace \mu_z = 0 $ or $\enspace \mu_y = \mu_z = 0 $.

Because it is tedious to implement, `pyalchemy` already provides these kernels. It does not provide any electron densities, or functions for numerical integration. Both must be handled with other libraries.

## Documentation

#### `kernel_1D(partial_v_A, partial_v_B, x, orders = [1,2,3], verbose = False)`
    
Parameters:
- `partial_v_A` (_callable_) : a function of the initial system's external potential in 1D. It expects two arguments, $\enspace n_x $ and $\enspace x $, such that `partial_v_A(n_x, x)` $= \frac{\partial^{n_x} }{\partial x^{n_x}} v_A(x) $
- `partial_v_B` (_callable_) : a function of the final system's external potential in 1D. It expects two arguments, $\enspace n_x $ and $\enspace x $, such that `partial_v_B(n_x, x)` $= \frac{\partial^{n_x} }{\partial x^{n_x}} v_B(x) $
- `x` (_float_): coordinate
- `orders` (_list_, _optional_) : a list of the orders $\enspace p $ in the kernel to be summed over. Recommended are at least `[1,2,3]`, precise is `[1,2,3,4,5]`. $\enspace p $ is implemented up to and including 9-th order which is ridiculous overkill.
- `verbose` (_bool_, _optional_) : when `True`, prints a warning if the naive convergence criterion $\enspace |1 - v_B(x)/v_A(x)| < 1 $ is violated. This does not imply divergence of the series but may hint towards too large differences between initial and final system.

Returns:
- `kernel_1D` (_float_) : the 1D kernel of AIT between systems $\enspace A $ and $\enspace B $ at $\enspace x $ for all orders in `orders`.

---
#### `kernel_2D(partial_v_A, partial_v_B, x,y, orders = [1,2,3], verbose = False)`

Parameters:
- `partial_v_A` (_callable_) : a function of the initial system's external potential in 2D. It expects four arguments, $\enspace n_x, n_y $ and $\enspace x, y $, such that `partial_v_A(n_x, n_y, x, y)` $= \frac{\partial^{n_x + n_y} }{\partial x^{n_x} \partial y^{n_y} } v_A(x,y) $
- `partial_v_B` (_callable_) : a function of the final system's external potential in 2D. It expects four arguments, $\enspace n_x, n_y $ and $\enspace x, y $, such that `partial_v_B(n_x, n_y, x, y)` $= \frac{\partial^{n_x + n_y} }{\partial x^{n_x} \partial y^{n_y} } v_B(x,y) $
- `x, y` (_float_): coordinates
- `orders` (_list_, _optional_) : a list of the orders $\enspace p $ in the kernel to be summed over. Recommended are at least `[1,2,3]`, precise is `[1,2,3,4,5]`. $\enspace p $ is implemented up to and including 9-th order which is ridiculous overkill.
- `verbose` (_bool_, _optional_) : when `True`, prints a warning if the naive convergence criterion $\enspace |1 - v_B(x,y)/v_A(x,y)| < 1 $ is violated. This does not imply divergence of the series but may hint towards too large differences between initial and final system.

Returns:
- `kernel_2D` (_float_) : the 2D kernel of AIT between systems $\enspace A $ and $\enspace B $ at $\enspace x, y $ for all orders in `orders`.

---
#### `kernel_3D(partial_v_A, partial_v_B, x,y,z, orders = [1,2,3], verbose = False)`

Parameters:
- `partial_v_A` (_callable_) : a function of the initial system's external potential in 3D. It expects six arguments, $\enspace n_x, n_y, n_z $ and $\enspace x, y, z $, such that `partial_v_A(n_x, n_y, n_z, x, y, z)` $= \frac{\partial^{n_x + n_y + n_z} }{\partial x^{n_x} \partial y^{n_y} \partial z^{n_z} } v_A(x,y,z) $
- `partial_v_B` (_callable_) : a function of the final system's external potential in 3D. It expects six arguments, $\enspace n_x, n_y, n_z $ and $\enspace x, y, z $, such that `partial_v_B(n_x, n_y, n_z, x, y, z)` $= \frac{\partial^{n_x + n_y + n_z} }{\partial x^{n_x} \partial y^{n_y} \partial z^{n_z} } v_B(x,y,z) $
- `x, y, z` (_float_): coordinates
- `orders` (_list_, _optional_) : a list of the orders $\enspace p $ in the kernel to be summed over. Recommended are at least `[1,2,3]`, precise is `[1,2,3,4,5]`. $\enspace p $ is implemented up to and including 9-th order which is ridiculous overkill.
- `verbose` (_bool_, _optional_) : when `True`, prints a warning if the naive convergence criterion $\enspace |1 - v_B(x,y,z)/v_A(x,y,z)| < 1 $ is violated. This does not imply divergence of the series but may hint towards too large differences between initial and final system.

Returns:
- `kernel_3D` (_float_) : the 3D kernel of AIT between systems $\enspace A $ and $\enspace B $ at $\enspace x, y, z $ for all orders in `orders`.

---
#### `partial_v_mol_3D(mole, n_x, n_y, n_z, x, y, z, nuc_rad = 0)`
A built-in function for the external potential in 3D of a given molecule and its spatial derivatives. These derivatives are analytical up to and including third order $\enspace n_x+n_y+n_z \leq 3$, and defined recursively via finite differences for higher orders.

Parameters:
- `mole` (_list_ of _list_) : a list of lists of the 4D coordinates (nuclear charge $\enspace Z $, coordinates $\enspace x, y, z $) of all atoms, i.e. `[[Z1, x1, y1, z1], [Z2, x2, y2, z2], ...]`
- `n_x, n_y, n_z` (_int_) : order of the derivative
- `x, y, z` (_float_): coordinates
- `nuc_rad` (_float_, _optional_) : an optional nuclear radius $\enspace \eta $, such that the Coulomb potential is rendered finite everywhere: $\enspace -Z_1 [(x - x_1)^2 + (y - y_1)^2 + (z - z_1)^2]^{-1/2} \rightarrow -Z_1 [(x - x_1)^2 + (y - y_1)^2 + (z - z_1)^2 + \eta^2]^{-1/2} $

Returns:
- the $\enspace n_x,n_y, n_z$-th derivative of the external potential of `mole` with nuclear radius `nuc_rad` at $\enspace x,y,z $, i.e. $\enspace \frac{\partial^{n_x + n_y + n_z} }{\partial x^{n_x} \partial y^{n_y} \partial z^{n_z} } v_{\text{mole}}(x,y,z) $


## Examples, tricks and functionality
The following examples are descirbed in the SI of the [paper](https://arxiv.org/abs/2203.13794). Except for the periodic systems and complemented by a multi-electron atom example, all codes can be found in the folder `examples`.

#### The hydrogen-like atom in 1D

#### The quantum harmonic oscillator in 1D

#### The Morse-potential in 1D

#### Periodic systems in nD
