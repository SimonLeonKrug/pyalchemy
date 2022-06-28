# PyAlchemy
A library which provides implementations of the kernel $\space\mathcal{K} $ of the Alchemical Integral Transform (AIT) in 1D, 2D, 3D.

Throughout this README and the code, [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are used.

## Introduction
Instead of calculating electronic energies of systems one at a time, this kernel provides a shortcut. By using an initial system's ($A$) electron density $\rho_A(\vec{r})$, one can calculate the energy difference to any other system ($B$) within the radius of convergence of AIT. A complete explanation and introduction of the concept can be found under https://arxiv.org/abs/2203.13794 .

Consider the two system's $A$ and $B$ with their external potentials $v_A$ and $v_B$. Then their electronic energy difference is given by

$$ E_B - E_A = \int_{\mathbb{R}^3} d\vec{r} \\, \rho_A \left( \vec{r} \right) \mathcal{K} \left( \vec{r}, v_A, v_B \right) $$

with 3D-kernel

$$ \mathcal{K} \left( \vec{r}, v_A, v_B \right) = \sum_{p = 1}^{\infty} \frac{1}{p} \left(1 - \frac{v_B(\vec{r})}{v_A(\vec{r})} \right)^{p-1} \sum_{S_p} \frac{\partial^{\mu_x + \mu_y + \mu_z} v_B(\vec{r}) - v_A(\vec{r})}{\partial x^{\mu_x} \partial y^{\mu_y} \partial z^{\mu_z}}
    \left[\prod_{i = 1}^{p-1} \frac{ \left( x + y + z \right)^{k_i}}{k_i!} \right]$$

$$S_p := \left\lbrace \mu_x, \mu_y, \mu_z, k_1, \dots, k_{p-1} \in \mathbb{N}\_{0}  \\, \Bigg\vert \\, p-1 = \sum_{i=1}^{p-1} i \cdot k_i , \\, \mu_x + \mu_y + \mu_z = \sum_{i=1}^{p-1} k_i\right\rbrace$$

Analytical expressions for 2D and 1D can be achieved by dropping $\\, z $- or $\\, y,z $-dependencies and forcing $\\, \mu_z = 0 $ or $\\, \mu_y = \mu_z = 0 $.

Because it is tedious to implement, `pyalchemy` already provides these kernels. It does not provide any electron densities, or functions for numerical integration. Both must be handled with other libraries.

## Documentation

#### `kernel_1D(partial_v_A, partial_v_B, x, orders = [1,2,3], verbose = False)`

Parameters:
- `partial_v_A` (_callable_) : a function of the initial system's external potential in 1D. It expects two arguments, $\\, n_x $ and $\\, x $, such that `partial_v_A(n_x, x)` $= \frac{\partial^{n_x} }{\partial x^{n_x}} v_A(x) $
- `partial_v_B` (_callable_) : a function of the final system's external potential in 1D. It expects two arguments, $\\, n_x $ and $\\, x $, such that `partial_v_B(n_x, x)` $= \frac{\partial^{n_x} }{\partial x^{n_x}} v_B(x) $
- `x` (_float_): coordinate
- `orders` (_list_, _optional_) : a list of the orders $\\, p $ in the kernel to be summed over. Recommended are at least `[1,2,3]`, precise is `[1,2,3,4,5]`. $\\, p $ is implemented up to and including 9-th order which is ridiculous overkill.
- `verbose` (_bool_, _optional_) : when `True`, prints a warning if the naive convergence criterion $\\, |1 - v_B(x)/v_A(x)| < 1 $ is violated. This does not imply divergence of the series but may hint towards too large differences between initial and final system.

Returns:
- `kernel_1D` (_float_) : the 1D kernel of AIT between systems $\\, A $ and $\\, B $ at $\\, x $ for all orders in `orders`.

---
#### `kernel_2D(partial_v_A, partial_v_B, x,y, orders = [1,2,3], verbose = False)`

Parameters:
- `partial_v_A` (_callable_) : a function of the initial system's external potential in 2D. It expects four arguments, $\\, n_x, n_y $ and $\\, x, y $, such that `partial_v_A(n_x, n_y, x, y)` $= \frac{\partial^{n_x + n_y} }{\partial x^{n_x} \partial y^{n_y} } v_A(x,y) $
- `partial_v_B` (_callable_) : a function of the final system's external potential in 2D. It expects four arguments, $\\, n_x, n_y $ and $\\, x, y $, such that `partial_v_B(n_x, n_y, x, y)` $= \frac{\partial^{n_x + n_y} }{\partial x^{n_x} \partial y^{n_y} } v_B(x,y) $
- `x, y` (_float_): coordinates
- `orders` (_list_, _optional_) : a list of the orders $\\, p $ in the kernel to be summed over. Recommended are at least `[1,2,3]`, precise is `[1,2,3,4,5]`. $\\, p $ is implemented up to and including 9-th order which is ridiculous overkill.
- `verbose` (_bool_, _optional_) : when `True`, prints a warning if the naive convergence criterion $\\, |1 - v_B(x,y)/v_A(x,y)| < 1 $ is violated. This does not imply divergence of the series but may hint towards too large differences between initial and final system.

Returns:
- `kernel_2D` (_float_) : the 2D kernel of AIT between systems $\\, A $ and $\\, B $ at $\\, x, y $ for all orders in `orders`.

---
#### `kernel_3D(partial_v_A, partial_v_B, x,y,z, orders = [1,2,3], verbose = False)`

Parameters:
- `partial_v_A` (_callable_) : a function of the initial system's external potential in 3D. It expects six arguments, $\\, n_x, n_y, n_z $ and $\\, x, y, z $, such that `partial_v_A(n_x, n_y, n_z, x, y, z)` $= \frac{\partial^{n_x + n_y + n_z} }{\partial x^{n_x} \partial y^{n_y} \partial z^{n_z} } v_A(x,y,z) $
- `partial_v_B` (_callable_) : a function of the final system's external potential in 3D. It expects six arguments, $\\, n_x, n_y, n_z $ and $\\, x, y, z $, such that `partial_v_B(n_x, n_y, n_z, x, y, z)` $= \frac{\partial^{n_x + n_y + n_z} }{\partial x^{n_x} \partial y^{n_y} \partial z^{n_z} } v_B(x,y,z) $
- `x, y, z` (_float_): coordinates
- `orders` (_list_, _optional_) : a list of the orders $\\, p $ in the kernel to be summed over. Recommended are at least `[1,2,3]`, precise is `[1,2,3,4,5]`. $\\, p $ is implemented up to and including 9-th order which is ridiculous overkill.
- `verbose` (_bool_, _optional_) : when `True`, prints a warning if the naive convergence criterion $\\, |1 - v_B(x,y,z)/v_A(x,y,z)| < 1 $ is violated. This does not imply divergence of the series but may hint towards too large differences between initial and final system.

Returns:
- `kernel_3D` (_float_) : the 3D kernel of AIT between systems $\\, A $ and $\\, B $ at $\\, x, y, z $ for all orders in `orders`.

---
#### `partial_v_mol_3D(mole, n_x, n_y, n_z, x, y, z, nuc_rad = 0)`
A built-in function for the external potential in 3D of a given molecule and its spatial derivatives. These derivatives are analytical up to and including third order $\\, n_x+n_y+n_z \leq 3$, and defined recursively via finite differences for higher orders.

Parameters:
- `mole` (_list_ of _list_) : a list of lists of the 4D coordinates (nuclear charge $\\, Z $, coordinates $\\, x, y, z $) of all atoms, i.e. `[[Z1, x1, y1, z1], [Z2, x2, y2, z2], ...]`
- `n_x, n_y, n_z` (_int_) : order of the derivative
- `x, y, z` (_float_): coordinates
- `nuc_rad` (_float_, _optional_) : an optional nuclear radius $\\, \eta $, such that the Coulomb potential is rendered finite everywhere: $\\, -Z_1 [(x - x_1)^2 + (y - y_1)^2 + (z - z_1)^2]^{-1/2} \rightarrow -Z_1 [(x - x_1)^2 + (y - y_1)^2 + (z - z_1)^2 + \eta^2]^{-1/2} $

Returns:
- the $\\, n_x,n_y, n_z$-th derivative of the external potential of `mole` with nuclear radius `nuc_rad` at $\\, x,y,z $, i.e. $\\, \frac{\partial^{n_x + n_y + n_z} }{\partial x^{n_x} \partial y^{n_y} \partial z^{n_z} } v_{\text{mole}}(x,y,z) $

---

## Examples, tricks and functionality
The following examples are descirbed in the SI of the [paper](https://arxiv.org/abs/2203.13794). Except for the periodic systems and complemented by a multi-electron atom example, all codes can be found in the folder `examples`.

---
#### The hydrogen-like atom in 1D

Consider the (purely radial) potential of the hydrogen-like atom:

$$ v(r) = -\frac{Z}{r}$$

The eigenenergies are given by

$$ E = -\frac{Z^2}{2n^2}$$

and the spherically-averaged electron density by

$$ \bar{\rho} (r, n, Z) = \frac{1}{n^2} \sum_{l = 0}^{n-1} \frac{2l+1}{4\pi} \left( \frac{2Z}{n} \right)^3 \left( \frac{2Z r}{n} \right)^{2l} \left( L_{n-l-1}^{(2l+1)}\left( \frac{2Zr}{n}\right) \right)^2 \frac{(n-l-1)!}{2n (n+l)!} \exp{\left(- \frac{2Zr}{n}\right)} $$

with generalized Laguerre polynomials $\\, L $. We neglected the change of the reduced mass $\\, \mu $ with increasing nuclear mass and chose $\\, \mu \approx m_e $, hence $\\, a_{\mu} \approx a_0 =1 $ in atomic units.

The radially-averaged hydrogen-like atom can be treated with the kernel in 1D, i.e. we find

$$ \Delta E_{BA} = -\frac{Z_B^2 - Z_A^2}{2n^2} = \int\limits_{0}^{\infty} dr \\, 4\pi r^2 \bar{\rho}_A(r,n,Z_A) \\, \mathcal{K}\left(r, v_A, v_B \right) $$

---
#### The quantum harmonic oscillator in 1D

Consider the potential of the one-dimensional harmonic oscillator

$$ v(x) = \frac{\omega^2}{2}x^2 $$

with eigenenergies
$$ E_n = \omega \\, (n+\frac{1}{2}) $$

and density

$$ \rho (x) = \frac{1}{2^n \\, n!} \sqrt{\frac{\omega}{\pi}} \exp \left(-\omega x^2 \right) \\, \left( H_n \left( \sqrt{\omega} x\right) \right)^2 $$

where $\\, H_n $ are the physicist's Hermite polynomials.

Using AIT to obtain the energy difference between two such systems A and B with frequencies $\\, \omega_A $ and $\\, \omega_B $ proves to be difficult numerically. These numerical difficulties come from the convergence behavior of the series in $\\, \mathcal{K}(x, v_B, v_A) $ and can be evaded by adding a regulatory energy constant $\\, \Lambda_{\text{reg}} \gg \Delta E_{BA} $ to initial and final potential. The energy difference between the systems $\\, \Delta E_{BA} $ and the desnity are unaffected by this but the convergence behavior of the kernel changes towards more favorable regimes.

$$ \Delta E_{BA} = (\omega_B - \omega_A) (n+\frac{1}{2}) = \int\limits_{-\infty}^{+\infty} dx \\, \rho_A(x) \\, \\, \mathcal{K} \left( x, v_A + \Lambda_{\text{reg}} , v_B + \Lambda_{\text{reg}} \right) $$

---
#### The Morse-potential in 1D

Consider the one-dimensional [Morse potential](https://backend.orbit.dtu.dk/ws/portalfiles/portal/3620619/Dahl.pdf) centered around $\\, x_0 $ with well depth $\\, D $ and range parameter $\\, a $:

$$ v(x) = D \\, \left( \exp (-2a (x - x_0)) - 2\exp (-a (x - x_0))\right) $$

with energy eigenvalue

$$ E_n = \sqrt{2D} \\, a \left(n+\frac{1}{2}\right)  \left(1 - \frac{a}{2\sqrt{2D}}\left(n+\frac{1}{2}\right) \right) - D $$

and wave function 

$$ \Psi_n (x) = N(z,n) \sqrt{a} \\, \xi^{z-n-1/2} e^{-\xi/2} L^{(2z-2n-1)}_n (\xi) $$

$$ z = \frac{2D}{a} $$

$$ \xi = 2z\cdot e^{-a(x-x_0)} $$

$$ N(z,n) = \sqrt{\frac{(2z-2n-1) \\, \Gamma (n+1)}{\Gamma (2z-n)}} $$

where $\\, L $ are the generalized Laguerre polynomials.

Again, adding a regulatory constant $\\, \Lambda_{\text{reg}} $ to initial and final potential in the kernel enables us to obtain the energy difference $\\, \Delta E_{BA} $ between to systems $\\, A $ and $\\, B $ with small numerical error.

---
#### Periodic systems in nD

For AIT in periodic systems, one replaces any one-cell-potential by an effective one $\\, v^{\text{eff}}(\vec{r}) $ and uses the borders of the cell $\\, \Omega^n $ as limits of integration:

$$ \Delta E^{\text{cell}}\_{BA} =  \\,\int_{\Omega^n} d\vec{r}\\, \rho_A \left( \vec{r} \right) \\, \\, \mathcal{K} \left( \vec{r}, v^{\text{eff}}_A, v^{\text{eff}}_B \right) $$

This gives the energy difference between two periodic systems per cell.

---
 
