# PyAlchemy
A library which provides implementations of the kernel $\space\mathcal{K} $ of the Alchemical Integral Transform (AIT) in 1D, 2D, 3D.

Throughout this README and the code, [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are used.

The library's code can be found in `src/pyalchemy/main.py`

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

## Examples, tricks and functionality
The following examples are described in the SI of the [paper](https://arxiv.org/abs/2203.13794). Except for the periodic systems and complemented by a multi-electron atom example, all codes can be found in the folder `examples`.

---
#### The hydrogen-like atom in 1D

Consider the (purely radial) potential of the [hydrogen-like atom](https://books.google.at/books?id=BT5RAAAAMAAJ):

$$ v(r) = -\frac{Z}{r}$$

The eigenenergies are given by

$$ E_n = -\frac{Z^2}{2n^2}$$

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

$$ \rho (x,n,\omega) = \frac{1}{2^n \\, n!} \sqrt{\frac{\omega}{\pi}} \exp \left(-\omega x^2 \right) \\, \left( H_n \left( \sqrt{\omega} x\right) \right)^2 $$

where $\\, H_n $ are the physicist's Hermite polynomials.

Using AIT to obtain the energy difference between two such systems A and B with frequencies $\\, \omega_A $ and $\\, \omega_B $ proves to be difficult numerically. These numerical difficulties come from the convergence behavior of the series in $\\, \mathcal{K}(x, v_B, v_A) $ and can be evaded by adding a regulatory energy constant $\\, \Lambda_{\text{reg}} \gg \Delta E_{BA} $ to initial and final potential. The energy difference between the systems $\\, \Delta E_{BA} $ and the desnity are unaffected by this but the convergence behavior of the kernel changes towards more favorable regimes.

$$ \Delta E_{BA} = (\omega_B - \omega_A) (n+\frac{1}{2}) = \int\limits_{-\infty}^{+\infty} dx \\, \rho_A(x) \\, \\, \mathcal{K} \left( x, v_A + \Lambda_{\text{reg}} , v_B + \Lambda_{\text{reg}} \right) $$

---
#### The Morse potential in 1D

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

Again, adding a regulatory constant $\\, \Lambda_{\text{reg}} $ to initial and final potential in the kernel enables us to obtain the energy difference $\\, \Delta E_{BA} $ between two systems $\\, A $ and $\\, B $ with small numerical error.

---
#### Periodic systems in nD

For AIT in periodic systems, one replaces any one-cell-potential by an effective one $\\, v^{\text{eff}}(\vec{r}) $ and uses the borders of the cell $\\, \Omega^n $ as limits of integration:

$$ \Delta E^{\text{cell}}\_{BA} =  \\,\int_{\Omega^n} d\vec{r}\\, \rho_A \left( \vec{r} \right) \\, \\, \mathcal{K} \left( \vec{r}, v^{\text{eff}}_A, v^{\text{eff}}_B \right) $$

This gives the energy difference between two periodic systems per cell.

---
