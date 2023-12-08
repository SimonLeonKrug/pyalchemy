# PyAlchemy

A library which provides implementations of the kernel $\mathcal{K}$ of the Alchemical Integral Transform (AIT) for general potentials in 1D. An introduction to the concept, further explanations and details can be found under https://arxiv.org/abs/2312.04458.

Throughout this README and the code, [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are used.

This repo includes two versions:

- pyalchemy 0.0.7 (old version) which also provides the code and the examples for the first [AIT paper](https://arxiv.org/abs/2203.13794). A full collection is published on [Zenodo](https://zenodo.org/records/10288996), too. This version's documentation can be found in `docs`. Run `firefox pyalchemy0.0.7/docs/_build/html/index.html`. Note that in this version only monoatomic systems have been proven to produce correct results. These shortcomings have been discussed in the [follow-up paper](https://arxiv.org/abs/2312.04458)

- The current version, pyalchemy 0.1.0, works for all systems as long as the expansion in $\mathcal{K}$ converges. It is available on [PyPI](https://pypi.org/project/pyalchemy/). Run `pip install pyalchemy`.

## Introduction

Instead of calculating electronic energies of systems one at a time, this kernel provides a shortcut. By using an initial system's $A$ electron density $\rho_A(\pmb{r}) $, one can calculate the energy difference to any other system $B$ within the radius of convergence of AIT.

Consider the two system's $A$ and $B$ with their external potentials $v_A$ and $v_B$. Then their electronic energy difference is given by

$$ E_B - E_A = \int_{\mathbb{R}^n} d\pmb{r}_A \\, \rho_A \left( \pmb{r}_A \right) \\, \mathcal{K} \left[ v_A, v_B \right] \left( \pmb{r}_A \right) $$

In 1D, only initial and final potentials $v_A, v_B$ are needed. In nD, the parametrization $\pmb{r}(\lambda)$ is necessary, too. $\pmb{r}(\lambda)$ is a solution of $v_A(\pmb{r}(\lambda)) = (v_B(\pmb{r}_A) - v_A(\pmb{r}_A)) \\, \lambda - v_A(\pmb{r}_A)$

Since this equation can be inverted uniquely for scalar functions (cf. [Lagrange inversion theorem](https://en.wikipedia.org/wiki/Lagrange_inversion_theorem)), no parametrization in the 1D case needs to be provided; the inversion of $v_A$ is handled internally.

`pyalchemy` provides the kernel $\mathcal{K}$ of AIT, as well as potentials, energies and electron densities of select systems. It does not provide functions for numerical integration or methods of computational chemistry such as (post-)Hartree-Fock methods or Density Functional Theory.

## Documentation

---

#### Kernels (`pyalchemy.kernels`)

---

`pyalchemy.kernels.kernel_1D(v_A, v_B, x, max_order = 4)`

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

---

`pyalchemy.kernels.param(v_A, v_B, x, Lambda, max_order=4)`

One-dimensional parametrization $x(\lambda)$ between two systems $A$ and $B$ with external potentials $v_A$ and $v_B$.

**Parameters:**

- `v_A` **: callable**
  A scalar function of the initial system's external potential in 1D. It expects two arguments, `k` and `x` such that `v_A(k, x)` $=\frac{\partial^k}{\partial x^k} v_A(x)$
- `v_B` **: callable**
  A scalar function of the final system's external potential in 1D. It expects two arguments, `k` and `x` such that `v_B(k, x)` $=\frac{\partial^k}{\partial x^k} v_B(x)$
- `x`**: float**
  coordinate $x$
- `Lambda` **: float**
  Interpolation parameter $\lambda$, where $\lambda = 0$ corresponds to system $A$, $\lambda = 1$ corresponds to system $B$
- `max_order` **: int, optional**
  Maximum order $p_{max}$ after which the inverted series will be truncated

**Returns:**

- **float**
  The 1D parametrization of AIT between systems $A$ and $B$ at $x$ up to order $p_{max}$.

---

#### Potentials (`pyalchemy.potentials`)

---

**class** `pyalchemy.potentials.QHO(omega)`

System ''Quantum harmonic oscillator'' (QHO) and its energy $E(n)$, $k$-th derivative of the external potential $v(x)$ and electron density $\rho(n, x)$.

**Parameters**

- `omega` **: float**
  Frequency $\omega$ of the QHO

**Attributes**

- `omega` **: float**
  Frequency $\omega$ of the QHO

**Methods**

- `E(self, n)`

  **Parameters**

  - `n` **: int**
    Number of the excited state, $n = \lbrace 0,1,2, \dots \rbrace$

  **Returns**

  - **float**
    The $n$-th eigenenergy of the system, $E = (n + 1/2) \omega$

- `v(self, k, x)`

  **Parameters**

  - `k` **: int**
    Order of spatial derivative

  - `x` **: float**
    Coordinate

  **Returns**

  - **float**
    The external potential of the system at coordinate $x$, $v(x) = \frac{\omega^2}{2} x^2$, and its $k$-th derivative

- `rho(self, n, x)`

  **Parameters**

  - `n` **: int**
    Number of the excited state, $n = \lbrace 0,1,2, \dots \rbrace$

  - `x` **: float**
    Coordinate

  **Returns**

  - **float**
    The electron density $\rho$ of the $n$-th excited state of the system at coordinate $x$

---

**class** `pyalchemy.potentials.Morse(D, a, r_e)`

System ''Morse potential'' and its energy $E(n)$, $k$-th derivative of the external potential $v(x)$ and electron density $\rho(n, x)$ as described [here](http://orbit.dtu.dk/files/3620619/Dahl.pdf).

**Parameters**

- `D` **: float**
  The depth of the well

- `a` **:float**
  The width of the well

- `r_e` **: float**
  The equilibrium distance

**Attributes**

- `D` **: float**
  The depth of the well

- `a` **:float**
  The width of the well

- `r_e` **: float**
  The equilibrium distance

**Methods**

- `E(self, n)`

  **Parameters**

  - `n` **: int**
    Number of the excited state, $n = \lbrace 0,1,2, \dots, \lfloor \frac{\sqrt{2D}}{a} - \frac{1}{2} \rfloor \rbrace$

  **Returns**

  - **float**
    The $n$-th eigenenergy of the system, $E = \frac{4D}{a^2}(n + 1/2) - (n + 1/2)^2$

- `v(self, k, x)`

  **Parameters**

  - `k` **: int**
    Order of spatial derivative

  - `x` **: float**
    Coordinate

    **Returns**

  - **float**
    The external potential of the system at coordinate $x$ and its $k$-th derivative

- `rho(self, n, x)`

  **Parameters**

  - `n` **: int**
    Number of the excited state, $n = \lbrace 0,1,2, \dots, \lfloor \frac{\sqrt{2D}}{a} - \frac{1}{2} \rfloor \rbrace$

  - `x` **: float**
    Coordinate

  **Returns**

  - **float**
    The electron density $\rho$ of the $n$-th excited state of the system at coordinate $x$

---

**class** `pyalchemy.potentials.hydlike(Z)`

System ''Hydrogen-like atom'' and its energy $E(n)$, $k$-th derivative of the external potential $v(r)$ and electron density $\rho(n, r)$.

**Parameters**

- `Z` **: float**
  Nuclear charge $Z$ of the atom

**Attributes**

- `Z` **: float**
  Nuclear charge $Z$ of the atom

**Methods**

- `E(self, n)`

  **Parameters**

  - `n` **: int**
    Number of the excited state, $n = \lbrace 0,1,2, \dots \rbrace$

  **Returns**

  - **float**
    The $n$-th eigenenergy of the system, $E = -\frac{Z^2}{2n^2}$

- `v(self, k, r)`

  **Parameters**

  - `k` **: int**
    Order of spatial derivative

  - `r` **: float**
    radius, must be greater 0

  **Returns**

  - **float**
    The $k$-th derivative of the external potential of the system at radius $r$, $\frac{\partial ^k}{\partial r^k} v(r) = -(-1)^k \\, k! \\, \frac{Z}{r^{k+1}}$

- `rho(self, n, r)`

  **Parameters**

  - `n` **: int**
    Number of the excited state, $n = \lbrace 0,1,2, \dots \rbrace$

  - `r` **: float**
    radius, must be greater 0

  **Returns**

  - **float**
    The electron density $\rho$ of the $n$-th excited state of the system at radius $r$

---

**class** `pyalchemy.potentials.Coulomb_3D(mol)`

Any Coulombic (multi-)atomic system in 3D with $N$ nuclei and its $\pmb{k}$-th derivative of the external potential $v(\pmb{x})$

**Parameters**

- `mol` **: array of shape (N,4)**
  $N$ 4-vectors of nuclear charge and 3D coordinates, , i.e. $\lbrace (Z_1, (\pmb{R}_1)_1, (\pmb{R}_1)_2, (\pmb{R}_1)_3), \\, \dots \rbrace$, e.g. $\text{N}_2$ = `[[7,0,0,0],[7,1.098/0.529,0,0]]`

**Attributes**

- `mol` **: array of shape (N,4)**
  $N$ 4-vectors of nuclear charge and 3D coordinates

**Methods**

- `v(self, k, x)`

  **Parameters**

  - `k` **: array of shape (3)**
    Integer order of spatial derivatives

  - `x` **: array of shape (3)**
    Coordinate

  **Returns**

  - **float**
    The external potential of the system at coordinate $\pmb{x}$, $v(\pmb{x}) = \displaystyle\sum^N_{i=1} \frac{-Z_i}{|| \pmb{x} - \pmb{R}_i ||_2}$, and its $\pmb{k}$-th spatial derivative. All derivatives $|\pmb{k}| < 4$ are defined analytically, all higher derivatives are computed recursively from the lower ones via central finite differences.

---
