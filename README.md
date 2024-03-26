# PyAlchemy

A library which provides implementations of the kernel of the Alchemical Integral Transform (AIT) for general potentials in $n$ dimensions. An introduction to the concept, further explanations and details can be found under https://arxiv.org/abs/2312.04458.

PyAlchemy uses [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units).

This repo includes three versions:

- pyalchemy 0.0.7 which also provides the code and the examples for the first [AIT paper](https://arxiv.org/abs/2203.13794). A full collection is published on [Zenodo](https://zenodo.org/records/10288996), too. This version's documentation can be found in `docs`. Run `firefox pyalchemy0.0.7/docs/_build/html/index.html`. Note that in this version only monoatomic systems have been proven to produce correct results. These shortcomings have been discussed in the [follow-up paper](https://arxiv.org/abs/2312.04458).

- pyalchemy 0.1.0 which includes the code, the examples and a plotting script for Fig. 2 in [the first version of the follow-up paper](https://arxiv.org/abs/2312.04458v1). Note, that this version accidentally ignored the constant term in the first Hohenberg-Kohn theorem when inverting $v_A$ and consequently performs meagerly.

- The current version, pyalchemy 0.2.0, works for all systems in $n$ dimensions if the problem coordinates of the problem statement can be expressed as the coordinates of the final system via an affine transformation. Check out the examples in the [paper](https://arxiv.org/abs/2312.04458). pyalchemy is available on [PyPI](https://pypi.org/project/pyalchemy/). Run `pip install pyalchemy`.

## Introduction

Instead of calculating electronic energies of systems one at a time, the kernel of AIT provides a shortcut. By using an initial system's $A$ electron density $\rho_A(\pmb{r}) $, one can calculate the energy difference and relative electron density to another system $B$ .

`pyalchemy` provides the kernel $\mathcal{K}$ of AIT, as well as potentials, energies and electron densities of select systems. It does not provide functions for numerical integration or methods of computational chemistry such as (post-)Hartree-Fock methods or Density Functional Theory.

## Documentation

---

#### Kernels (`pyalchemy.kernels`)

---

`pyalchemy.kernels.kernel_nD()`

n-dimensional kernel of the Alchemical Integral Transform

**Parameters:**
- `Delta_v` **: callable**
  The difference in external potentials, i.e. $v_B(x) - v_A(x)$ which takes the nD position as argument
- `x` **: array of size n**
  The nD position
- `A` **: callable**
  invertible matrix
- `b` **: callable**
  vector offset
- `rtol` **: float, optional**
  The relative tolerance of the kernel. It determines the number of steps used in the midpoint rule of the $\lambda$-integration

**Returns:**
- **float**
  The kernel in nD at position $x$

---

#### Potentials (`pyalchemy.potentials`)

---

**class** `pyalchemy.potentials.QHO(omega)`

System ''Quantum harmonic oscillator'' (QHO) and its energy $E(n)$, the external potential $v(x)$ and electron density $\rho(n, x)$.

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

- `v(self, x)`

  **Parameters**

  - `x` **: float**
    Coordinate

  **Returns**

  - **float**
    The external potential of the system at coordinate $x$, $v(x) = \frac{\omega^2}{2} x^2$

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

System ''Morse potential'' and its energy $E(n)$, the external potential $v(x)$ and electron density $\rho(n, x)$ as described [here](http://orbit.dtu.dk/files/3620619/Dahl.pdf).

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

- `v(self, x)`

  **Parameters**

  - `x` **: float**
    Coordinate

  **Returns**

  - **float**
    The external potential of the system at coordinate $x$

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

System ''Hydrogen-like atom'' and its energy $E(n)$, the external potential $v(r)$ and electron density $\rho(n, r)$.

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

- `v(self, r)`

  **Parameters**

  - `r` **: float**
    radius, must be greater 0

  **Returns**

  - **float**
    The external potential of the system at radius $r$

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

Any Coulombic (multi-)atomic system in 3D with $N$ nuclei and its external potential $v(\pmb{x})$

**Parameters**

- `mol` **: array of shape (N,4)**
  $N$ 4-vectors of nuclear charge and 3D coordinates, i.e. $\lbrace (Z_1, (\pmb{R}_1)_1, (\pmb{R}_1)_2, (\pmb{R}_1)_3), \\, \dots \rbrace$, e.g. $\text{N}_2$ = `[[7,0,0,0],[7,1.098/0.529,0,0]]`

**Attributes**

- `mol` **: array of shape (N,4)**
  $N$ 4-vectors of nuclear charge and 3D coordinates

**Methods**

- `v(self, x)`

  **Parameters**

  - `x` **: array of shape (3)**
    Coordinate

  **Returns**

  - **float**
    The external potential of the system at coordinate $\pmb{x}$, $v(\pmb{x}) = \displaystyle\sum^N_{i=1} \frac{-Z_i}{|| \pmb{x} - \pmb{R}_i ||_2}$

---
