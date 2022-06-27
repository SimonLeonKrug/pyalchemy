# pyalchemy
A library which provides implementations of the kernel of the Alchemical Integral Transform (AIT) in 1D, 2D, 3D. Instead of calculating electronic energies of systems one at a time, this kernel provides a shortcut. By using an initial system's electron density, one can calculate the energy difference to any other system within the radius of convergence of AIT. A complete explanation and introduction of the concept can be found under https://arxiv.org/abs/2203.13794 .

$$\left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)$$

'pyalchemy' does not provide any electron densities, nor the functionology for numerical integration.
