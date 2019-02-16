- Documentation: https://gefica.readthedocs.io
- Github page: https://github.com/jintonic/gefica

# Basics of numeric solution of Poisson's Equation
The electric potential \f$\varphi\f$ can be calculated using Poisson's equation \f$ \nabla^2 \varphi(\mathbf{x}) = - \rho(\mathbf{x})/\epsilon_0/\epsilon_R \f$, where \f$\rho\f$ is the ionized impurity concentration, \f$\mathbf{x}\f$ denotes the coordinates, and \f$\epsilon_0,\epsilon_R\f$ are the permittivity of vacuum and relative permittivity of Ge, respectively. The equation can be simplified to \f$ \mathrm{d}^2\varphi/\mathrm{d}x^2 = - \rho/\epsilon_0/\epsilon_R \f$ in 1D. The second-order derivative \f$ \mathrm{d}^2\varphi/\mathrm{d}x^2 \f$ can be calculated numerically as \f$ \mathrm{d}^2 \varphi / \mathrm{d} x^2 = [\varphi(x_{i-1}) - 2\varphi(x_i) + \varphi(x_{i+1})]/\Delta x^2 \f$, where \f$ x_i \f$ is the coordinate of the _i_-th point in a grid, and \f$ \Delta x \f$ is the step length of the grid. The numerical calculation of the second-order derivative with higher orders of accuracy can be done using the [table][] of Finite difference coefficient in Wikipedia.

[table]:https://en.wikipedia.org/wiki/Finite_difference_coefficient

# Coding convention

Classes named as a type of Ge detectors, such as Planar1D, TrueCoaxial2D, Sphere, etc., are used to setup boundary conditions for the field calculation. Classes named as a type of coordinates, such as XYZ, RThetaPhi, RhoPhi, etc., are used to update field grid points based on boundary conditions and Poisson's equation. Calculations for different dimensions are split to different classes because of performance considerations - it is slow to use a class designed for 3D calculation to calculate 1D fields. Classes dealing with higher dimensions inherit from lower dimension classes to reuse calculation codes. Voltages applied to a detector are saved as public members of class X, which are inherited by all other classes. Constant impurity concentration is saved the same way. Classes named as a type of coordinates are hence not limited to coordinates-related things only, they take care of Ge crystal properties as well. While classes named as a type of Ge detectors takes care of dimensions of electrodes and initialize the boundary conditions accordingly.

We follow ROOT coding convertion with the following exceptions:

- Classes do not start with `T`, instead, they are all in a namespace called `GeFiCa`
- Public member variables do not start with `f`
