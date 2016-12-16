## Introduction

GeFiCa stands for *Ge detector Field Calculator*. It provides classes to calculate static electric fields and potentials in Ge detectors using various coordinates in up to three dimensions. The field potential values together with their grid coordinates are saved in a [ROOT][] [tree][] to take the advantage of the file compression and the [Draw][] function provided by the [ROOT][] [TTree][] class. GeFiCa is provided as a shared library that can be directly loaded by [ROOT][]. All GeFiCa classes can be used directly in a [ROOT][] interactive session.

The electric potential \f$\varphi\f$ can be calculated using Poisson's equation \f$ \nabla^2 \varphi(\mathbf{x}) = - \rho(\mathbf{x})/\epsilon_0/\epsilon_R \f$, where \f$\rho\f$ is the ionized impurity concentration, \f$\mathbf{x}\f$ denotes the coordinates, and \f$\epsilon_0,\epsilon_R\f$ are the permittivity of vacuum and relative permittivity of Ge, respectively. The equation can be simplified to \f$ \mathrm{d}^2\varphi/\mathrm{d}x^2 = - \rho/\epsilon_0/\epsilon_R \f$ in 1D. The second-order derivative \f$ \mathrm{d}^2\varphi/\mathrm{d}x^2 \f$ can be calculated numerically as \f$ \mathrm{d}^2 \varphi / \mathrm{d} x^2 = [\varphi(x_{i-1}) - 2\varphi(x_i) + \varphi(x_{i+1})]/\Delta x^2 \f$, where \f$ x_i \f$ is the coordinate of the _i_-th point in a grid, and \f$ \Delta x \f$ is the step length of the grid. The numerical calculation of the second-order derivative with higher orders of accuracy can be done using the [table][] of Finite difference coefficient in Wikipedia.

Classes named as a type of Ge detectors, such as Planar1D, TrueCoaxial2D, Sphere, etc., are used to setup boundary conditions for the field calculation. Classes named as a type of coordinates, such as XYZ, RThetaPhi, RhoPhi, etc., are used to update field grid points based on boundary conditions and Poisson's equation. Calculations for different dimensions are split to different classes because of performance considerations - it is slow to use a class designed for 3D calculation to calculate 1D fields. Classes dealing with higher dimensions inherit from lower dimension classes to avoid code duplication. Voltages applied to a detector are saved as public members of class X, which are inherited by all other classes. Constant impurity concentration is saved the same way. Classes named as a type of coordinates are hence not limited to coordinates-related things only, they take care of Ge crystal properties as well. While classes named as a type of Ge detectors takes care of dimensions of electrodes and initialize the boundary conditions accordingly.

## Links

- Doxygen documentation: http://www.physino.xyz/gefica
- Github page: https://github.com/jintonic/gefica

## Directories

Directory | Contents
----------|-----------
core      | source code
docs      | documentation generated using [Doxygen][]
macro     | [ROOT][] macros showing usages of some classes

## Get started

1. Make sure that [ROOT][] is installed.
2. Execute the following commands in a terminal:

~~~{.sh}
git clone https://github.com/jintonic/gefica.git
cd gefica/core
make install
export LD_LIBRARY_PATH=~/lib:$LD_LIBRARY_PATH
cd ../macro
root planar1d.C
~~~

## References

- Overview: https://mediatum.ub.tum.de/node?id=701884
- Field calculation: https://mediatum.ub.tum.de/node?id=969435
- Numerical methods: https://www.mppmu.mpg.de/~jingliu/ECPI/, Lecture 4 and 5

[ROOT]:https://root.cern.ch
[tree]:https://root.cern.ch/root/htmldoc/guides/users-guide/Trees.html
[Draw]:https://root.cern.ch/doc/master/classTTree.html#a73450649dc6e54b5b94516c468523e45
[TTree]:https://root.cern.ch/doc/master/classTTree.html
[Doxygen]:http://www.stack.nl/~dimitri/doxygen/index.html
[table]:https://en.wikipedia.org/wiki/Finite_difference_coefficient
