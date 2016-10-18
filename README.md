GeFiCa stands for Ge detector Field Calculator. It provides classes to calculate static electric field and potential in various Ge detectors using various coordinates in up to three dimensions. The field potential values together with their grid coordinates are saved in a [ROOT][] [tree][] to take the advantage of the file compression and the [Draw][] function provided by the [ROOT][] [TTree][] class. GeFiCa is provided as a shared library that can be directly loaded by [ROOT][]. All GeFiCa classes can be used directly in a [ROOT][] interactive session.

Directory | Contents
----------|-----------
core      | source code
docs      | documentation generated using [Doxygen][]
macro     | [ROOT][] macros showing usages of some classes

The calculation is done using Poisson's equation:

\f[
\nabla^2\varphi(\mathbf{x}) = - \rho(\mathbf{x})/\epsilon_0/\epsilon_R,
\f]

where \f$\varphi\f$ is the electric potential, \f$\rho\f$ is the ionized impurity concentration, \f$\mathbf{x}\f$ denotes the coordinates, and \f$\epsilon_0,\epsilon_R\f$ are the permittivity of vacuum and relative permittivity of Ge, respectively.

# Documentation

- http://www.physino.xyz/gefica

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
