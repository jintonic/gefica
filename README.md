[![Doxygen](https://codedocs.xyz/jintonic/gefica.svg)](https://codedocs.xyz/jintonic/gefica/)
[![ReadTheDocs](https://readthedocs.org/projects/gefica/badge)](https://gefica.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

GeFiCa stands for *Ge detector Field Calculator*. It provides classes to calculate static electric fields and potentials in Ge detectors using various coordinates in up to three dimensions. The field potential values together with their grid coordinates are saved in a [ROOT][] [tree][] to take the advantage of the file compression and the [Draw][] function provided by the [ROOT][] [TTree][] class. GeFiCa is provided as a shared library that can be directly loaded by [ROOT][]. All GeFiCa classes can be used directly in a [ROOT][] interactive session. This means that users can modify and run their calculation codes without compilation.

## Directories

Directory | Contents
----------|-----------
src       | source code
docs      | documentation
examples  | [ROOT][] & Python scripts demonstrating usage of GeFiCa

## Get started

1. Make sure that [ROOT][] (version 6 and above) is installed.
2. Execute the following commands in a terminal:

~~~sh
git clone https://github.com/jintonic/gefica.git
cd gefica/src
make
export LD_LIBRARY_PATH=$(PWD):$LD_LIBRARY_PATH
# change LD_LIBRARY_PATH to DYLD_LIBRARY_PATH for MAC
cd ../examples
root -l planar1d.C
~~~

## References

- Overview: https://mediatum.ub.tum.de/node?id=701884
- Field calculation: https://mediatum.ub.tum.de/node?id=969435
- Numerical methods: https://www.mppmu.mpg.de/~jingliu/ECPI/, Lecture 4 and 5

[ROOT]:https://root.cern.ch
[tree]:https://root.cern.ch/root/htmldoc/guides/users-guide/Trees.html
[Draw]:https://root.cern.ch/doc/master/classTTree.html#a73450649dc6e54b5b94516c468523e45
[TTree]:https://root.cern.ch/doc/master/classTTree.html

