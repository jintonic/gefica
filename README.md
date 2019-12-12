[![Doxygen](https://codedocs.xyz/jintonic/gefica.svg)](https://codedocs.xyz/jintonic/gefica/)
[![ReadTheDocs](https://readthedocs.org/projects/gefica/badge)](https://gefica.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Code size](https://img.shields.io/github/languages/code-size/jintonic/gefica.svg?style=flat)
![](https://img.shields.io/github/languages/top/jintonic/gefica.svg?style=flat)

## Introduction

GeFiCa stands for *Ge detector Field Calculator*. It provides classes to calculate static electric fields and potentials in Ge detectors using various coordinates in up to three dimensions. The field potential values together with their grid coordinates are saved in a [ROOT][] [tree][] to take the advantage of the file compression and the [Draw][] function provided by the [ROOT][] [TTree][] class. GeFiCa is provided as a shared library that can be directly loaded by [ROOT][]. All GeFiCa classes can be used directly in a [ROOT][] interactive session or a Jupyter notebook. Users can modify and run their calculation codes without compilation.

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
cd ../examples/segmented
root -l drawSliceInPhi.cc
~~~

## References

- Overview: <https://mediatum.ub.tum.de/node?id=701884>
- Field calculation: <https://mediatum.ub.tum.de/node?id=969435>
- Numerical methods: <https://www.mppmu.mpg.de/~jingliu/ECPI/>, Lecture 4 and 5

## Similar projects

- <https://github.com/JuliaPhysics/SolidStateDetectors.jl>
- <https://radware.phy.ornl.gov/MJ/mjd_siggen/>

## Git submodule of MJD fieldgen

GeFiCa contains a Git submodule of [MJD fieldgen](https://github.com/jintonic/siggen). It is a mirror of the MJD siggen (including fieldgen) subversion repository in GitHub. The mirroring is done by

```sh
apt install git-svn
git svn clone https://radware.phy.ornl.gov/MJ/mjd_siggen/ siggen
git remote add origin git@github.com:jintonic/siggen.git
git push origin master
```

Adding it as a git submodule to GeFiCa is done by

```sh
cd /path/to/GeFiCa
git submodle add https://github.com/jintonic/siggen.git examples/pointContact/fieldgen
git commit -m "added siggen submodule"
git push
```

Cloning GeFiCa with fieldgen can be done by doing

```sh
git clone https://github.com/jintonic/gefica.git
cd gefica
git submodule init
git submodule update
```

Updating GeFiCa and fieldgen can be done by doing

```sh
cd gefica
git pull
cd examples/pointContact/fieldgen
git pull
```

[ROOT]:https://root.cern.ch
[tree]:https://root.cern.ch/root/htmldoc/guides/users-guide/Trees.html
[Draw]:https://root.cern.ch/doc/master/classTTree.html#a73450649dc6e54b5b94516c468523e45
[TTree]:https://root.cern.ch/doc/master/classTTree.html

