GeFiCa stands for Ge detector Field Calculator.

# Documentation

- http://www.physino.xyz/gefica

## Get started

1. Make sure that [ROOT](https://root.cern.ch) is installed.
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
