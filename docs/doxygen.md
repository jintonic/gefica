- Homepage: http://physino.xyz/gefica
- Documentation: https://gefica.readthedocs.io
- Github: https://github.com/jintonic/gefica

# Coding convention
We follow ROOT coding convertion with the following exceptions:

- Classes do not start with `T`, instead, they are all in a namespace called `GeFiCa`
- Public member variables do not start with `f`

## Variables in TTree
short variable names in TTree.

# Class structure
Classes named as a type of Ge detectors, such as Planar1D, TrueCoaxial2D, Sphere, etc., are used to setup boundary conditions for the field calculation. Classes named as a type of coordinates, such as XYZ, RThetaPhi, RhoPhi, etc., are used to update field grid points based on boundary conditions and Poisson's equation. Calculations for different dimensions are split to different classes because of performance considerations - it is slow to use a class designed for 3D calculation to calculate 1D fields. Classes dealing with higher dimensions inherit from lower dimension classes to reuse calculation codes. Voltages applied to a detector are saved as public members of class X, which are inherited by all other classes. Constant impurity concentration is saved the same way. Classes named as a type of coordinates are hence not limited to coordinates-related things only, they take care of Ge crystal properties as well. While classes named as a type of Ge detectors takes care of dimensions of electrodes and initialize the boundary conditions accordingly.

# ROOT macros and Python scripts as examples
No need to compile before run.

# Limitation
- clover cannot be handled. We need something as general as FineCS
