- Homepage: http://physino.xyz/gefica
- Documentation: https://arxiv.org/abs/2001.02762
- Source code: https://github.com/jintonic/gefica

### Coding convention
The coding convention is similar to that of ROOT. For example,

- Classes and functions all start with capital letters. Word boundaries are indicated by *CamelCase*.
- Classes names are all nouns.
- Function names are all verbs.
- Private member variables all start with letter _f_.
- Boolean variables/functions start with *Is/Are*.
- Indentation is made by three spaces instead of a hard tab to ensure the same appearance of the codes in different editors.

The following exceptions are used to increase the readability of the codes to the user:
- Class names do not have prefix letters, such as _T_ in ROOT. Instead, the name space GeFiCa is used to avoid name collision should GeFiCa be used together with other libraries.
- Configurable member variables are made public to avoid trivial getters and setters. Their first letters are capitalized. Unlike private member variables, they do not have letter _f_ prefixed.

### Class structure
Most of the GeFiCa classes belong to two categories: *grid* and *detector*. Those that are derived from class Grid are used to describe grid setups. Those that are derived from class Detector are used to describe detector configurations. The Grid class inherits a set of arrays from the Points class to describe variables associated with individual grid points, such as coordinates and field values. Names of its derived classes indicate the dimension and coordinates used to construct the grid. For example, X is used for a one dimensional grid in Cartesian coordinates, RhoZ is used for a two dimensional grid in cylindrical coordinates. The Detector class inherits impurity setup from the Crystal class. Its derived classes, such as PointContact, inherit from it the common detector setups, such as bias voltages. A grid class can get boundary conditions and the impurity distribution from a corresponding detector class through a virtual function interface defined in the Grid class:

```cpp
virtual void Grid::SetupWith(Detector&);
```

This is demonstrated in the following code snippet:

```cpp
RhoZ grid; //create 2D Cylindrical grid
PointContact detector; //create detector
//setup grid with detector configuration
grid.SetupWith(detector);
```

The data flow can be the other way around, that is, a detector class gets grid setups from a *grid* class. However, since it is the grid that the SOR process updates instead of the detector configuration, this is a less natural choice. With the current data flow direction, the SOR can then be performed by simply calling

```cpp
  grid.SuccessiveOverRelax();
```

Another choice would be to combine the *detector* and *grid* classes. For example, instead of having both PointContact and RhoZ, we can create a single class called PointContactRhoZ. The advantage of this approaches is that there is no need to pass information from the latter to the former through some interface functions. The drawback is the lack of clarity, the same class object will be used for both detector configuration and grid operation. Considering the main purpose of GeFiCa is to demonstrate the logic, methods, and techniques for field calculation, we chose not to use this approach.

To its root, this is actually a question of to what extend we want to utilize the object-oriented (OO) coding style. Think about two extreme cases. First, we can write everything in a single main function. Second, we can create a class for each individual functionality, such as the impurity profile and the bias voltage. The first approach relies on careful documentation to clarify its internal logic. The second introduces many trivial interfaces to pass information between classes. A balanced approach in between is adopted for GeFiCa.

###  Detector Configurations
Two pieces of information are needed for electric field calculation: first, boundary conditions, and second, the space charge distribution.

Boundary conditions can be set through the detector geometry and voltages on electrodes. Take the previously defined PointContact detector as an example, its basic dimensions can be set as
```cpp
  detector.Radius=3.45*cm;
  detector.Height=5.05*cm;
  detector.PointContactR=1.4*mm;
  detector.PointContactH=0.1*mm;
```

A full list of geometry parameters that can be set for a PointContact detector can be found in the class. Its bias voltages can be set as an array:

```cpp
  detector.Bias[0] = - 2.5*kV; // point-contact voltage
  detector.Bias[1] = 0*volt; // surface contact voltage
```

In case of a segmented detector, the bias voltage array can have more than two elements. The index of an element can be kept the same as the corresponding segment identification number.

It is reasonable to use a first-order polynomial to approximate the space charge distribution in a HPGe crystal. With this simplification, we just need to specify the impurities at the top and the bottom of a crystal given by the manufacturer. For example,

```cpp
  detector.BottomImpurity=3e9/cm3;
  detector.TopImpurity=7e9/cm3;
```

The impurity level at a specific axial position is interpolated in GeFiCa based on these two numbers.
In case of a small crystal, the impurity can be regarded as a constant. Its average impurity can be set as

```cpp
  detector.SetAverageImpurity(3e9/cm3);
```


### Units and Constants
We have seen in previous code snippets that an input parameter in GeFiCa is a product of a number and a unit. Common units and constants for field calculation, together with their conversion rules, are defined in `GeFiCa/src/Units.h`. The following is a snippet of the file:

```cpp
namespace GeFiCa {
  static const double C=1; // Coulomb
  static const double cm=1;
  static const double cm3=cm*cm*cm;
  static const double mm=0.1*cm;
  static const double volt=1;
  static const double kV=1000*volt;
  // vacuum permittivity [C/volt/cm]
  static const double epsilon0 = 8.854187817e-14*C/volt/cm;
  // dielectric constant of Ge
  static const double epsilonR=16;
}
```

The advantage of this unit system is three-fold. First, the code is self-explainable, there is no ambiguity in the unit of an input value. Second, the user has freedom to choose units, such as *mm* instead of *cm* or *kV* instead of *volt* Otherwise, he or she has to use the set of units used for internal calculation. Third, since the unit conversion rules are pre-defined, there is no need to worry about them when using input parameters for internal calculations. The programmer can focus on the logic instead of unit conversion.  This way of handling units is adopted from Geant4.  Most of the units and constants have been defined in Geant4 already. However, since only a small subset of the units are useful for field calculation, they are re-defined in GeFiCa to avoid unnecessary dependence on Geant4.

### Macros and Scripts

A modern C++ interpreter, cling, has been created and adopted as the back-end of the interactive session of ROOT since the version 6 of it. A user can run C++ snippets, sometimes called ROOT macros or scripts, interactively in cling without writing and compiling the *main* function. With immediate feedback after the execution of each line of a script, a user can learn and experiment a new C++ class, a function, or simply a syntax easily. To fulfill its educational purpose, GeFiCa is compiled as a ROOT library. All snippets in previous sections demonstrating the configuration of a detector or the operation of a grid can be run as they are in cling.

ROOT also provides a Python extension module, PyROOT, that allows the user to interact with any ROOT class from the Python interpreter. For users who prefer the Python interpreter to cling, they can call GeFiCa classes with Python syntax directly in the standard Python interpreter.

It is worth noting that cling comes with a Jupyter kernel, which makes it possible to run GeFiCa scripts in a Jupyter notebook with either C++ or python syntax.

All concrete grid and detector classes in GeFiCa inherit the capability to inspect themselves from the TObject class in ROOT. Some standard functions in TObject, such as Dump(), can be used to check the default or user-specified configurations of a grid or detector object. The first column of the output are the member variables of the GeFiCa::X class. The second are their current values. The last are explanations of those variable. These explanations are written as C++ comments after the member variables. They can be parsed by both Doxygen and ROOT to generate code documentation in various formats and contexts.

```cpp
  root [] GeFiCa::X x
  (GeFiCa::X &) Name: x Title: 1D Cartesian coordinate
  root [] x.Dump()
  ==> Dumping object at: 0x00007f76e5d80150, name=x, class=GeFiCa::X
  Src                      ->7f76e5d802a8      -(net impurity concentration)x|Qe|/epsilon
  N1                       101                 number of points along the 1st coordinate
  N2                       0                   number of points along the 2nd coordinate
  N3                       0                   number of points along the 3rd coordinate
  MaxIterations            5000                maximal iterations of SOR to be performed
  RelaxationFactor         1.95                within (0,2), used to speed up convergence
  Tolerance                1e-07               SOR stops when error<Tolerance
  ...
```

Macros are organized in sub-folders in `GeFiCa/examples/` to demonstrate the usage of GeFiCa classes. The `planar/`, `trueCoaxial/`, `hemispherical/`, `pointContact/`, and `segmented/` folders are used to show how to configure specific types of HPGe detectors and then calculate the fields in them. The `analytic/` and the `fenics/` folders contains macros that are independent of the GeFiCa libraries. The macros in the former demonstrate how to calculate and visualize the field distribution in simple HPGe detectors using ROOT. The latter shows Python codes to calculate and visualize the field distribution in a simplified point-contact geometry using FEniCS. All field distributions shown in this work are generated using these macros. A user can learn the topics by at first running these macros to reproduce plots in this work, and then modifying them to meet his/her own needs.

\example analytic/CV.cc
\example analytic/coaxial.cc
\example analytic/getVd.cc
\example analytic/planar.cc
\example analytic/spherical.cc

\example fenics/minippc.py

\example hemispherical/compare2analytic.cc

\example planar/calculateC.cc
\example planar/checkInitialization.cc
\example planar/compare2CG.cc
\example planar/compare2analytic.cc
\example planar/optimizeRelaxationFactor.cc
\example planar/showConvergingSteps.cc

\example pointContact/calculateFields.cc
\example pointContact/checkInitialization.cc
\example pointContact/comapare2RhoZ.cc
\example pointContact/comapare2analytic.cc
\example pointContact/comapare2fenics.cc
\example pointContact/comapare2fieldgen.cc
\example pointContact/drawFields.cc
\example pointContact/drawFields.py
\example pointContact/g2f.cc
\example pointContact/optimizeRelaxationFactor.cc
\example pointContact/search4Vd.cc

\example segmented/drawSliceInPhi.cc
\example segmented/drawSliceInZ.cc

\example trueCoaxial/checkInitialization.cc
\example trueCoaxial/compare2analytic.cc
\example trueCoaxial/verifyCV.cc
