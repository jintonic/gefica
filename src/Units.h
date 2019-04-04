/** @file Units.h
 * A file defining commonly used units & constants.
 */
#ifndef GeFiCa_Units
#define GeFiCa_Units
namespace GeFiCa { 
   static const double ns=1; ///< nano second
   static const double us=1e3*ns; ///< micro second
   static const double ms=1e6*ns; ///< mini second
   static const double sec=1e9*ns; ///< second
   static const double C=1; ///< charge unit Coulomb
   static const double cm=1; ///< centimeter
   static const double cm2=cm*cm; ///< centimeter squared
   static const double cm3=cm*cm*cm; ///< centimeter cubed
   static const double mm=0.1*cm; ///< minimeter
   static const double mm2=mm*mm; ///< minimeter squared
   static const double mm3=mm*mm*mm; ///< minimeter cubed
   static const double nm=1e-6*mm; ///< nanometer
   static const double um=1e-3*mm; ///< micrometer
   static const double volt=1;
   static const double kV=1000*volt; ///< kilo volt
   static const double pF=C/volt*1e-12; ///< pico farad
   static const double nF=C/volt*1e-9; ///< nano farad
   static const double Qe=1.6e-19*C;  ///< electron charge in Coulomb [C]
   static const double epsilon0=8.854187817e-14*C/volt/cm; ///< vacumm permittivity [C/volt/cm]
   static const double epsilonR=16; ///< dielectric constant of Ge
   static const double epsilon=epsilonR*epsilon0; ///< permittivity of Ge [C/volt/cm]

   static const double Pi=3.14159265;
}
#endif

