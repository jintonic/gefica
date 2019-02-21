/** @file Units.h
 * A file defining commonly used units & constants.
 */
#ifndef GeFiCa_UNITS_H
#define GeFiCa_UNITS_H

namespace GeFiCa { 
   static const double C=1;
   static const double cm=1;
   static const double cm3=cm*cm*cm;
   static const double volt=1;
   static const double pF=C/volt*1e-12;
   static const double Qe=1.6e-19*C;  ///< electron charge in Coulomb [C]
   static const double epsilon=16*8.854187817e-14*C/volt/cm; ///< dielectric constant of Ge [C/volt/cm]
}

#endif

