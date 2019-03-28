#include "Units.h"
#include "TrueCoaxial.h"
using namespace GeFiCa;

TrueCoaxial::TrueCoaxial(const char *name, const char *title)
   : Detector(name, title), Radius(3*cm), BoreR(0.5*cm)
{ Bias.push_back(1*kV); }
