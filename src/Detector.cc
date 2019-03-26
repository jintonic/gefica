#include "Units.h"
#include "Detector.h"
using namespace GeFiCa;

Crystal::Crystal() : Height(1*cm),
   TopImpurity(1e10/cm3), BottomImpurity(1e10/cm3) {};

Detector::Detector() : Crystal() { Bias.push_back(0*volt); }
