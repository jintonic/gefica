#include "Units.h"
#include "Planar.h"
using namespace GeFiCa;

Planar::Planar() : Detector(), TNamed("planar","planar detector"),
   Width(1*cm), Depth(1*cm) {};
