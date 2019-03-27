#include "Units.h"
#include "Planar.h"
using namespace GeFiCa;

Planar::Planar(const char *name, const char *title)
   : Detector(name,title), Width(1*cm), Depth(1*cm) {};
