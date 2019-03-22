#include "RTheta.h"
#include "Units.h"
using namespace GeFiCa;

RTheta::RTheta(size_t n_r, size_t n_theta) : R(n_r)
{ N2=n_theta; SetName("rt"); SetTitle("2D spherical coordinates"); }
