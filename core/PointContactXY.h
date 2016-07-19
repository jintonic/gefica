#ifndef GEFICA_POINTCONTACTXY_H
#define GEFICA_POINTCONTACTXY_H

#include "XY.h"

namespace GEFICA { class PointContactXY; }

class GEFICA::PointContactXY : public GEFICA::XY
{
   public :
      PointContactXY(int ix,int iy) : XY(ix,iy) {};
      void SetVoltage(double voltage,double topbegin,double topend); 
      ClassDef(PointContactXY,1);
};

#endif
