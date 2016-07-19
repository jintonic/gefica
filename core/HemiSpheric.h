#ifndef IDNTNSHAPE_H
#define IDNTNSHAPE_H

#include "halfball.h"

namespace GEFICA { class idntnshape; }

class GEFICA::idntnshape : public GEFICA::halfball
{
  public :
    idntnshape(int ix,int iy,int iz) : halfball(ix,iy,iz) {};
    void SetVoltage(double voltage,double r); 
    ClassDef(idntnshape,1);
};

#endif
