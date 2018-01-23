#ifndef SOLENOIDMAP_H
#define SOLENOIDMAP_H

#include "ElementMap.h"
#include <Eigen/LU>
class SolenoidMap:public ElementMap
{
public:
  SolenoidMap();
  SolenoidMap(double length,double field,double gamma=2000);
  void Map();
  void setfield(double field);
  void print();
};

#endif
