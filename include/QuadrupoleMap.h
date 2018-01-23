#ifndef QUADRUPOLEFOCUSMAP_H
#define QUADRUPOLEFOCUSMAP_H

#include "ElementMap.h"
class QuadrupoleMap:public ElementMap
{
public:
  virtual ~QuadrupoleMap();
  QuadrupoleMap();
  QuadrupoleMap(double length,double gradient,double gamma=2000);//(mm,T/m)
  virtual void Map();
  void setgradient(double gradient);
  virtual void print();
};

#endif
