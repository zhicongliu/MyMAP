#ifndef QUADRUPOLEDEFOCUS_H
#define QUADRUPOLEDEFOCUS_H

#include "ElementMap.h"

class QuadrupoleDefocusMap:public ElementMap
{
public:
  virtual ~QuadrupoleDefocusMap();
  QuadrupoleDefocusMap();
  QuadrupoleDefocusMap(double,double,double gamma=2000);//(mm,T/m)
  void setgradient(double);
  virtual void Map();
  virtual void print();
};

#endif
