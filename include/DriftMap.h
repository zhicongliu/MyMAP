#ifndef DRIFTMAP_H
#define DRIFTMAP_H

#include "ElementMap.h"
class DriftMap:public ElementMap
{
public:
  virtual ~DriftMap();
  DriftMap();
  DriftMap(double length);
  void Map();
  void print();
};
#endif
