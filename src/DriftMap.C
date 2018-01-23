#include "DriftMap.h"
DriftMap::~DriftMap(){}
DriftMap::DriftMap(){}
DriftMap::DriftMap(double length):ElementMap(length)
{
  Map();
}
void DriftMap::Map()
{
//xx
R(0,0)=1;
R(0,1)=length/1000;
R(1,0)=0;
R(1,1)=1;
//yy
R(2,2)=1;
R(2,3)=length/1000;
R(3,2)=0;
R(3,3)=1;
//zz
R(4,4)=1;
R(4,5)=length/gamma/gamma/1000;
R(5,4)=0;
R(5,5)=1;
//zx


//xz

//
}

void DriftMap::print()
{
cout<<"The transfer matrix of drift(length = "<<length<<" mm) is\n"<<R<<endl;
}
