#include "SolenoidMap.h"

SolenoidMap::SolenoidMap(){}
SolenoidMap::SolenoidMap(double length,double field,double gamma)
{
  setlength(length);
  setfield(field);
  setgamma(gamma);
  Map();
}
void SolenoidMap::Map()
{
  Bro=_MASS/_Q*gamma*beta*C_light*1000;//T mm     //
  k=field/(2*Bro);
  cout<<"bro & k:   "<<Bro<<"    "<<k<<endl;
//xx
R(0,0)=cos(k*length)*cos(k*length);
R(0,1)=sin(k*length)*cos(k*length)/k/1000;
R(1,0)=-k*sin(k*length)*cos(k*length)*1000;
R(1,1)=cos(k*length)*cos(k*length);
//yy
R(2,2)=cos(k*length)*cos(k*length);
R(2,3)=sin(k*length)*cos(k*length)/k/1000;
R(3,2)=-k*sin(k*length)*cos(k*length)*1000;
R(3,3)=cos(k*length)*cos(k*length);
//zz
R(4,4)=1;
R(4,5)=length/gamma/gamma/1000;
R(5,4)=0;
R(5,5)=1;
//xy
R(0,2)=sin(k*length)*cos(k*length);
R(0,3)=sin(k*length)*sin(k*length)/k/1000;
R(1,2)=-k*sin(k*length)*sin(k*length)*1000;
R(1,3)=sin(k*length)*cos(k*length);
//yx
R(2,0)=-sin(k*length)*cos(k*length);
R(2,1)=-sin(k*length)*sin(k*length)/k/1000;
R(3,0)=k*sin(k*length)*sin(k*length)*1000;
R(3,1)=-sin(k*length)*cos(k*length);
//zx
R(4,0)=0;
R(4,1)=0;
R(5,0)=0;
R(5,1)=0;
//xz
R(0,4)=0;
R(0,5)=0;
R(1,4)=0;
R(1,5)=0;
}

void SolenoidMap::setfield(double field)
{
this->field = field;
}

void SolenoidMap::print()
{
cout<<"The transfer matrix of Solenoid(length = "<<length<<" mm,Field = "<<field<<" T) is\n"<<R<<endl;
cout<<R.block<4,4>(0,0).determinant()<<endl;
cout<<"R.block<2,2>(0,0) "<<R.block<2,2>(0,0).determinant()<<endl;
cout<<"R.block<2,2>(0,2) "<<R.block<2,2>(0,2).determinant()<<endl;
cout<<"R.block<2,2>(2,0) "<<R.block<2,2>(2,0).determinant()<<endl;
cout<<"R.block<2,2>(2,2) "<<R.block<2,2>(2,2).determinant()<<endl;
}
