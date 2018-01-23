#include "QuadrupoleMap.h"

QuadrupoleMap::~QuadrupoleMap(){}

QuadrupoleMap::QuadrupoleMap(){}

QuadrupoleMap::QuadrupoleMap(double length,double gradient,double gamma)//(mm,T/m)
{
  setlength(length);
  setgradient(gradient);
  setgamma(gamma);
  Map();
}
void QuadrupoleMap::Map()
{
  Bro=_MASS*C_light/_Q*gamma*beta*1000;//T mm
  //cout<<"BRO  "<<Bro<<"   "<<Bro2<<endl;
  k=sqrt(abs(gradient/Bro));
  //cout<<k*k*1e6<<endl;
  if(gradient>0)
  {
//xx
	R(0,0)=cos(k*length);
	R(0,1)=sin(k*length)/k/1000;
	R(1,0)=-k*sin(k*length)*1000;
	R(1,1)=cos(k*length);
//yy
	R(2,2)=cosh(k*length);
	R(2,3)=sinh(k*length)/k/1000;
	R(3,2)=k*sinh(k*length)*1000;
	R(3,3)=cosh(k*length);
//zz
	R(4,4)=1;
	R(4,5)=length/gamma/gamma/1000;
	R(5,4)=0;
	R(5,5)=1;
}
else if(gradient<=0)
{
//xx
	R(0,0)=cosh(k*length);
	R(0,1)=sinh(k*length)/k/1000;
	R(1,0)=k*sinh(k*length)*1000;
	R(1,1)=cosh(k*length);
//yy
	R(2,2)=cos(k*length);
	R(2,3)=sin(k*length)/k/1000;
	R(3,2)=-k*sin(k*length)*1000;
	R(3,3)=cos(k*length);
//zz
	R(4,4)=1;
	R(4,5)=length/gamma/gamma/1000;
	R(5,4)=0;
	R(5,5)=1;
}
//
}

void QuadrupoleMap::setgradient(double gradient)
{
this->gradient=gradient/1000;//T/m transfer to T/mm
}


void QuadrupoleMap::print(){
cout<<"The transfer matrix of quadrupole(length = "<<length<<" mm,gradient = "<<gradient*1000<<" T/m) is\n"<<R<<endl;
}




