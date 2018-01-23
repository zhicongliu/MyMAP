
#include "QuadrupoleDefocusMap.h"

QuadrupoleDefocusMap::~QuadrupoleDefocusMap(){}

QuadrupoleDefocusMap::QuadrupoleDefocusMap(){}

QuadrupoleDefocusMap::QuadrupoleDefocusMap(double length,double gradient,double gamma):ElementMap(length)
{
	this->gradient=gradient/1000;
	this->gamma=gamma;
	Map();
}

void QuadrupoleDefocusMap::setgradient(double gradient)
{
	this->gradient=gradient/1000;//   T/m to T/mm
}


void QuadrupoleDefocusMap::Map()
{
	Bro=0.511*gamma/1000/0.29979*1000;
	k=sqrt(abs(gradient/Bro));
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
};



void QuadrupoleDefocusMap::print(){
	cout<<"The transfer matrix of defocusing quadrupole(length = "<<length<<" mm,gradient = "<<gradient*1000<<" T/m) is\n"<<R<<endl;
}




