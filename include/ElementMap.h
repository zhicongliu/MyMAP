#ifndef ELEMENTMAP_H
#define ELEMENTMAP_H

#define HUGE			1
#define _MASS		 	(1.672621777e-27*10000)	//kg
#define _Q			(1.602176565e-19*10000)		//C
#define C_light			299792458		//m/s
#define Dielectric_const	8.854187817e-12		//F/m

#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace Eigen;
using namespace std;
class ElementMap
{
public:
  virtual ~ElementMap();
  ElementMap();
  ElementMap(double);
  ElementMap(double,double);
//  ElementMap(const ElementMap &C);
  virtual void Map();
  virtual MatrixXd getmap();
  void setlength(double);
  void setbeta(double);
  void setgamma(double);
  virtual void print();
  int sign(double x)
{
  if(x < 0) return -1;
  else if(x==0) return 0;
  else return 1;
}
protected:
  MatrixXd R=MatrixXd::Identity(6,6);
  double beta,gamma;
  double gradient,field,k,kx,ky,h,Bro,aperture;
  double length,location;
};
#endif
