#include "Allelement.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <math.h>
#include <stdlib.h>
int main(int argc,char *argv[])
{
BendingMap a(100,100,10);
a.print();

a.setbeta(0.8);
a.print();
a.setgamma(234);
a.print();

QuadrupoleFocusMap b(10,0.3);
b.print();

QuadrupoleDefocusMap c(10,0.03);
c.print();

DriftMap d(100);
d.print();

ThinlensMap e(123,52.44);
e.twissmap();
e.print();

SolenoidMap f(100,0.8);
f.print();

BeamRotationMap g(30,45,77);
g.print();

GapMap h(1e7,20,100);
h.print();

MatrixXd i=g.getmap();
cout<<"getmap\n"<<i<<endl;

EdgeangleMap j(10,1000,200);
j.print();

MatrixXd k=EdgeangleMap(10,100,200).getmap()*BendingMap(15,100,10).getmap()*EdgeangleMap(10,100,200).getmap();

cout<<"The transfer matrix of a bending magnet with edge effect is\n"<<k<<endl;



return 0;
}
