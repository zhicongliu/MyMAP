#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h> 
#include <math.h>
#include <Eigen/Dense>
#include "MyFunc.h"

using namespace std;
using namespace Eigen;

void KVdistributionA(double emitx,double betax,double emity,double betay,double alpx,double alpy,int particlenumber=10000,char *p="particles.dat",double sig_z=1,double sig_pz=0.1);

void KVdistributionS(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,int particlenumber=10000,char *p="particles.dat",double sig_z=1,double sig_pz=0.1);

void KVdistributionT(const vector<double> &twiss,double emitx,double emity,int particlenumber=10000,char *p="particles.dat",double sig_z=1,double sig_pz=0.1);

void KVdistribution6A(double emitx,double betax,double emity,double betay,double emitz,double betaz,double alpx,double alpy,double alpz,int particlenumber=10000,char *p="particles.dat");

void KVdistribution6S(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,double sig_z,double sig_pz,double sig_zpz,int particlenumber=10000,char *p="particles.dat");

void KVdistribution6T(const vector<double> &twiss,double emitx,double emity,double emittz,int particlenumber=10000,char *p="particles.dat");

void distribution(double a,double b,double c,double d,double rx=0,double ry=0,double rxy=0,double rxtyt=0,double rxyt=0,double rxty=0,int particlenumber=10000,char *p="particles.dat");

void WBdistribution(double x1,double x2,double alpx,double y1,double y2,double alpy,int particlenumber=10000);

#endif
