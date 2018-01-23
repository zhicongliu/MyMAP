#ifndef MYFUNC_H
#define MYFUNC_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h> 
#include <math.h>
#include <Eigen/Dense>//!!! how is this one related to string.remove???

using namespace std;
using namespace Eigen;
void StringSplit(string s,char splitchar,vector<string> &vec);
int getlinenumber(char *p);
double gaussrand();
double bessi0(double x);
double bessi1(double x);
MatrixXd twissmapx(const MatrixXd &R);//beta alpha gamma,in x direction test
MatrixXd twissmapy(const MatrixXd &R);//beta alpha gamma,in y direction test
MatrixXd twissmapz(const MatrixXd &R);//beta alpha gamma,in z direction test

#endif
