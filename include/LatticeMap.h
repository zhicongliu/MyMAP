#ifndef LATTICEMAP_H
#define LATTICEMAP_H

#include "AllElement.h"
#include "MyFunc.h"

#include <vector>
#include <map>
#include <string>

#include <fstream>
#include <iostream>
#include <stdio.h> 

#include <mpich/mpi.h>

//#include <boost/any.hpp>
using namespace std;
using namespace Eigen;

struct Element
{
  bool operator== (const Element& e)
  {
    if(name==e.name&&p.size()==e.p.size())
    {
      for(int i=0;i<p.size();i++)
      {
        if(p[i]!=e.p[i])
          return false;
      }
      return true;
    }
    return false;
  }
  bool operator!= (const Element& e)
  {
    if(name==e.name&&p.size()==e.p.size())
    {
      for(int i=0;i<p.size();i++)
      {
        if(p[i]!=e.p[i])
          return true;
      }
      return false;
    }
    return true;
  }
  string name;
  string path;
  vector<double> p;
};
class LatticeMap
{
public:
  virtual ~LatticeMap();
  LatticeMap();
  LatticeMap(const LatticeMap& );
  LatticeMap& operator= (const LatticeMap&);
  bool operator== (const LatticeMap&);
  void callelement(const Element &v,ifstream &infile);
  void add(const int &num,const Element &v);
  void modify(const int &num,const Element &v);
  void read(char *p,double EGamma,int numberdivided=1,int elemlimit=0);

  void del(int Nelement);
  void deleteall();
  void print();
  void prelattice(const char *in,const char *out,int n);
  void prelattice1(const char *in,const char *out);
  void addEGamma(const double num);

  MatrixXd getMat(int L,double &x,double &dx,double &y,double &dy,double &z,double &dp);
  vector<MatrixXd>		Mat;
  vector<Element>		ELE;

  vector<double> EBeta;
  vector<double> EGamma;
  double fitness;
protected:
  //Element
  Element elemtemp;

  vector<QuadrupoleMap>		FQ;
  vector<BendingMap>		BM;
  vector<DriftMap>		DL;
  vector<SolenoidMap>		SN;
  vector<RFQMap>		RFQ;
  double rfqtypepre,rfqtypenow,rfqtypenext;
  
  //Common
  int NDivide;

  int Mfail=0,particle_count=0,totalparticalnumber;
  int i=0,j,k;
  
  //MPI
  int myid,numprocs;
  
  //RFQ
  int positionpre,positionnow,positionnext;

};

#endif
