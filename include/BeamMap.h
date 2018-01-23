#ifndef BEAM_MAP
#define BEAM_MAP

#include "ParticleMap.h"
#include "iostream"
#include "ElementMap.h"
#include <math.h>
#include <vector>

using namespace std;

class BeamMap
{
public:
  virtual ~BeamMap();
  BeamMap(char,double,int);
  BeamMap();
  void initial(char,double,int);
  void settwissx(double alpha,double beta);
  void settwissy(double alpha,double beta);
  void setbeta_beam_energy(double a);
  void setgamma_beam_energy(double a);
  void caculate_emittance();
  void output(char*);
  void clear();
  
  double getsigmax();
  double getsigmay();
  double getsigmaz();
  double getsigmadx();
  double getsigmady();
  double getsigmadz();
  double getsigmaxdx();
  double getsigmaydy();
  double getsigmazdz();
  
  double getTBetax();
  double getTAlphax();
  double getTGammax();
  double getTBetay();
  double getTAlphay();
  double getTGammay();
  double getTBetaz();
  double getTAlphaz();
  double getTGammaz();
  
  double getbeta_beam_energy();
  double getgamma_beam_energy();
  double getemittancex();
  double getemittancey();
  double getcorrelationx();
  double getcorrelationy();
  double getparticlenumber();
  void setparticlenumber(int a);
  ParticleMap getparticle(int i);
  void setparticle(int i,double a,double b,double c,double d,double e,double f);
  void setlost(int,int);
  
  vector<ParticleMap> particle;
  int particleLIVEnumber=0;
  
protected:
double beam_staticenergy,beam_energy;
double beta_beam_energy;
double gamma_beam_energy;
double TBetax,TBetay,TBetaz,TAlphax,TAlphay,TAlphaz,TGammax,TGammay,TGammaz;
double x_sigma,y_sigma,z_sigma,dx_sigma,dy_sigma,dz_sigma,xdx_sigma,ydy_sigma,zdz_sigma;
double emittancex,emittancey,emittancez,correlationx,correlationy,correlationz;


};


#endif
