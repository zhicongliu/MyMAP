#include "BeamMap.h"
#include <fstream>
#include <iostream>
BeamMap::~BeamMap()
{
clear();
}
BeamMap::BeamMap(){}
BeamMap::BeamMap(char particletype,double gamma,int number)
{
  initial(particletype,gamma,number);
}

void BeamMap::initial(char particletype,double gamma,int number)
{
  if(particletype=='e')
  {
    beam_staticenergy=_MASS*C_light*C_light;
  }
  else if(particletype=='p')
  {
    beam_staticenergy=_MASS*C_light*C_light;
  }
  else
  {
    cout<<"error! particletype not defined"<<endl;
  }
  beam_energy=beam_staticenergy*gamma;
  setgamma_beam_energy(gamma);
  particleLIVEnumber=number;
  clear();
  particle.resize(particleLIVEnumber);
}
void BeamMap::clear()
{
particle.clear();
}

void BeamMap::settwissx(double alpha,double beta)
{
if(beta>=0)
{
  TAlphax=alpha;
  TBetax=beta;
  TGammax = (1+TAlphax*TAlphax)/beta;
}
else{cout<<"Error! betax= "<<beta<<endl;}
}
void BeamMap::settwissy(double alpha,double beta)
{
if(beta>=0)
{
  TAlphay=alpha;
  TBetay=beta;
  TGammay = (1+TAlphay*TAlphay)/beta;
}
else{cout<<"Error! betay= "<<beta<<endl;}
}

void BeamMap::setbeta_beam_energy(double a)
{
if(a>0&&a<1)
{
  beta_beam_energy=a;
  gamma_beam_energy=1/sqrt(1-beta_beam_energy*beta_beam_energy);
  beam_energy = gamma_beam_energy * beam_staticenergy;
}
else if(beta_beam_energy=1)
{
  beta_beam_energy=a;}
else
{
  cout<<"beta_energy wrong!"<<endl;
}
}

void BeamMap::setgamma_beam_energy(double a)
{
if(a>=1)
{
  gamma_beam_energy=a;
  beta_beam_energy=sqrt(1-1/(gamma_beam_energy*gamma_beam_energy));
  beam_energy= gamma_beam_energy* beam_staticenergy;
}
else
{
  cout<<"gamma_energy wrong!"<<endl;
}
}

void BeamMap::caculate_emittance()
{
double x_total=0,dx_total=0;
x_sigma =0;
dx_sigma=0;
xdx_sigma=0;
double y_total=0,dy_total=0;
y_sigma=0;
dy_sigma=0;
ydy_sigma=0;
double z_total=0,dz_total=0;
z_sigma=0;
dz_sigma=0;
zdz_sigma=0;
double x_average,y_average,z_average,dx_average,dy_average,dz_average;
for(int i=0;i<particleLIVEnumber;i++)
{
	x_total=x_total+particle[i].getlocationx();
	dx_total=x_total+particle[i].getdirectionx();
	y_total=y_total+particle[i].getlocationy();
	dy_total=y_total+particle[i].getdirectiony();
	z_total=z_total+particle[i].getlocationz();
	dz_total=dz_total+particle[i].getdirectionz();
}
x_average=x_total/particleLIVEnumber;
dx_average=dx_total/particleLIVEnumber;
y_average=y_total/particleLIVEnumber;
dy_average=dy_total/particleLIVEnumber;
z_average=z_total/particleLIVEnumber;
dz_average=dz_total/particleLIVEnumber;
for(int i=0;i<particleLIVEnumber;i++)
{
	x_sigma=x_sigma+(particle[i].getlocationx()-x_average)*(particle[i].getlocationx()-x_average);
	dx_sigma=dx_sigma+(particle[i].getdirectionx()-dx_average)*(particle[i].getdirectionx()-dx_average);
	xdx_sigma=xdx_sigma+(particle[i].getlocationx()-x_average)*(particle[i].getdirectionx()-dx_average);
	y_sigma=y_sigma+(particle[i].getlocationy()-y_average)*(particle[i].getlocationy()-y_average);
	dy_sigma=dy_sigma+(particle[i].getdirectiony()-dy_average)*(particle[i].getdirectiony()-dy_average);
	ydy_sigma=ydy_sigma+(particle[i].getlocationy()-y_average)*(particle[i].getdirectiony()-dy_average);
	z_sigma=z_sigma+(particle[i].getlocationz()-z_average)*(particle[i].getlocationz()-z_average);
	dz_sigma=dz_sigma+(particle[i].getdirectionz()-dz_average)*(particle[i].getdirectionz()-dz_average);
	zdz_sigma=zdz_sigma+(particle[i].getlocationz()-z_average)*(particle[i].getdirectionz()-dz_average);
}
x_sigma=x_sigma/particleLIVEnumber;
dx_sigma=dx_sigma/particleLIVEnumber;
xdx_sigma=xdx_sigma/particleLIVEnumber;
	correlationx=xdx_sigma/(sqrt(x_sigma*dx_sigma));
y_sigma=y_sigma/particleLIVEnumber;
dy_sigma=dy_sigma/particleLIVEnumber;
ydy_sigma=ydy_sigma/particleLIVEnumber;
	correlationy=ydy_sigma/(sqrt(y_sigma*dy_sigma));

z_sigma=z_sigma/particleLIVEnumber;
dz_sigma=dz_sigma/particleLIVEnumber;
zdz_sigma=zdz_sigma/particleLIVEnumber;
emittancex=sqrt(x_sigma*dx_sigma-xdx_sigma*xdx_sigma);
emittancey=sqrt(y_sigma*dy_sigma-ydy_sigma*ydy_sigma);
emittancez=sqrt(z_sigma*dz_sigma-zdz_sigma*zdz_sigma);

TBetax=x_sigma/emittancex;
TBetay=y_sigma/emittancey;
TBetaz=z_sigma/emittancez;
TGammax=dx_sigma/emittancex;
TGammay=dy_sigma/emittancey;
TGammaz=dz_sigma/emittancez;
if(TBetax*TGammax>1)
  TAlphax=sqrt(TBetax*TGammax-1);
else
  TAlphax=0;
if(TBetay*TGammay>1)
  TAlphay=sqrt(TBetay*TGammay-1);
else
  TAlphay=0;
if(TBetaz*TGammaz>1)
  TAlphaz=sqrt(TBetaz*TGammaz-1);
else
  TAlphaz=0;

}
void BeamMap::output(char *p)
{
ofstream out(p);
for(int i=0;i<particleLIVEnumber;i++)
{
	out<<particle[i].getlocationx()<<"	 ";
	out<<particle[i].getdirectionx()<<"	 ";
	out<<particle[i].getlocationy()<<"	 ";
	out<<particle[i].getdirectiony()<<"	 ";
	out<<particle[i].getlocationz()<<"	 ";
	out<<particle[i].getdirectionz()<<endl;
}
out.close();
}

//=========================================
double BeamMap::getsigmax()
{
return x_sigma;
}
double BeamMap::getsigmay()
{
return y_sigma;
}
double BeamMap::getsigmaz()
{
return z_sigma;
}
double BeamMap::getsigmadx()
{
return dx_sigma;
}
double BeamMap::getsigmady()
{
return dy_sigma;
}
double BeamMap::getsigmadz()
{
return dz_sigma;
}
double BeamMap::getsigmaxdx()
{
return xdx_sigma;
}
double BeamMap::getsigmaydy()
{
return ydy_sigma;
}
double BeamMap::getsigmazdz()
{
return zdz_sigma;
}


double BeamMap::getTBetax()
{
return TBetax;
}
double BeamMap::getTAlphax()
{
return TAlphax;
}
double BeamMap::getTGammax()
{
return TGammax;
}
double BeamMap::getTBetay()
{
return TBetay;
}
double BeamMap::getTAlphay()
{
return TAlphay;
}
double BeamMap::getTGammay()
{
return TGammay;
}

double BeamMap::getTBetaz()
{
return TBetaz;
}
double BeamMap::getTAlphaz()
{
return TAlphaz;
}
double BeamMap::getTGammaz()
{
return TGammaz;
}

double BeamMap::getbeta_beam_energy()
{
return beta_beam_energy;
}
double BeamMap::getgamma_beam_energy()
{
return gamma_beam_energy;
}
double BeamMap::getemittancex()
{
return emittancex;
}
double BeamMap::getemittancey()
{
return emittancey;
}
double BeamMap::getcorrelationx()
{
return correlationx;
}
double BeamMap::getcorrelationy()
{
return correlationy;
}
double BeamMap::getparticlenumber()
{
return particleLIVEnumber;
}
void BeamMap::setparticlenumber(int a)
{
  particleLIVEnumber=a;
}

ParticleMap BeamMap::getparticle(int i)
{
return particle[i];
}

void BeamMap::setparticle(int i,double a,double b,double c,double d,double e,double f)
{
particle[i].setparticle(a,b,c,d,e,f);
}
void BeamMap::setlost(int i,int COUNT)
{
swap(particle[i],particle[particleLIVEnumber-1]);
particleLIVEnumber--;
particle[particleLIVEnumber-1].lost=COUNT;
}



