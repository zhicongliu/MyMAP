#include "Distribution.h"

void KVdistributionT(const vector<double> &twiss,double emitx,double emity,int particlenumber,char *p,double sig_z,double sig_pz)
{
  if(twiss.size()<4)
    return;
  double sigx	= (twiss[0]*emitx)/4;
  double sigpx	= ((twiss[1]*twiss[1]+1)/twiss[0]*emitx)/4;
  double sigxpx	= sqrt(sigx*sigpx-emitx*emitx/16);
  if(twiss[1]==0)
    sigxpx=0;
  double sigy	= (twiss[2]*emity)/4;
  double sigpy	= ((twiss[3]*twiss[3]+1)/twiss[2]*emity)/4;
  double sigypy	= sqrt(sigy*sigpy-emity*emity/16);
  if(twiss[3]==0)
    sigypy=0;
  //  cout<<"sig : "<<sigx<<"  "<<sigpx<<"  "<<sigxpx<<endl;
  KVdistributionS(sigx,sigpx,sigxpx,sigy,sigpy,sigypy,particlenumber,p);
}


void KVdistributionS(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,int particlenumber,char *p,double sig_z,double sig_pz)
{
  if(sig_x*sig_px<sig_xpx*sig_xpx)
  {
    cout<<"error sigma_x_px"<<endl;
    return;
  }
  MatrixXd CM = MatrixXd::Zero(2,2);
  //  MatrixXd DM = MatrixXd::Zero(3,3);
  VectorXd eivals;
  double x1=2*sig_xpx/sqrt(sig_px);
  double px1=2*sqrt(sig_px);
  double x2=2*sqrt(sig_x);
  double px2=2*sig_xpx/sqrt(sig_x);
  //  double Ax=px2/x1/(x2*px1-x1*px2);
  double Bx=-2*px2/px1/(x2*px1-x1*px2);
  double Cx=x2/px1/(x2*px1-x1*px2);
  double Ax=(1+px2*px2*Cx)/x2/x2;
  CM<<Ax,   Bx/2,
    Bx/2, Cx;
  /*
     DM<<Ax,   Bx/2, 0,
     Bx/2, Cx  , 0,
     0  ,  0  ,  1;
   */
  eivals = CM.selfadjointView<Lower>().eigenvalues();
  if(sig_x<sig_px)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ax=sqrt(abs(1/eivals(0)));
  double bx=sqrt(abs(1/eivals(1)));
  double phix;
  if(Ax!=Cx)
    phix=1.0/2*atan(Bx/(Ax-Cx))*180/M_PI;
  else
    phix=0;
  double emix=ax*bx;
  double y1=2*sig_ypy/sqrt(sig_py);
  double py1=2*sqrt(sig_py);
  double y2=2*sqrt(sig_y);
  double py2=2*sig_ypy/sqrt(sig_y);
  double By=-2*py2/py1/(y2*py1-y1*py2);
  double Cy=y2/py1/(y2*py1-y1*py2);
  double Ay=(1+py2*py2*Cy)/y2/y2;
  CM<<Ay,   By/2,
    By/2, Cy;
  eivals = CM.selfadjointView<Lower>().eigenvalues();
  if(sig_y<sig_py)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ay=sqrt(1/eivals(0));
  double by=sqrt(1/eivals(1));
  double phiy;
  if(Ay!=Cy)
    phiy=1.0/2*atan(By/(Ay-Cy))*180/M_PI;
  else
    phiy=0;
  double emiy=ay*by;
  /*
     cout<<"x1 = "<<x1<<" , px1= "<<px1<<endl;
     cout<<"x2 = "<<x2<<" , px2= "<<px2<<endl;
     cout<<Ax<<"  "<<Bx<<"  "<<Cx<<"  ABC"<<endl;
     cout<<phix<<"  phi  "<<phiy<<endl;
     cout<<emix<<"  "<<ax<<"  "<<emiy<<"  "<<ay<<"  "<<phix<<"  "<<phiy<<"  "<<particlenumber<<"  "<<p<<endl;;
     cout<<"emit x exp : "<<4*sqrt(sig_x*sig_px-sig_xpx*sig_xpx)<<endl;
   */
  KVdistributionA(emix,ax,emiy,ay,phix,phiy,particlenumber,p);
  /*
     particleread(p);
     cout<<beam.getsigmax()<<endl;
     cout<<beam.getsigmadx()<<endl;
     cout<<beam.getsigmaxdx()<<endl;
   */
}

void KVdistributionA(double emitx,double betax,double emity,double betay,double alpx,double alpy,int particlenumber,char *p,double sig_z,double sig_pz)
{
  ofstream par(p);
  double ax,axs,ay,ays,emitratio,F;
  double x,xs,y,ys,z,zs;
  double anglebetax,anglebetay,zetax,zetay;
  double alphax=alpx/180*M_PI;
  double alphay=alpy/180*M_PI;
  double x1=betax;
  double x2=emitx/betax;
  double y1=betay;
  double y2=emity/betay;
  double limity_x,limity_y,limitx_x,limitx_y;
  if(alpx!=0&&alpx!=90)
  {
    //limity_x=sqrt(1.0/pow(((2.0*cos(alphax)*cos(alphax)/x1/x1-cos(2.0*alphax)/x2/x2)/(1.0/x1/x1-1.0/x2/x2)/cos(alphax))/x1,2)+pow(((2.0*sin(alphax)*sin(alphax)/x2/x2+cos(2.0*alphax)/x1/x1)/(1.0/x1/x1-1.0/x2/x2)/sin(alphax))/x2,2));
    //limity_y=limity_x*(cos(alphax)*cos(alphax)/x1/x1+sin(alphax)*sin(alphax)/x2/x2)/(1.0/x1/x1-1.0/x2/x2)/sin(alphax)/cos(alphax);
    limitx_x=abs(sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax)))/(1+1/(tan(alphax)*tan(alphax)))));
    limitx_y=abs(((1/tan(alphax)*sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2/(tan(alphax)*tan(alphax))))-abs(limitx_x*cos(alphax)))/sin(alphax));
    /*    alphax=M_PI/2-alphax;
          limity_y=abs(sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax)))/(1+1/(tan(alphax)*tan(alphax)))));
          limity_x=abs(((1/tan(alphax)*sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2/(tan(alphax)*tan(alphax))))-abs(limitx_x*cos(alphax)))/sin(alphax));
     */
    limity_y=abs(sqrt((x2*x2+x1*x1*tan(alphax)*tan(alphax))/(1+tan(alphax)*tan(alphax))));
    limity_x=abs(((tan(alphax)*sqrt((x2*x2+x1*x1*(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2*(tan(alphax)*tan(alphax))))-abs(limity_y*sin(alphax)))/cos(alphax));
  }
  else if(alpx==0)
  {
    limitx_x=x1;
    limitx_y=0;
    limity_y=x2;
    limity_x=0;
  }
  else
  {
    limitx_x=x2;
    limitx_y=0;
    limity_y=x1;
    limity_x=0;
  }
  /*
     if(myid==0)
     {
     cout<<"x1re = "<<limity_x<<" , px1re= "<<limity_y<<endl;
     cout<<"x2re = "<<limitx_x<<" , px2re= "<<limitx_y<<endl;
     cout<<limitx_x*limitx_y<<"  "<<limity_x*limity_y<<endl;
     }
   */
  ax=sqrt(x1/x2*pow(cos(alphax),2)+x2/x1*pow(sin(alphax),2));
  axs=(1/ax/2)*((x1/x2)-(x2/x1))*sin(2*alphax);
  ay=sqrt(y1/y2*pow(cos(alphay),2)+y2/y1*pow(sin(alphay),2));
  ays=(1/ay/2)*(y1/y2-y2/y1)*sin(2*alphay);
  emitx=x1*x2;
  //emitx=sqrt(limitx_x*limity_y*limitx_x*limity_y-limitx_x*limitx_y*limitx_x*limitx_y);
  emity=y1*y2;
  emitratio=emitx/emity;
  F=emitx;
  //  cout<<"F! "<<F<<endl;
  for(int num=0;num<particlenumber;num++)
  {
    zetax=sqrt((double)rand() / RAND_MAX*F);
    zetay=sqrt((F-zetax*zetax)/emitratio);
    anglebetax=(double)rand() / RAND_MAX*2*M_PI;
    anglebetay=(double)rand() / RAND_MAX*2*M_PI;
    x=zetax*ax*cos(anglebetax);
    xs=zetax*(axs*cos(anglebetax)-sin(anglebetax)/ax);
    y=zetay*ay*cos(anglebetay);
    ys=zetay*(ays*cos(anglebetay)-sin(anglebetay)/ay);
    z=((double)rand() / RAND_MAX-0.5)*1;
    zs=gaussrand()*1;
    par<<setw(10)<<x<<"   "<<setw(10)<<xs<<"   "<<setw(10)<<y<<"   "<<setw(10)<<ys<<"   "<<setw(10)<<z<<"   "<<setw(10)<<zs<<endl;
  }
  par.close();
}


void KVdistribution6T(const vector<double> &twiss,double emitx,double emity,double emitz,int particlenumber,char *p)
{
  if(twiss.size()<6)
    return;
  double sigx	= (twiss[0]*emitx)/4;
  double sigpx	= ((twiss[1]*twiss[1]+1)/twiss[0]*emitx)/4;
  double sigxpx	= sqrt(sigx*sigpx-emitx*emitx/16);
  if(twiss[1]==0)
    sigxpx=0;
  double sigy	= (twiss[2]*emity)/4;
  double sigpy	= ((twiss[3]*twiss[3]+1)/twiss[2]*emity)/4;
  double sigypy	= sqrt(sigy*sigpy-emity*emity/16);
  if(twiss[3]==0)
    sigypy=0;
  double sigz	= (twiss[4]*emitz)/4;
  double sigpz	= ((twiss[5]*twiss[5]+1)/twiss[4]*emitz)/4;
  double sigzpz	= sqrt(sigz*sigpz-emitz*emitz/16);
  if(twiss[5]==0)
    sigzpz=0;
  //  cout<<"sig : "<<sigx<<"  "<<sigpx<<"  "<<sigxpx<<endl;
  KVdistribution6S(sigx,sigpx,sigxpx,sigy,sigpy,sigypy,sigz,sigpz,sigzpz,particlenumber,p);
}

void KVdistribution6S(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,double sig_z,double sig_pz,double sig_zpz,int particlenumber,char *p)
{
  if(sig_x*sig_px<sig_xpx*sig_xpx||sig_y*sig_py<sig_ypy*sig_ypy||sig_z*sig_pz<sig_zpz*sig_zpz)
  {
    cout<<"error sigma minu"<<endl;
    return;
  }
  MatrixXd CM = MatrixXd::Zero(2,2);
  //  MatrixXd DM = MatrixXd::Zero(3,3);
  VectorXd eivals;
  double x1=2*sig_xpx/sqrt(sig_px);
  double px1=2*sqrt(sig_px);
  double x2=2*sqrt(sig_x);
  double px2=2*sig_xpx/sqrt(sig_x);
  //  double Ax=px2/x1/(x2*px1-x1*px2);
  double Bx=-2*px2/px1/(x2*px1-x1*px2);
  double Cx=x2/px1/(x2*px1-x1*px2);
  double Ax=(1+px2*px2*Cx)/x2/x2;
  CM<<Ax,   Bx/2,
    Bx/2, Cx;
  /*
     DM<<Ax,   Bx/2, 0,
     Bx/2, Cx  , 0,
     0  ,  0  ,  1;
   */
  eivals = CM.selfadjointView<Lower>().eigenvalues();
  if(sig_x<sig_px)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ax=sqrt(abs(1/eivals(0)));
  double bx=sqrt(abs(1/eivals(1)));
  double phix;
  if(Ax!=Cx)
    phix=1.0/2*atan(Bx/(Ax-Cx))*180/M_PI;
  else
    phix=0;
  double emix=ax*bx;
  
  double y1=2*sig_ypy/sqrt(sig_py);
  double py1=2*sqrt(sig_py);
  double y2=2*sqrt(sig_y);
  double py2=2*sig_ypy/sqrt(sig_y);
  double By=-2*py2/py1/(y2*py1-y1*py2);
  double Cy=y2/py1/(y2*py1-y1*py2);
  double Ay=(1+py2*py2*Cy)/y2/y2;
  CM<<Ay,   By/2,
    By/2, Cy;
  eivals = CM.selfadjointView<Lower>().eigenvalues();
  if(sig_y<sig_py)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ay=sqrt(1/eivals(0));
  double by=sqrt(1/eivals(1));
  double phiy;
  if(Ay!=Cy)
    phiy=1.0/2*atan(By/(Ay-Cy))*180/M_PI;
  else
    phiy=0;
  double emiy=ay*by;
  
  double z1=2*sig_zpz/sqrt(sig_pz);
  double pz1=2*sqrt(sig_pz);
  double z2=2*sqrt(sig_z);
  double pz2=2*sig_zpz/sqrt(sig_z);
  double Bz=-2*pz2/pz1/(z2*pz1-z1*pz2);
  double Cz=z2/pz1/(z2*pz1-z1*pz2);
  double Az=(1+pz2*pz2*Cz)/z2/z2;
  CM<<Az,   Bz/2,
    Bz/2, Cz;
  eivals = CM.selfadjointView<Lower>().eigenvalues();
  if(sig_z<sig_pz)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double az=sqrt(1/eivals(0));
  double bz=sqrt(1/eivals(1));
  double phiz;
  if(Az!=Cz)
    phiz=1.0/2*atan(Bz/(Az-Cz))*180/M_PI;
  else
    phiz=0;
  double emiz=az*bz;
  
  /*
     cout<<"x1 = "<<x1<<" , px1= "<<px1<<endl;
     cout<<"x2 = "<<x2<<" , px2= "<<px2<<endl;
     cout<<Ax<<"  "<<Bx<<"  "<<Cx<<"  ABC"<<endl;
     cout<<phix<<"  phi  "<<phiy<<endl;
     cout<<emix<<"  "<<ax<<"  "<<emiy<<"  "<<ay<<"  "<<phix<<"  "<<phiy<<"  "<<particlenumber<<"  "<<p<<endl;;
     cout<<"emit x exp : "<<4*sqrt(sig_x*sig_px-sig_xpx*sig_xpx)<<endl;
   */
  KVdistribution6A(emix,ax,emiy,ay,emiz,az,phix,phiy,phiz,particlenumber,p);
  /*
     particleread(p);
     cout<<beam.getsigmax()<<endl;
     cout<<beam.getsigmadx()<<endl;
     cout<<beam.getsigmaxdx()<<endl;
   */
}

void KVdistribution6A(double emitx,double betax,double emity,double betay,double emitz,double betaz,double alpx,double alpy,double alpz,int particlenumber,char *p)
{
  ofstream par(p);
  double co1,co2;
  double ax,axs,ay,ays,az,azs,emitratio,emitratio2,F;
  double x,xs,y,ys,z,zs;
  double anglebetax,anglebetay,anglebetaz,zetax,zetay,zetaz;
  double alphax=alpx/180*M_PI;
  double alphay=alpy/180*M_PI;
  double alphaz=alpz/180*M_PI;
  double x1=betax;
  double x2=emitx/betax;
  double y1=betay;
  double y2=emity/betay;
  double z1=betaz;
  double z2=emitz/betaz;
  double limity_x,limity_y,limitx_x,limitx_y;
  if(alpx!=0&&alpx!=90)
  {
    //limity_x=sqrt(1.0/pow(((2.0*cos(alphax)*cos(alphax)/x1/x1-cos(2.0*alphax)/x2/x2)/(1.0/x1/x1-1.0/x2/x2)/cos(alphax))/x1,2)+pow(((2.0*sin(alphax)*sin(alphax)/x2/x2+cos(2.0*alphax)/x1/x1)/(1.0/x1/x1-1.0/x2/x2)/sin(alphax))/x2,2));
    //limity_y=limity_x*(cos(alphax)*cos(alphax)/x1/x1+sin(alphax)*sin(alphax)/x2/x2)/(1.0/x1/x1-1.0/x2/x2)/sin(alphax)/cos(alphax);
    limitx_x=abs(sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax)))/(1+1/(tan(alphax)*tan(alphax)))));
    limitx_y=abs(((1/tan(alphax)*sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2/(tan(alphax)*tan(alphax))))-abs(limitx_x*cos(alphax)))/sin(alphax));
    /*    alphax=M_PI/2-alphax;
          limity_y=abs(sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax)))/(1+1/(tan(alphax)*tan(alphax)))));
          limity_x=abs(((1/tan(alphax)*sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2/(tan(alphax)*tan(alphax))))-abs(limitx_x*cos(alphax)))/sin(alphax));
     */
    limity_y=abs(sqrt((x2*x2+x1*x1*tan(alphax)*tan(alphax))/(1+tan(alphax)*tan(alphax))));
    limity_x=abs(((tan(alphax)*sqrt((x2*x2+x1*x1*(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2*(tan(alphax)*tan(alphax))))-abs(limity_y*sin(alphax)))/cos(alphax));
  }
  else if(alpx==0)
  {
    limitx_x=x1;
    limitx_y=0;
    limity_y=x2;
    limity_x=0;
  }
  else
  {
    limitx_x=x2;
    limitx_y=0;
    limity_y=x1;
    limity_x=0;
  }
  /*
     if(myid==0)
     {
     cout<<"x1re = "<<limity_x<<" , px1re= "<<limity_y<<endl;
     cout<<"x2re = "<<limitx_x<<" , px2re= "<<limitx_y<<endl;
     cout<<limitx_x*limitx_y<<"  "<<limity_x*limity_y<<endl;
     }
   */
  ax=sqrt(x1/x2*pow(cos(alphax),2)+x2/x1*pow(sin(alphax),2));
  axs=(1/ax/2)*((x1/x2)-(x2/x1))*sin(2*alphax);
  ay=sqrt(y1/y2*pow(cos(alphay),2)+y2/y1*pow(sin(alphay),2));
  ays=(1/ay/2)*(y1/y2-y2/y1)*sin(2*alphay);
  az=sqrt(z1/z2*pow(cos(alphaz),2)+y2/y1*pow(sin(alphaz),2));
  azs=(1/az/2)*(z1/z2-z2/z1)*sin(2*alphaz);
  emitx=x1*x2;
  //emitx=sqrt(limitx_x*limity_y*limitx_x*limity_y-limitx_x*limitx_y*limitx_x*limitx_y);
  emity=y1*y2;
  emitz=z1*z2;
  emitratio=emitx/emity;
  emitratio2=emitx/emitz;
  F=emitx;
  //  cout<<"F! "<<F<<endl;
  for(int num=0;num<particlenumber;num++)
  {
    co1=(double)rand() / RAND_MAX;
    co2=(double)rand() / RAND_MAX;
/*    if(co1>co2)
    {
      double t=co1;
      co1=co2;
      co2=t;
    }*/
    zetax=sqrt(co1*F);
    zetay=sqrt((1-co1)*F/emitratio);
    zetaz=sqrt((co2)*F/emitratio2);
    anglebetax=(double)rand() / RAND_MAX*2*M_PI;
    anglebetay=(double)rand() / RAND_MAX*2*M_PI;
    anglebetaz=(double)rand() / RAND_MAX*2*M_PI;
    x=zetax*ax*cos(anglebetax);
    xs=zetax*(axs*cos(anglebetax)-sin(anglebetax)/ax);
    y=zetay*ay*cos(anglebetay);
    ys=zetay*(ays*cos(anglebetay)-sin(anglebetay)/ay);
    z=zetaz*az*cos(anglebetaz);
    zs=zetaz*(azs*cos(anglebetaz)-sin(anglebetaz)/az);
    par<<setw(10)<<x<<"   "<<setw(10)<<xs<<"   "<<setw(10)<<y<<"   "<<setw(10)<<ys<<"   "<<setw(10)<<z<<"   "<<setw(10)<<zs<<endl;
  }
  par.close();
}

void WBdistribution(double emitx,double betax,double emity,double betay,double alpx,double alpy,int particlenumber)
{
  ofstream par("particles.dat");
  double x1=2*sqrt(emitx*betax);
  double x2=2*sqrt(emitx/betax);
  double y1=2*sqrt(emity*betay);
  double y2=2*sqrt(emity/betay);
  double ax,axs,ay,ays,emitratio,F,G;
  double x,xs,y,ys;
  double anglebetax,anglebetay,zetax,zetay;
  double alphax=alpx/180*M_PI;
  double alphay=alpy/180*M_PI;
  ax=sqrt(x1/x2*pow(cos(alphax),2)+x2/x1*pow(sin(alphax),2));
  axs=(1/ax/2)*((x1/x2)-(x2/x1))*sin(2*alphax);
  ay=sqrt(y1/y2*pow(cos(alphay),2)+y2/y1*pow(sin(alphay),2));
  ays=(1/ay/2)*(y1/y2-y2/y1)*sin(2*alphay);
  emitx=x1*x2;
  emity=y1*y2;
  emitratio=emitx/emity;
  for(int num=0;num<particlenumber;num++)
  {
    G=(double)rand() / RAND_MAX;
    F=3*emitx/2*G;
    zetax=sqrt((double)rand() / RAND_MAX*F);
    zetay=sqrt((F-zetax*zetax)/emitratio);
    anglebetax=(double)rand() / RAND_MAX*2*M_PI;
    anglebetay=(double)rand() / RAND_MAX*2*M_PI;
    x=zetax*ax*cos(anglebetax);
    xs=zetax*(axs*cos(anglebetax)-sin(anglebetax)/ax);
    y=zetay*ay*cos(anglebetay);
    ys=zetay*(ays*cos(anglebetay)-sin(anglebetay)/ay);
    par<<setw(10)<<x<<"   "<<setw(10)<<xs<<"   "<<setw(10)<<y<<"   "<<setw(10)<<ys<<"   "<<setw(10)<<0<<"   "<<setw(10)<<0<<endl;
  }
  par.close();
}

void distribution(double emitx,double betax,double emity,double betay,double rx,double ry,double rxy,double rxtyt,double rxyt,double rxty,int particlenumber,char *p)
{
  //x,x'   y,y'
  emitx=emitx;
  emity=emity;
  double a=sqrt(emitx*betax);
  double b=sqrt(emitx/betax);
  double c=sqrt(emity*betay);
  double d=sqrt(emity/betay);
  ofstream par(p);
  ofstream par1("particletests.dat");
  double U1,U2,U3,U4,X1,X2,X3,X4;
  MatrixXd C=MatrixXd::Identity(4,4);
  MatrixXd m = MatrixXd::Identity(4,4);
  m(0,0)=a*a;
  m(1,1)=b*b;
  m(0,1)=rx*a*b;
  m(1,0)=rx*a*b;
  m(2,2)=c*c;
  m(3,3)=d*d;
  m(2,3)=ry*c*d;
  m(3,2)=ry*c*d;
  m(0,2)=rxy*a*c;
  m(2,0)=rxy*a*c;
  m(0,3)=rxyt*a*d;
  m(3,0)=rxyt*a*d;
  m(1,2)=rxty*b*c;
  m(2,1)=rxty*b*c;
  m(1,3)=rxtyt*b*d;
  m(3,1)=rxtyt*b*d;		//or m(3,1)=rxtyt*rxtyt*b*d;?
  C = m.llt().matrixL();
  //cout<<C<<endl;
  for(int number_of_particles=0;number_of_particles<particlenumber;number_of_particles++)//
  {
    double x[4]={0,0,0,0};
    double mu[4]={0,0,0,0};//:expectation
    double u[4];
    for(int i=0;i<4;i++)
    {
      u[i]=gaussrand();
      par1<<setw(10)<<u[i]<<"   ";
    }
    par1<<setw(10)<<0<<"   "<<setw(10)<<0<<endl;
    for(int k=0;k<4;k++)
    {
      for(int i=0;i<=k;i++)
      {
        x[k]=x[k]+C(k,i)*u[i];
      }
      x[k]=x[k]+mu[k];
      par<<setw(10)<<x[k]<<"   ";
    }
    double z=(double)rand() / RAND_MAX*1;
    double zs=gaussrand()*1;
    par<<setw(10)<<z<<"   "<<setw(10)<<zs<<endl;
  }
  par.close();
  par1.close();
}
