#include "LatticeMap.h"

LatticeMap::~LatticeMap()
{
}
LatticeMap::LatticeMap()
{
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
}
LatticeMap::LatticeMap(const LatticeMap &e)
{
  ELE=e.ELE;
  //  cout<<"copy"<<endl;
  Mat=e.Mat;
  fitness=e.fitness;
  EBeta=e.EBeta;
  EGamma=e.EGamma;
  this->elemtemp=e.elemtemp;
  this->FQ=e.FQ;
  this->BM=e.BM;
  this->DL=e.DL;
  this->SN=e.SN;
  this->RFQ=e.RFQ;
  this->rfqtypepre=e.rfqtypepre;
  this->rfqtypenow=e.rfqtypenow;
  this->rfqtypenext=e.rfqtypenext;
  // this->infile=e.infile;
  this->NDivide=e.NDivide;
  this->Mfail=e.Mfail;
  this->particle_count=e.particle_count;
  this->totalparticalnumber=e.totalparticalnumber;
  this->i=e.i;
  this->j=e.j;
  this->k=e.k;
  this->myid=e.myid;
  this->numprocs=e.numprocs;
  this->positionpre=e.positionpre;
  this->positionnow=e.positionnow;
  this->positionnext=e.positionnext;

}

LatticeMap& LatticeMap::operator =(const LatticeMap& e)
{
  ELE=e.ELE;
  //  cout<<"assign"<<endl;
  Mat=e.Mat;
  fitness=e.fitness;
  EBeta=e.EBeta;
  EGamma=e.EGamma;
  this->elemtemp=e.elemtemp;
  this->FQ=e.FQ;
  this->BM=e.BM;
  this->DL=e.DL;
  this->SN=e.SN;
  this->RFQ=e.RFQ;
  this->rfqtypepre=e.rfqtypepre;
  this->rfqtypenow=e.rfqtypenow;
  this->rfqtypenext=e.rfqtypenext;
  // this->infile=e.infile;
  this->NDivide=e.NDivide;
  this->Mfail=e.Mfail;
  this->particle_count=e.particle_count;
  this->totalparticalnumber=e.totalparticalnumber;
  this->i=e.i;
  this->j=e.j;
  this->k=e.k;
  this->myid=e.myid;
  this->numprocs=e.numprocs;
  this->positionpre=e.positionpre;
  this->positionnow=e.positionnow;
  this->positionnext=e.positionnext;    
  return *this;
}

bool LatticeMap::operator==(const LatticeMap& l)
{
  if(ELE.size()!=l.ELE.size())
    return false;
  for(int i=0;i<ELE.size();i++)
  {
    if(ELE[i]!=l.ELE[i])
      return false;
  }
  return true;
}


void LatticeMap::callelement(const Element &v,ifstream &infile)
{
  if(v.name=="drift")
  {
    DriftMap drift;
    drift.setlength(v.p[2]);
    addEGamma(EGamma.back());
    drift.setgamma(EGamma.back());
    drift.Map();
    Mat.push_back(drift.getmap());
    //cout<<i<<"  drift \n"<<Mat[i]<<endl;
  }
  else if(v.name=="quadm")
  {
    QuadrupoleMap qf;
    qf.setlength(v.p[2]);
    qf.setgradient(v.p[4]);
    addEGamma(EGamma.back());
    qf.setgamma(EGamma.back());
    qf.Map();
    Mat.push_back(qf.getmap());
    //cout<<i<<"  quadm \n"<<Mat[i]<<endl;
  }
  else if(v.name=="bendm")
  {
    BendingMap bendm;
    bendm.setlength(v.p[2]);
    bendm.setradius(v.p[4]);
    addEGamma(EGamma.back());
    bendm.setgamma(EGamma.back());
    bendm.Map();
    Mat.push_back(bendm.getmap());
  }
  else if(v.name=="solen")
  {
    SolenoidMap solenoid(v.p[2],  v.p[4],  EGamma.back());
    addEGamma(EGamma.back());
    solenoid.Map();
    Mat.push_back(solenoid.getmap());
  }
  else if(v.name=="gap")
  {
    //            if(myid==0)
    //              cout<<"BEFORE   "<<EGamma.back()<<endl;
    addEGamma((EGamma.back()*_MASS*C_light*C_light+v.p[4]*1e6*cos(v.p[5]*M_PI/180.0)*abs(_Q))/(_MASS*C_light*C_light));
    //            if(myid==0)
    //              cout<<"AFTER    "<<EGamma.back()<<endl;
    if(v.p[1]==1)
    {
      GapMap gap(v.p[2],v.p[4],v.p[5],v.p[6],EGamma.back());
      gap.Map();
      Mat.push_back(gap.getmap());
    }
    else if(v.p[1]==2)
    {
      char temp[99];
      strcpy(temp,v.path.c_str());
      //cout<<temp<<endl;
      GapMap gap(v.p[2],temp,v.p[5],v.p[6],EGamma.back());
      Mat.push_back(MatrixXd::Identity(6,6));
      gap.getEf();
    }
    else
    {
      cout<<"GAP TYPE error!"<<endl;
      Mat.push_back(MatrixXd::Identity(6,6));
    }
    //cout<<i<<"  drift \n"<<Mat[i]<<endl;
  }
  else if(v.name=="rfqce")
  {
    //(rfqce type length radius voltage phase )
    RFQMap rfq;
    rfq.clear();
    addEGamma(EGamma.back());
    rfq.initial(v.p[2],v.p[4],v.p[5],v.p[6],EGamma.back());
    string rfqtemp;
    infile.seekg(positionpre);
    getline(infile,rfqtemp);
    vector<string> parainput;
    if(rfqtemp.size()>5)
    {
      StringSplit(rfqtemp,' ',parainput);
      if(parainput[0]=="rfqce")
      {
        rfqtypepre=atof(parainput[1].c_str());
        //cout<<"typepre=\t  "<<para[1].c_str()<<endl;
      }
    }
    infile.seekg(positionnow);
    getline(infile,rfqtemp);
    if(rfqtemp.size()>5)
    {
      StringSplit(rfqtemp,' ',parainput);
      if(parainput[0]=="rfqce")
      {
        rfqtypenow=atof(parainput[1].c_str());
        //cout<<"typenow=\t  "<<para[1].c_str()<<endl;
      }
    }
    infile.seekg(positionnext);
    getline(infile,rfqtemp);
    if(rfqtemp.size()>5)
    {
      StringSplit(rfqtemp,' ',parainput);

      if(parainput[0]=="rfqce")
      {
        rfqtypenext=atof(parainput[1].c_str());
        //cout<<"typenext=\t  "<<para[1].c_str()<<endl;
      }
    }
    infile.seekg(positionnext);

    rfq.settype(rfqtypepre,rfqtypenow,rfqtypenext);
    rfq.Map();
    Mat.push_back(rfq.getmap());
    //cout<<i<<"  rfq \n"<<Mat[i]<<endl;
    //rfq.print();
  }
  else
  {
    ELE.pop_back();
    i--;
    Mfail++;
  }
}

void LatticeMap::modify(const int &num,const Element &v)
{
  elemtemp=ELE[num];
  if(num<Mat.size())
  {
    ELE[num]=v;
    if(v.name=="drift")
    {
      DriftMap drift;
      drift.setlength(v.p[2]);
      drift.setgamma(EGamma[num]);
      drift.Map();
      Mat[num]=drift.getmap();
    }
    else if(v.name=="quadm")
    {
      QuadrupoleMap qf;
      qf.setlength(v.p[2]);
      qf.setgradient(v.p[4]);
      qf.setgamma(EGamma[num]);
      qf.Map();
      Mat[num]=qf.getmap();
    }
    else if(v.name=="bendm")
    {
      BendingMap bendm;
      bendm.setlength(v.p[2]);
      bendm.setradius(v.p[4]);
      bendm.Map();
      Mat[num]=bendm.getmap();
    }
    else if(v.name=="solen")
    {
      SolenoidMap solenoid;
      solenoid.setlength(v.p[2]);
      solenoid.setfield(v.p[4]);
      solenoid.Map();
      Mat[num]=solenoid.getmap();
    }
    else if(v.name=="rfqce")
    {/*
      //(rfqce type length radius voltage phase )
      RFQMap rfq;
      rfq.clear();
      rfq.initial(v.p[2],v.p[4],v.p[5],v.p[6],EGamma);
      string rfqtemp;
      infile.seekg(positionpre);
      getline(infile,rfqtemp);
      vector<string> parainput;
      if(rfqtemp.size()>5)
      {
      StringSplit(rfqtemp,' ',parainput);
      if(parainput[0]=="rfqce")
      {
      rfqtypepre=atof(parainput[1].c_str());
      //cout<<"typepre=\t  "<<para[1].c_str()<<endl;
      }
      }
      infile.seekg(positionnow);
      getline(infile,rfqtemp);
      if(rfqtemp.size()>5)
      {
      StringSplit(rfqtemp,' ',parainput);
      if(parainput[0]=="rfqce")
      {
      rfqtypenow=atof(parainput[1].c_str());
      //cout<<"typenow=\t  "<<para[1].c_str()<<endl;
      }
      }
      infile.seekg(positionnext);
      getline(infile,rfqtemp);
      if(rfqtemp.size()>5)
      {
      StringSplit(rfqtemp,' ',parainput);

      if(parainput[0]=="rfqce")
      {
      rfqtypenext=atof(parainput[1].c_str());
      //cout<<"typenext=\t  "<<para[1].c_str()<<endl;
      }
      }
      infile.seekg(positionnext);

      rfq.settype(rfqtypepre,rfqtypenow,rfqtypenext);
      rfq.Map();
      Mat[num]=rfq.getmap();
      */
    }
    else
    {
      ELE[num]=elemtemp;
    }
  }
}
void LatticeMap::read(char *p,double EGam,int numberdivided,int elemlimit)
{
  addEGamma(EGam);
  char lp[50]="latticeprecision19920523.liuzc";
  ifstream infile;
  ofstream outlat("latticereadin.dat");
  prelattice(p,lp,numberdivided);
  NDivide=numberdivided;
  int preciseflag=0;
  string line;
  positionpre=0;
  infile.open(lp);
  i=0;
  Mfail=0;
  vector<string> parainput;
  int elementcount=0,elementlimit;
  if(elemlimit==0)
    elementlimit=1e8;
  else
    elementlimit=elemlimit;
  while(infile.peek()!=EOF&&elementcount<=elementlimit){
    elementcount++;
    positionpre=positionnow;
    positionnow=infile.tellg();
    getline(infile,line);
    positionnext=infile.tellg();
    if(line.size()>5)
    {
      StringSplit(line,' ',parainput);
      elemtemp.name=parainput[0];
      elemtemp.path=parainput[4];
      elemtemp.p.clear();
      for(int r=0; r<parainput.size()&&r<12;r++)
      {
        elemtemp.p.push_back(atof(parainput[r].c_str()));
      }
      ELE.push_back(elemtemp);
      callelement(ELE.back(),infile);
      outlat<<line<<"  "<<EGamma.back()<<endl;
    }else 
    {
      i--;
      Mfail++;
    }
    i++;
  }//while end
  if(myid==0)
    cout<<"number of element CAN/CAN'T read in  "<<Mat.size()<<" and "<<Mfail<<" in "<<p<<endl;
  infile.close();
  remove(lp);
}

void LatticeMap::prelattice(const char *in,const char *out,int n)
{
  if(myid==0)
  {
    ifstream infilep(in);
    ofstream outfile(out);
    string line;
    vector<string> parainput;
    while(infilep.peek()!=EOF)
    {
      getline(infilep,line);
      string str = line;
      remove(str.begin(), str.end(), ' ');
      remove(str.begin(), str.end(), '\t');
      if(str.size()>1)
      {
        StringSplit(line,' ',parainput);
        double length_prelattice = atof(parainput[2].c_str());
        if((!(parainput[0]=="gap"))&&length_prelattice>0)
        {
          for(int lnum=0;lnum<n;lnum++)
          {
            for(int i1=0;i1<parainput.size();i1++)
            {
              if(i1!=2)
                outfile<<parainput[i1]<<"  ";
              else
                outfile<<length_prelattice/n<<"  ";
            }
            outfile<<endl;
          }
        }
        else
        {
          for(int i1=0;i1<parainput.size();i1++)
          {
            outfile<<parainput[i1]<<"  ";
          }
          outfile<<endl;
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void LatticeMap::prelattice1(const char *in1,const char *out1)
{
  char in[99],out[99];
  for(int i=50001;i<=50200;i++)
  {
    sprintf(in, "/home/ct/fperturbationevolution/sigma=87KVbeam/pho%d", i);
    sprintf(out, "/home/ct/fperturbationevolution/KVNEW/pho%d", i);
  if(myid==0)
  {
    ifstream infilep(in);
    ofstream outfile(out);
    string line;
    vector<string> parainput;
    double length_prelattice=1e5;
    getline(infilep,line);
    getline(infilep,line);
    while(infilep.peek()!=EOF)
    {
      getline(infilep,line);
      string str = line;
      remove(str.begin(), str.end(), ' ');
      remove(str.begin(), str.end(), '\t');
      remove(str.begin(), str.end(), '\n');
      if(str.size()>2)
      {
        StringSplit(line,' ',parainput);
        double temp=abs(length_prelattice-atof(parainput[1].c_str()));
        if(length_prelattice<1e4&&temp>1e-5)
        {
          outfile<<endl;
        }
        for(int i1=1;i1<parainput.size();i1++)
        {
          if(atof(parainput[1].c_str())==4&&i1==3)
            outfile<<2<<"  ";
          else if(atof(parainput[1].c_str())==-4&&i1==3)
            outfile<<-2<<"  ";
          else
          outfile<<atof(parainput[i1].c_str())<<"  ";
        }
        outfile<<endl;
      length_prelattice = atof(parainput[1].c_str());
      }
    }

  }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


void LatticeMap::print()
{
  for(int ii=0;ii<ELE.size();ii++)
  {
    cout<<ii<<"/"<< ELE.size()<<"  ";
    cout<<ELE[ii].name<<"  ";
    for(int i2=1;i2<ELE[ii].p.size();i2++)
      cout<<ELE[ii].p[i2]<<"  ";
    cout<<endl;
    //cout<<Mat[ii]<<endl;
  }
}
void LatticeMap::addEGamma(const double num)
{
  if(num<1)
  {
    cout<<"EGamma Wrong in lattice"<<endl;
    EGamma.push_back(1);
    EBeta.push_back(0);
    return;
  }
  EGamma.push_back(num);
  EBeta.push_back(sqrt(1-1/num/num));
}

MatrixXd LatticeMap::getMat(int L,double &x,double &dx,double &y,double &dy,double &z,double &dp)
{
  if(ELE[L].name=="gap")
  {
    if(ELE[L].p[1]==1)
    {
      GapMap gap(ELE[L].p[2],ELE[L].p[4],ELE[L].p[5],ELE[L].p[6],EGamma[L]);
      return(gap.getmap(x,dx,y,dy,z,dp));
    }
    else if(ELE[L].p[1]==2)
    {
      char temp[99];
      strcpy(temp,ELE[L].path.c_str());
      GapMap gap(ELE[L].p[2],temp,ELE[L].p[5],ELE[L].p[6],EGamma[L]);
      //cout<<z<<" BEF "<<dp<<endl;
      gap.FieldMethod(x,dx,y,dy,z,dp);
      //cout<<z<<" AFT "<<dp<<endl;
      return MatrixXd::Identity(6,6);
    }
    else
    {
      cout<<"GAP TYPE error!"<<endl;
      return MatrixXd::Identity(6,6);
    }
  }
  else
    return Mat[L];
}








