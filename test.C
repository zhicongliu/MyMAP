#include "Mapmain.h"

int main(int argc,char **argv)
{
  double t0,t1,t2,t3;
  MPI_Init(&argc,&argv);
  t0=MPI_Wtime();
  int myid,numprocs,namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Get_processor_name(processor_name,&namelen);  
  srand(time(0));
  cout<<setprecision(8)<<fixed;
  Mapmain Acc;
  LatticeMap Ltemp;
  Ltemp.prelattice1("pho50090","pho50090e");
  /*
  if(myid==0)
  {
    //distribution(10.931334*0.049039921,10.931334,1.2421721*0.00059654248,1.2421721,0.5,0);
    //double a[]={10,0,1.2421721,1.0421721};
    double a[]={10.38592243,-0.00000000,0.77894515,-0.00000000,1000,0};
    vector<double> twiss1(a,a+6);
    KVdistribution6T(twiss1,10,10,90,10000);
    KVdistributionT(twiss1,10,10,10000,"particle2.dat");
    //KVdistributionS(1e-5,1e-5,0,100,200,2,10000);
    Acc.preparticle("particles.dat");
  }
  */
  MPI_Barrier(MPI_COMM_WORLD);
  t1=MPI_Wtime();
  Acc.particleread("particles.dat");//include energy initial
  Acc.latticeinput("lattice.txt",1,0);//need enengy
  
  /*
  GapMap a;
  a.MakeMap(5,10,1,1e9,1024);
  */
  
  /*
  LatticeMap e1,e2;
  e1.read("lattice1.txt",1.3);
  MPI_Barrier(MPI_COMM_WORLD);
  e2.read("lattice2.txt",1.3);
  MPI_Barrier(MPI_COMM_WORLD);
  double v11[]={10.38607658,0.0,0.77897733,0.0};
  vector<double> v1(v11,v11+5);
  double v22[]={10.38607658,0.0,0.77897733,0.0};
  vector<double> v2(v22,v22+5);
  Acc.Glatticescan(e1,e2,v1,v2,2);
  */
  
  /*
  LatticeMap e1,e2;
  e1.read("lattice1.txt");
  e2.read("lattice2.txt");
  Acc.Glatticescan(e1,e2,1.898811,2.898811);
  */
  
  /*
  LatticeMap e1,e2;
  e1.read("lattice1.txt");
  e2.read("lattice2.txt");
  double v11[]={10.931333,0.0,1.242172,0.0};
  vector<double> v1(v11,v11+5);
  Acc.Glatticescan(e1,e2,v1,39.7,75.2);
*/
  
  
  //Acc.twissscan(7,5,1e-5,1.24,1e-5);
  //Acc.Gtwissscan(8,1e-8,50,-0.1,0.1,1e-8,10,-0.1,0.1);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0) 
  {
  Acc.preparticle("particles.dat");
  Acc.caculate();
  Acc.twisstrack(Acc.getbetax(),Acc.getalphax(),Acc.getbetay(),Acc.getalphay());
  Acc.ptrackall(0,"particles.dat");
  //Acc.twisstrack(3.5929,-0.1939,1.61326,-0.45533);
  //Acc.twisstrack(11.254229 ,  0.078000  , 1.261411   ,0.043102);
  //Acc.ptrackonce();
  }
  
  //Acc.ptrackall(0);
  MPI_Barrier(MPI_COMM_WORLD);
  t2=MPI_Wtime();
  Acc.postparticle();
  t3=MPI_Wtime();
  if(myid==0)
  {
    cout<<"pre time = "<<t1-t0<<endl;
    cout<<"track time = "<<t2-t1<<endl;
    cout<<"post time = "<<t3-t2<<endl;
  }
  MPI_Finalize();
  return 0;
}
