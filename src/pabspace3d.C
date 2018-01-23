
#include "pabspace3d.h"
//#include "driver3d.h"
#include <memory>
#include <iostream>

using std::cout;
using std::endl;
using std::auto_ptr;

PABspace3D::~PABspace3D(){ };

void PABspace3D::setComm(MPI_Comm comm_){
  Parallel::setComm(comm_);
//  bunch.setComm(comm_);
  //  bunchold.setComm(comm_);
};
  

void PABspace3D::initspace(){

  return;

};


void PABspace3D::pprepose(){

#if NODEBUG
  printf0("in pprepose\n");
#endif

  allocatememory();

  prepose();

  return;

};

