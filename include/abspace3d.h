#ifndef ABSPACE3D_H
#define ABSPACE3D_H

#include        "Devices.h"
#include        "Monitor.h"
#include        "Field.h"
#include        "Lattice.h"
#include        "Bunch.h"
#include        "util.h"
#include        "ReadFile.h"
//#include        "picGrid.h"
#include        <vector>


class ABSpace3D
{
 public:
  ABSpace3D();
  ~ABSpace3D();

  /* Preprocessing */
  void prepose();

  /* Allocate memory */
  
  void allocatememory();
  void initfields();
  void initField(DomainContent & dc);
  void gridsize(DomainContent &dc);

  /* General functions */

  void interpolatepbzf(double *BZ);
  int findxid(int id);
  int findyid(int id);
  int findzid(int id);

  /* General variables and handlers */

  std::vector<Monitor>  *monitor;
  Lattice lattice;
  std::vector<Bunch> bunchs;
//  SCsolver * solve;

//  Bunch bunch,bunchold;
  ParamFileContent projectContent;
//  BunchInitParam_Beampath initParam;
//  int inInitParam_Beampath(std::string & file_name);

  double **pfds;
  // Devices list
  Devices devs;

  // Synchronized particle

  Phase *synpart;

};


#endif// ABSPACE3D_H
