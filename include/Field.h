#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;

class Field
{
public:
  double getField(int,int,int);
  double getField(int);
  void read(char* );
  void addField(double);
  
protected:
  vector<double> 	 field;
  
};


#endif
