#include "Field.h"

double Field::getField(int z)
{
  return field[z];
}

void Field::read(char *p)
{
  ifstream infile(p);
  string line;
  while(infile.peek()!=EOF){
    getline(infile,line);
    if(line.size()>0)
    {
      field.push_back(atof(line.c_str())*1e0);
    }
  }
}

void Field::addField(double a)
{
  field.push_back(a);
}



