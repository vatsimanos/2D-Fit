


#ifndef TTIMESERIES_H
#define TTIMESERIES_H

//#include <stdlib.h>
#include <fstream>	//.h weg
#include <cstdlib>
#include "patchio.h"
#include "err.h"
//#include "glib.h"
#include "daten.h"


#define ASCII 0
#define BINARY 1


class TTimeseries
{
 private:
  unsigned short int *Value;
  unsigned long int Num_Samples;
  double delta_t;
  char *FileName;
  double f_3dB;
 
 public:
  Terror err;

  TTimeseries();

  ~TTimeseries();

  char *get_filename(char *Filename);

  unsigned long int get_num_samples();

  void set_num_samples(unsigned long int num_samples);

  void set_value(unsigned short int *value,unsigned long int length);

  void set_delta_t(double deltat);

  void set_filename(char *name);

  void set_f_3dB(double f3dB);

  inline unsigned long int timeseries(unsigned long int location);  

  int load(char *FName);

  void save(char *FName);
  
  void save_amplhis(char *FName);

  void clear();

  //void simulate_timeseries(TSetfile *setfile);

  inline unsigned short int& operator[](const int& index);

  friend std::ostream& operator<< (std::ostream& os, const TTimeseries& ts);

};

// friend functions must be global
std::ostream& operator<< (std::ostream& os, const TTimeseries& ts);

// inline functions must be declared here

inline unsigned long int  TTimeseries::timeseries(unsigned long int location)
{
  if (Value!=NULL)
    return ((Value[location]) & 0x0fff); //the upper 4 Bits contain some marker
  else
    {
      err.error(NO_VALUES,(char *)"TTimeseries::operator[]",FileName);
      return 0;
    }
}

inline unsigned short int& TTimeseries::operator[](const int& index)
{
  if (Value==NULL)
    err.error(NO_VALUES,(char *)"TTimeseries::operator[]",FileName);
  else if ((index<0)||(index>(int)Num_Samples)){
    char buffer[32];
    sprintf(buffer,"index=%i",index);
    err.error(OUT_OF_RANGE,(char *)"",buffer);}
  //err.warning(VALUE_NOT_MASKED,"TTimeseries::operator[]");
  return Value[index];
}



#endif /* TTIMESERIES_H */ 







