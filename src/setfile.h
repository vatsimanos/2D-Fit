#ifndef SETFILE_H
#define SETFILE_H

#include <iostream>	//.h weg
#include <fstream>	//.h weg
#include <ctype.h>
#include <string.h>
#include <cstdlib>
#include "patchio.h"
#include "timeseries.h"
#include "channel.h"
#include "err.h"
#include "daten.h"
#include "rausch.h"
#include <chrono>


class TSetfile
{
 private:
  double delta_t;
  double sigma;
  double zero_current;
  double f_3dB_digfi;
  int no_channel_types;
  int *no_channels_type;
  int no_channels_total;
  int line_counter;
  unsigned long no_samples;
  
  //unsigned short int *noiseseries;
  //double *rauschcopy;
  
  //BEGIN KARSTENFIT (Zu "public" ge�ndert)
  //char datafile[256];
  //END KARSTENFIT  
  
  char rausch_file[256];
  char stepresponse_file[256];
  Tchannel_type	**channel_type;
  
  int **int_matrix(int rows, int colums);
  double **double_matrix(int rows, int colums);
  char *next_line(FILE *Fsetfile);
  double **read_dMatrix(FILE *Fsetfile, int *Zeilen, int *Spalten);
  int **read_iMatrix(FILE *Fsetfile, int* Zeilen, int *Spalten);

 public:
 
 //BEGIN KARSTENFIT (von "private" zu "public")
char datafile[256];
//END KARSTENFIT

  Terror err;

  TSetfile();

  ~TSetfile();
  
  int iterations;  //Tobias: f�r den 2D-Fit
  double total_nr_events; // Tobias f�r den 2D-Fit
  
  void init(double deltat,double sigm, double zerocurr, double f3dB, int nochanneltypes, int *nochanneltype, int nosamples, char *dataf, char *rauschfile, char *stepresp, Tchannel_type **channeltype); 

  double read_setfile(char *setfilename);
  
  long int get_no_samples();  //Tobias

  double get_sigma();
    
  void set_sigma (double noise);				//Tobias 2021 set noise lvl 	
  
  double get_zero_current();
  
  double get_single_channel_current();                      //Tobias
  
  int get_no_rconst();		                    //Tobias: die Anzahl der Ratenkonstanten in der Matrix <> 0
  
  double get_rconst(unsigned int nr_rconst);                //Tobias: die Ratenkonstante n<>0 (aus der Matrix erst x, dann y) wird ausgelesen;

  void set_rconst(unsigned int nr_rconst, double rconst);   //Tobias: die Ratenkonstante n<>0 (aus der Matrix erst x, dann y) wird gesetzt;
  
  double get_rconstxy (int x, int y);			    //Tobias: liefert die Ratenkonstante x,y aus der Matrix	
  
  int get_matrix_dimension();				    //Tobias
		
  double get_delta_t();

  double get_f3dB();

  int get_no_channels_total();

  int get_no_channel_types();

  char *get_timeseries_filename(char *Filename);

  char *get_noiseseries_filename(char *Filename);

  char *get_rausch_filename(char *Filename);

  char *get_stepresponse_filename(char *Filename);  

  void write_setfile(char *name);

  void clear();

  TTimeseries *simulate_timeseries(TTimeseries *ts, TRausch  *rausch, char *signal_filename=(char *)"", double dT_filter=0);

};

#endif /* SETFILE_H */








