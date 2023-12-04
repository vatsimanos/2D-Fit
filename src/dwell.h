/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#ifndef _DWELL_H_
#define _DWELL_H_

#include <vector> //.h weg

//#include "glib.h"

#include "matrix.h"
//#include "windows.h"

//Tobias

// ergibt fuer wert aus Intervall einen farbcode 
void ZToRGB(const double z, const double zmin, const double zmax, 
	        unsigned int &R, unsigned int &G, unsigned int &B, bool grey);

class TDwell_lb
{
// abstrakte klasse fuer dwells mit log binning

  // name des objektes, wird an verschiedenen stellen ausgegeben!
  char      dwell_name[50];       
  // intervalle fuer die achsen: wenn [0:0] dann auto range!
  short    dwell_min_close_log;
  short    dwell_max_close_log;     
  short    dwell_min_open_log;   
  short    dwell_max_open_log;   // in log10(x us)
  double   dwell_min_e;          // events: 1d y, 2d z
  double   dwell_max_e;          // wegen normierung bei differenz nicht ganzzahlig 
  // mit bins intervallen pro log
  unsigned short   dwell_bins_per_log;    
  bool  dwell_sqrt_e;         // muss noch in load, store, << und >>!

  virtual double min() = 0;      // minimale anzahl von events
  virtual double max() = 0;      // maximale anzahl von events
                                  // hier werden alle bins verglichen

  virtual void calc(std::vector<TSprung>* vec_ptr, double maxevents = 0) = 0; // spruenge addieren

public:
  // zum auslesen der elemente
  virtual char*    name();		//Tobias 2017 now virtual
  virtual short   min_close_log();	//Tobias 2017 now virtual
  virtual short   max_close_log();	//Tobias 2017 now virtual
  virtual short   min_open_log();	//Tobias 2017 now virtual
  virtual short   max_open_log();	//Tobias 2017 now virtual

  // auch mit autorange, wenn sym, min_z == -max_z
  virtual double  min_e(bool sym = true);	//Tobias 2017 now virtual
  virtual double  max_e(bool sym = true);		//Tobias 2017 now virtual
  virtual unsigned short  bins_per_log();	//Tobias 2017 now virtual
  virtual bool sqrt_e();  //Tobias 2017 now virtual

  // zum setzen der elemente
  void name(char* aname);
  void close_range(short amin, short amax);
  void open_range(short amin, short amax);
  virtual void e_range(double amin, double amax);  // zum setzten der werte //Tobias 2017 now virtual
  virtual void bins_per_log(unsigned short abins);			//Tobias 2017 now virtual
  virtual void sqrt_e(bool asqrt);	//Tobias 2017 now virtual

  // weitere element funktionen
  bool e_auto_range();                   // autorange bei events?

  virtual void init(char* aname, unsigned short abins = 0);

  virtual void init(char* aname,
	            short log_min_close, short log_max_close, // nur ganzzahlige grenzen!
                    short log_min_open, short log_max_open,    
                    unsigned short abins = 0);

// grenzen aus spruengen berechnen
  void calcrange(short &log_min_close, short &log_max_close,
                 short &log_min_open, short &log_max_open);

  virtual void GNUplot(bool dosqrt = true) = 0;
  virtual ~TDwell_lb(){};

  unsigned short bins_open();
  unsigned short bins_close();

  void out(std::ostream& out);
  void in(std::istream& in);

  virtual void store(std::ostream& out);
  virtual void load(std::istream& in);   // fuer jobfile 
};

class TDwell_1d : public TDwell_lb
{
// 1d mit log binning
public:
  std::vector<double> Open;     // hier wird open histogramm gespeichert 
  std::vector<double> Close;    // hier close
  // double wegen normierung bei differenz 
  TDwell_1d(unsigned short abins = 50);

  virtual void init(char* aname, unsigned short abins = 0);

  virtual void init(char* aname,
	                short log_min_close, short log_max_close,
				    short log_min_open, short log_max_open,
				    unsigned short abins = 0);

  virtual void calc(std::vector<TSprung>* vec_ptr, double maxevents = 0);

  virtual double min(); // machen nicht wirklich sinn, solange nicht
  virtual double max(); // Open und Close getrennt sind! 

  TDwell_1d& operator+=(TDwell_1d& Dwell_1d);

  virtual void GNUplot(bool dosqrt = true);

  // abspeichern ascii
  friend std::ostream& operator<<(std::ostream& out, TDwell_1d& Dwell);
  // einlesen ascii
  friend std::istream& operator>>(std::istream& in, TDwell_1d& Dwell);   

  virtual void store(std::ostream& out);
  virtual void load(std::istream& in);   // fuer jobfile 
};

class TDwell_2d : public TDwell_lb
{
// 2d mit log binning
public:
  matrix<double> CO;   // hier wird close - open histogramm gespeichert 
  matrix<double> COD;   //Tobias: Matrix f�r Dependency-Plots
  // double wegen normierung bei differenz 

  TDwell_2d(unsigned short abins = 10);

  virtual void init(char* aname, unsigned short abins = 0);

  virtual void init(char* aname, 
	                short log_min_close, short log_max_close,
  			short log_min_open, short log_max_open,
 			unsigned short abins = 0);

  virtual void calc(std::vector<TSprung>* vec_ptr, double maxevents = 0);

  // sowohl 1d als auch 2d aus spruengen addieren

  void addjumps(TDwell_1d* dwell_1d = NULL, double maxevents = 0);
  
  //Tobias:Procedur zum Schreiben eines Log-Files, das die Sprungl�ngen sortiert nach Abfolge bei potentialgesteuerten Zeitreihen enth�lt
  void Hinkley_Log(TDwell_1d* dwell_1d = NULL, double maxevents = 0);
  
  //Tobias: Procedur zum berechnen der Dependency-Plots nach Megleby & Song
  void calc_dependency();

  TDwell_2d& operator+=(TDwell_2d& B);

  virtual void GNUplot(bool dosqrt = true);

  virtual double min(); // minimale anzahl von events in einem bin
  virtual double max(); // maximale anzahl von events in einem bin

  double  events();     // gesamtanzahl der events

  double lnlikelihood(TDwell_2d &simul);  // bildet loglikelihood aus messwerten und simulierten werten!

  double checkrange(TDwell_2d &simul);    // wieviel prozent liegen im simulations bereich?

  TDwell_2d& operator-(TDwell_2d &B); // a.operator-(B) -> A - B

  // abspeichern ascii
  friend std::ostream& operator<<(std::ostream& out, TDwell_2d& Dwell);
  // einlesen ascii
  friend std::istream& operator>>(std::istream& in, TDwell_2d& Dwell);   

  virtual void store(std::ostream& out);
  virtual void load(std::istream& in);   // fuer jobfile 
};

#endif /* _DWELL _H_ */
