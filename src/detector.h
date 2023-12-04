/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#ifndef	_DETECTOR_H_
#define	_DETECTOR_H_

//#include "glib.h"
//#include "gtk/gtk.h"
//#include "gdk/gdk.h"

#include "daten.h"

/******************************************************************************
 * Tsingledetector: Detector 4. Ordung, speichert das zu testende Niveau ab      
 * Methoden: -Init: Erzeugt Object}
 *                  Parameter des Init constructors:                          
 *                  open:             zu testender Kanal                    
 *                  extra_open:       zu testender Sub - Kanal              
 *                  m0:               aktueller Mittelwert                  
 *                  alastmin:         letztes Minimum                       
 *                  aparent:          Datenobject                           
 *           -Zero: L�scht Testwert                                              
 *           -Cumsum: Cumulative Summe f�r Filter h�here Ordungen                
 ******************************************************************************/

const int filter_ord = 4;

struct Tsingledetector
{   
  TDaten   *parent;
  double  h[9];
  double  greatest_h;   // Nur f�r SHD -> Diagnosis_mode
  double  p;
  double  m;
  int   n_open;
  int   n_extra_open;
  int   lastmin;
  int   lastmax;      // Nur f�r SND -Diagnosis_mode
  double  lambda, sigma;
  bool jump;

  Tsingledetector (int open, int extra_open, double m0, int alastmin, TDaten *aparent);
  void zero();
  void cumsum();
};

/******************************************************************************
 * Thohdetector: HOHD Nach Draber/Schultz                                                   
 ******************************************************************************/

struct Thohdetector
{
  TDaten   *parent;
  Tsingledetector *detector;
  std::vector<Tsingledetector> detectors;
  int   n_open, n_extra_open;
  double  m;
  int   new_model;
  double  SNR;
  bool jump;
    
  Thohdetector(TDaten  *aparent,
               int  &i,
               int  &an_open, int &an_extra_open,
               double ai_null, double ai_channel, double ai_extra_channel,
               int  aSamples);

  virtual ~Thohdetector(){};

  virtual void CalcHinkley(double e, int &i);
  virtual void CheckJump(int &i); //Ueberschreitet irgendein Detektor Schranke
  virtual void Create_N(int &i);

private:
  int Vorzeichen(double x);
  void CalcH(Tsingledetector &item, int i, double e);
  void CheckJ(Tsingledetector &item); 
};

/*******************************************************************************
 * tshdetector: Sublevel-Detector nach Draber/Schultz                            
 *******************************************************************************/

struct Tshdetector : Thohdetector
{
  bool diagnosis_mode;
  int   diagnosis_counter;
  double  epsilon;  // Es gibt nur eine Schranke f�r alle Niveaus

  Tshdetector(TDaten  *aparent,
              int  &i,
              int  &an_open, int &an_extra_open,
              double ai_null, double ai_channel, double ai_extra_channel,
              int  aSamples);

  virtual void Create_N(int &i);
  virtual void CalcHinkley(double e, int &i);
  virtual void CheckJump(int &i);
private:
  void CalcH(Tsingledetector &item, int i, double e);
  void CheckJ(Tsingledetector &item);
  void Greatest_H(Tsingledetector &item, int &i, bool &positiv_h);
  void Decide(Tsingledetector &item, double hmax);
};

/******************************************************************************
 * Tdetector: neuer Hinkley Detector, erkennt Sublevel und kann auf h�here      
 *            Ordungen eingestellt werden                                       
 * Methoden: -Init: Erzeugt Object                                              
 *                    Parameter des Init constructors:                          
 *                      aparent:          Elternobject                          
 *                      var i:            momentane Position                    
 *                      var jump_pos:     Position des letzten Sprunges         
 *                      var n_open:       aktueller Kanal                       
 *                      var n_extra_open: aktueller Sub - Kanal                 
 *           -Vorzeichen: Berechnet Vorzeichen eines Real Parameters            
 *           -Create_Neighbor: Erzeugt Collection (Liste) benachbarter Niveaus  
 ******************************************************************************/

struct Tdetector
{
TDaten   *parent;
int   new_model;
double  m,sigma;
Tsingledetector *detector;
std::vector<Tsingledetector> detectors;
int   n_open, n_extra_open;
int   jump_pos;

  Tdetector(TDaten *aparent,
            int &i,
            int &ajump_pos,
            int &an_open, int &an_extra_open);

  int Vorzeichen(double x);
  void Create_Neighbors();
private:
  void calc_sigma(Tsingledetector &item, int &t_res);
  void FindMinSigma(Tsingledetector &item);
  void FirstJump(Tsingledetector &item);
  void CheckTresTooLong(Tsingledetector &item, int &t_res);
  void CalcHinkley(Tsingledetector &item, double e, int &i);
};

#endif	/* _DETECTOR_H_ */
