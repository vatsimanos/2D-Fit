/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/
#ifndef	_FIT_H_
#define	_FIT_H_


#include <ext/slist>  //tobias



//#include "glib.h"

#include "declare.h"
#include "matu.h"

// Einiges f�r Fit im Amplitudenhistogram mit Simplex/Gau� (oder Exp.Fit)

const short ndatap = 4096;                // 4096 Datenpunke maximal

const short fit_dim = NMaxTimeConstants;  // maximale Dimension bei Amplitudenfit: Nullstrom + Kanalstrom + Sublevelstrom + Sigma = 4
                                           // bei Exponentialfit: NMaxTimeConstants = 12
                                           // bei Alg. Modellfit: (N_states_max-1)*2 = Zahl der Ratenkonstanten

const short map = AnzLevel; // Zahl der Level

typedef double RealArrayNData[ndatap+1];
typedef double RealArrayMA[map+1];
typedef short  IntegerArrayMFIT[map+1];
typedef double RealArrayMAbyMA[map+1][map+1];
typedef double RealArrayNPbyMA[map+1][2]; //[1..map,1..1]

const short np = (n_channels_max+1) * (n_extra_channels_max+1);
const short mp = 1;

typedef RealArrayMAbyMA RealArrayNPbyNP;
typedef double RealArrayNPbyMP[np+1][mp+1];
typedef IntegerArrayMFIT IntegerArrayNP;
typedef double glmpnp[fit_dim+2][fit_dim+1];
typedef double glmp[fit_dim+2];
typedef double glnp[fit_dim+1];
typedef short  glnpint[fit_dim+1];

/*
 * Tlfit: Least Square Object, basierend auf N.R.                         
 * Methoden: 
 *           -init: Speicher belengen                                     
 *           -done: Speicher freigeben                                    
 *           -lfit: Least Square Routine aus N.R.                         
 *           -funcs: Ist eine Dummy-Funktion (virtual), die in anderen    
 *                   Objecten z.B durch e-Fkt oder Gauss-Fkt ersetzt wird 
 *           -func: Fehlerfunktion, Quadratische Abweichung durch Aufruf  
 *                  von lfit                                              
 *           -gaussj: Gauss-Jordan zur L�sung linearer Gleichungen        
 * Variablen:
 *           -singularMatrix, timeconstantSero: Fehlerflags f�r Abruch des
 *                                              Algorithmus               
 *           -ap: Fitparameter                                            
 *           -NumPoints: Anzahl der Datenpunkte                           
 */

struct Tlfit               // least Squares
{
  bool singularMatrix;
  bool timeconstantSero;

  RealArrayNData   ax;     // x-Datenwerte
  RealArrayNData   ay;     // y-Datenwerte
  RealArrayNData   asig;   // Standartabweichungen
  
  RealArrayMA      aa;     // Amplituden der Gausskurven (Ergebnis)
  IntegerArrayMFIT alista; 
  RealArrayMAbyMA  acovar; 
  
  double  achisq;
  short   ama;            // Zahl der Basisvektoren (Gausskurven)
  short   amfit;          // Basisvektoren, deren Ampl. gefitet wird, meist = ama
  short   NumPoints;      // Anzahl der Datenwerte
  // HWindow: HWnd;        // Window Handle (f�r Ausgabe von Fehlermeldungen)
  
  glnp     ap;             // enth�lt: i_null,i_channel,i_extra_channel und sigma
                           // f�r die Berechnung der Gaussfunktion und als Startwert f�r Simplex

  Tlfit();                 // init

  virtual ~Tlfit(){};      // virtueller destruktor, da virtuelle fkt vorhanden

  // Least Square Routine
  void lfit(RealArrayNData    x, 
            RealArrayNData    y,
            RealArrayNData    sig,             
            short            ndata, 
            RealArrayMA       a,                           
            short            ma,
            IntegerArrayMFIT  lista,
            short            mfit,
            RealArrayMAbyMA   covar,
            double           &chisq);

  // mehrere Gaussfunktionen
  virtual void funcs(double z, RealArrayMA afunc); 

  // Fehlerfunkion,
  // Quadratische Abweichung der Gaussumme von den Messdaten
  virtual double func(glnp pr);                
                       
  // Gauss-Jordan zur L�sung lineare Gl.-> wird von lfit aufgerufen
  void gaussj(RealArrayNPbyNP a,   
              short          n,
              RealArrayNPbyMP b,      
              short          m);

  void covsrt(RealArrayMAbyMA  covar,
              short           ma,
              IntegerArrayMFIT lista,
              short           mfit);
};

/*
 * TAmoeba: Downhill Simplex Object, basierend auf N.R.
 * Methoden:
 *           -amoeba: Downhill Simplex                                    
 *           -fit: kombinierter Fit: Simplex f�r Exponenten, lfit f�r     
 *                 Amplituden der e- /Gauss- Fkt.                         
 */

struct Tamoeba : public Tlfit 
{
  bool ende;
  glmpnp   bp;      // dim+1 Startwert der Dimension dim
  glmp     by;      // Funktionswerte an der Stelle der Startwerte
  short   bndim;   // Dimension dim eines (Start-)wertes
  double  bftol;   // Abbruchgenauigkeit
  short   biter;   // Zahl der ben�tigten Iterationen (Ergebnis)

  Tamoeba();        // init

// Downhill Simplex-Routine
private:
  void startwerte();
  void funktionswerte(glnp cp);
  void ergebnis_sichern();
  bool abbruch(double ftol);
public:
  void amoeba(glmpnp  p,   
              glmp    y,
              short  ndim,
              double ftol,
              short  &iter);

  // kombinierter Fit: Gauss-> Amplituden Simplex-> Kanalstr�me u. Abweichung
  double fit(double ftol, short &iter);    
};

/*
 * TGaussFit: Fit mit Summe von Gaussfkt.                                
 * Methoden: 
 *            -funcs: Funktionswerte jeder Gaussfkt, erzetzt Dummy ftk. von
 *                    oben                                                 
 *            -summe_gauss_glocken: Summe aller Gaussfktionen, wird i.a nur
 *                                  zum Darstellen der Gausskurven auf     
 *                                  Schirm verwendet                       
 *            -determine_base: Bestimmung der Basisvektoren, gleiche       
 *                             Niveaus werden weggelassen                   
 */
struct TGaussFit : public Tamoeba
{ 
  class __gnu_cxx::slist<Tniveau> niveaus;        // Liste mit Kanalst�men und deren Kanal-Indices, f�r die Berechung der Gausskurven; Tobias oldlib:slist <Tniveau> niveaus;
  short  n_channels;            // Anzahl der Kan�le
  short  n_extra_channels;      // und Sublevel

  TGaussFit(Tampl_fit_type &af); // Speicher belegen
  void done(Tampl_fit_type &af);
  

// mehrere Gaussfunktionen
private:
  void calculate(Tniveau niveau, 
                 double z, 
                 RealArrayMA afunc,
                 short  &k, 
                 double &m, 
                 double &b);
  bool IsAlreadyThere(double akt_niveau);

public:
  void funcs(double z, RealArrayMA afunc);

// Summe aller Gaussfunktionen
  double summe_gauss_glocken(RealArrayMA  afunc);        
// eine gaussfunktion
  double gauss_glocke(short niveau, RealArrayMA  afunc);

// Bestimmung der Basisvektoren in Liste, R�ckgabe: Anz. der Vektoren
// gleiche Niveaus werden weggelassen! ->"Basis" d.h. l.u.
  short  determine_base();                      
};


/* TExpFit: Fit mit Summe von e-Funtionen
 * Methoden: -funcs: Funtionswerte jeder e-Fkt., ersetzt Dummy fkt. von 
 *                   oben                                                 
 *           -summe_expo: Summe aller e-Fkt., wird i.a. nur zum Darstellen
 *                        der e-Fkt. auf Schirm verwendet                 
 */

struct TExpFit : public Tamoeba
{
  TExpFit(TTau &tau);
  void done(TTau &tau);
  ~TExpFit();
  void funcs(double z, RealArrayMA afunc); // mehrere Exponentialfunktionen
  double summe_expo(double x);  // Summe aller Exp. funktionen
};

/* TTargetFit: Rikards Target Fit                                         
 * Methoden: -func: Die func-Methode von oben ersetzen, da lfit hier      
 *                  nichts zu suchen hat. Die Amplituden werden aus Modell
 *                  bestimmt und nicht mit lfit.                          
 *           -Rates_to_TimeConstants: aus Raten Amplituden und Zeitkonst  
 *                                    bestimmen                           
 *           -Factor: Vorfactor vor allen e-Fkt.                          
 *           -calc_matrix: Gesamtzustandsmatrix aufstellen                
 *           -calc_row_dim: Summe f�r Indizierung der Level               
 * Variablen: -n_channels: Zahl der Kan�le                                
 *            -tres: Aufl�sungszeit                                       
 *            -samplelaenge: Zahl der Samples in Zeitreihe                
 *            -data: Histogram-Daten                                      
 *            -events: Zahl der Events, wird nicht ben�tigt, da sie auch  
 *                     schon aus der Theoie folgen, wenn man die Zahl der 
 *                     Samples kennt, nur zum Vergleich                   
 */

struct TTargetFit : public Tamoeba
{
  short  n_channels;
  double tres;
  double samplelaenge;
  TLeveldaten    data;
  int  events;

  TTargetFit(TRates &k, short an_channels, double atres, double asamplelaenge, int aevents);
  TTargetFit(); // weil konstruktor vom TTargetFitAlg unkompatibel
  void done(TRates &k);
  ~TTargetFit();
  void Check_Sero_Rates();
  virtual double func(glnp pr);
  double Fehlersumme(Vector lambda[], Vector ampli[]);
  void  rates_To_TimeConstants(Vector lambda[], Vector ampli[]);
  double Factor(short level, Vector lambda[], Vector ampli[]);
  void calc_matrix(double &k12, double &k21, double &k23, double &k32, short ch, Tmama mama); // Gesamtmatrix aufstellen
  void calc_c_matrix(short k, Tcmatrix cmatrix, Tmama mama); // Reduzierte Matrix, deren EW und EV bestimmt werden sollen
  void calc_row_dim();
  void show_big_matrix(Tmama mama, short dim);
};


/* TUntermatrix: Dimension und Position einer Untermatrix gesamten Matrix
 */

struct TUntermatrix
{
  short  start_index;
  short  stop_index;
  short  dimension;

  TUntermatrix(short astart_index, short astop_index, short adimension);
};



#endif	/* _FIT_H_ */
