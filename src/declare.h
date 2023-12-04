/*******************************************************************
 *  Kiel-Patch 
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest    
 *  CAU KIEL
 *******************************************************************/  
#ifndef	_DECLARE_H_
#define	_DECLARE_H_

#include <fstream> //.h weg
#include <vector>  //.h weg


//#include "glib.h"
//#include "gtk/gtk.h"
//#include "gdk/gdk.h"
#include "mpi.h"	//Tobias 2013 multiprocessor
#include "matrix.h"


const short fsPathName = 1024; //FILENAME_MAX;  

const char Progname[]="2D-Fit Kiel/Erlangen";
const char Proginfo[]="Version 4.00";

//BEGIN KARSTEN1
const char Proginfo1[]="Zentrum f�r Biochemie und";
const char Proginfo2[]="Molekularbiologie der CAU Kiel";
const char Proginfo3[]="Institut f�r Physiologie und";
const char Proginfo4[]="Pathophysiologie der FAU Erlangen";
//END KARSTEN1

const unsigned short MaxTextLen = 255;         // Textl�nge in Eingabefeldern und char[]

const short n_channels_max = 13;      // 13 Kan�le
const short n_extra_channels_max = 2; // 2 Sublevel
const short n_bins_max = 2000;        // im Dwell-Time-Histogram bis t=2000+Ts
const short quenched_bins_max = n_bins_max / 8; // Hoffentlich nie mehr als 250 Datenpunkte im Histogramm
 
// Define colors (0..255) fuer GdkColor * 256
const unsigned short BACKCOLOR_R = 0;
const unsigned short BACKCOLOR_G = 0;
const unsigned short BACKCOLOR_B = 0;

const short WINBORDERCOLOR_R = 0;
const short WINBORDERCOLOR_G = 191;
const short WINBORDERCOLOR_B = 191;

const short DOTCOLOR_R = 255;
const short DOTCOLOR_G = 255;
const short DOTCOLOR_B = 255;

const unsigned short CURRENTCOLOR_R = 255;
const unsigned short CURRENTCOLOR_G = 0;
const unsigned short CURRENTCOLOR_B = 0;

const unsigned short MARKCOLOR_R = DOTCOLOR_R;
const unsigned short MARKCOLOR_G = DOTCOLOR_G;
const unsigned short MARKCOLOR_B = DOTCOLOR_B;

const unsigned short HISFITCOLOR_R = WINBORDERCOLOR_R;
const unsigned short HISFITCOLOR_G = WINBORDERCOLOR_G;
const unsigned short HISFITCOLOR_B = WINBORDERCOLOR_B;

const unsigned short HOHD_default_Color_R = 255;
const unsigned short HOHD_default_Color_G = 191;
const unsigned short HOHD_default_Color_B = 0;

const unsigned short SHD_default_Color_R = 0;
const unsigned short SHD_default_Color_G = 0;
const unsigned short SHD_default_Color_B = 255;

const unsigned short HOSHD_default_Color_R = 0;
const unsigned short HOSHD_default_Color_G = 255;
const unsigned short HOSHD_default_Color_B = 0;

const short R1 = 0;   // Farbe Gruen f�r Punkte mit hoher Gewichtung
const short G1 = 255;
const short B1 = 0;

const short R2 = 255; // Farbe Rosa f�r Punkte mit kleiner Gewichtung
const short G2 = 0;
const short B2 = 255;

const short R3 = 0;   // Farbe Blau f�r Punkte, die nicht gefittet werden -> wg. Anstieg
const short G3 = 0;
const short B3 = 255;

const short MaxLen=20;                // Statuszeile?

// Verstaerkungsfaktoren fuer Dagan und EPC7
const double VDagan=1.0;   //ganz neu:1.0 f�r pclamp neu:1.257 Dagan: Maike; Silke/Roland(= 2.2897)  
const double VEPC7=2.5;       

// struktur fuer Bitmustertabelle:
// hier soll jede 'neue' 4bit markierung gespeichert werden
struct TMarkBits
{
  unsigned short    mark;    // welcher marker? nur 4bit benoetigt, 
                     // 0-F; 0: kein marker 
  int   startpos; // an welcher pos ist der marker neu
 //Tobias begin
  bool    chosen;   //fuer die Auswahl einer Markierung in der Zeitreihe 
 //Tobias end; 
  TMarkBits *next;   // einfach verkettete liste

  TMarkBits(unsigned short aMark, int aStartPos, bool achosen); //achosen von Tobias
};

// typ der mittelung (nur anzeige; gleitend nicht implementiert)
enum Mittelwertbildung {keine, normal, gleitend};
std::ostream& operator<<(std::ostream& out, Mittelwertbildung mit);

enum Verstaerk {Dagan, EPC7};
std::ostream& operator<<(std::ostream& out, Verstaerk ver);

// fit-parameter, der gesetzt/ver�ndert werden soll

enum commandtype {cnone,cbaseline, cccurrent, cucurrent,
                  cgreatchannels, csmallchannels, csigma};

// Struktur, die in TSettings auftaucht
struct TMSet
{
  short            StackedWindows;
  Mittelwertbildung AvarageFilter;
  short            LengthOfAvarageFilter;
};

// F�r Settings 
struct TSettings
{
  bool  Zoom;

  Verstaerk Ampl;
  double   VFactor;
  double   Gain;
  double   Samplininterval;

  bool  PositiveMembranPotential;
  int    SamplesToSkip;
  int    SamplesToProcess;

  unsigned short   CurrentFrom;
  unsigned short   CurrentTo;
  TMSet     Zoomed;
  TMSet     UnZoomed;
  bool  enlarge_current;

  int    multi_samples;
  int    x, y, w, h; // nur um fensterpos laden zu koennen

// rueckgabewert abhaengig von zoom (nicht redundant, wie bei day!)
  short     StackedWindows();
  Mittelwertbildung     AvarageFilter();
  short    LengthOfAvarageFilter();

// hier stehen randbedingungen als klassenkonstanten
  static const short LengthOfAverageFilterMin=1;
  static const short LengthOfAverageFilterMax=3000;

  TSettings();
  void load(std::istream& in);    // zum laden und 
  void store(std::ostream& out);  // abspeichern ins jobfile

  // abspeichern ascii
  friend std::ostream& operator<<(std::ostream& out, TSettings& Settings);
  // einlesen ascii
  friend std::istream& operator>>(std::istream& in, TSettings& Settings);   
};
std::ostream& operator<<(std::ostream& out, TSettings& Settings);
std::istream& operator>>(std::istream& in, TSettings& Settings);   

// Amplituden Histogram und Fit-Parameter
const short AnzLevel = (n_channels_max+1) * (n_extra_channels_max+1);  // Zahl der Level
 
typedef int Tampl_histo_type[65536]; 


struct Tampl_fit_type
{
  double i_null;
  double i_channel;
  double i_extra_channel;
  double sigma; // in Wandlereinheiten

  double a[AnzLevel+1]; //Amplituden der Gau�kurven [1..AnzLevel]
  short  n_channels;
  short  n_extra_channels;
};

// was wird im fenster dargestellt? z.B. zeitreihe oder ampl.histo.
enum show_graph {nothing, growing_ampl_histo, ampl_histo, time_series, dwell_time, start_ampl_histo, multi_time, dwell_2d, dwell_diff};

// was ist in datapix gespeichert? unknown immer wenn aenderung noetig
enum show_data {unknown, single, multi};

// hoehe der statuszeile im fenster
const int Status_Height = 20;

// Erkannter Sprung, speichert Pos. und absolute Sprungh�hen vor/nach Sprung
struct TSprung 
{
  double sprung_i_null;
  double sprung_i_channel;
  double sprung_i_extra_channel;  // wegen altem design, muss in jedem sprung verfuegbar sein!
	
  unsigned short  sprung_C_1;
  unsigned short  sprung_U_1;
  unsigned short  sprung_C_2;
  unsigned short  sprung_U_2;              // so kann ich die nivaus wieder zurueckbekommen!

  int  sprung_position;
  
public:
  unsigned short  C_1();
  unsigned short  U_1();
  unsigned short  C_2();
  unsigned short  U_2();
  int  position();

  void    level_1(unsigned short aC_1, unsigned short aU_1);
  void    level_2(unsigned short aC_2, unsigned short aU_2);
  void    position(int aposition);
  void    setnvs(double anull, double achannel, double aextra_channel); 

  double  level_1();
  double  level_2();
  
  // fuer plot
  static short i_oldwindow;
  static short x_pixel_old;
  static short y_pixel_old;
  static short y_pixel_max; // von wo bis wo sind
  static short y_pixel_min; // spruenge schon gezeichnet
  // reset setzt die klassenvariablen zur�ck und zeichnet das niveau nach
  // letztem sprung bis zum ende (aufruf sollte also nach 'plotSpruenge' sein)
  static void reset(bool draw=true);
 	
  // fuer calcStD
  static int p;
  static short j;
  static int pos1, pos2;
  static bool YetStarted;
  static double SD;
  static double calcSD(double level, std::vector<TSprung> &spruenge);
 	
  TSprung();
  void load(std::istream& in);
  void store(std::ostream& out);
  void plot(short strom_min, short strom_max);

  void calcStD(double level); // fuer einen Sprung wird von calcSD aufgerufen
};


// das Gequetschetes Dwell-Time-Histogram wird in Liste (PCollection) gespeichert.
struct TDwell
{
  double  x;
  double  y;
  
  TDwell(double ax, double ay);
  void load(std::istream& in);
  void store(std::ostream& out);
};

// Benachbarte Niveaus f�r Hinkley Detektor (veraltet, wird sobald wie m�glich vegrationalisiert)
struct Tneighbor_type
{
  short  n_neighbors;
  short  n_open[AnzLevel+1];
  short  n_extra_open[AnzLevel+1];
  double m[AnzLevel+1]; // aktuelles Niveau
  //Kenngroessen die der DHD braucht bzw. bildet:
  double p[AnzLevel+1]; // halbe Sprungh�he
  double h[AnzLevel+1]; // Testwert
  double greatest_h[AnzLevel+1]; // gr��ter Testwert
  double lambda[AnzLevel+1];   
  int  lastmin[AnzLevel+1]; // letztes Minimum( Testwert = 0)
  int  lastmax[AnzLevel+1]; // Maximum des Testwertes
};

enum detector_type {_SHD, _HOHD, _HOSHD, _SHD_AH};

//Ergebnisse der Sprungerkennung
struct Tresult
{
  short   n_bins;
  short   candidate_restriction;
  short   filter_ord;
  matrix<int>   dwell_time_histo_array;
  int   detected_level_histo[AnzLevel+1];
  int   transitions[AnzLevel+1][AnzLevel+1];
  detector_type  detector;
  
  //void load(istream& in);
  //void store(ostream& out);
};

typedef double g_type[9];

struct Tniveau
{
  double  n;
  short   channel;
  short   extra_channel;

  Tniveau(double f, short i, short j);
  void load(std::istream& in);
  void store(std::ostream& out);
};

//Farben f�r Color And Order Dialog
struct TColorAndOrderRec
{
  short   DrawingOrder;
  bool HOHD_Showflag;
  unsigned short   HOHD_R;
  unsigned short   HOHD_G;
  unsigned short   HOHD_B;
  bool SHD_Showflag;
  unsigned short   SHD_R;
  unsigned short   SHD_G;
  unsigned short   SHD_B;
  bool HOSHD_Showflag;
  unsigned short   HOSHD_R;
  unsigned short   HOSHD_G;
  unsigned short   HOSHD_B;
  
  TColorAndOrderRec();
};

// F�r History-Dialog
struct THistorytext
{
  int  detected_jumps;
  double sigma_of_channel;
  short  level;
  short  sublevel;
  detector_type  detector;
};

// Fuer Dialog zum Festlegen des Grow-Factors bzw. der Zahl des Samples
// und zur Wahl des darzustellenden Dwell-Time-Histograms

struct TSDwellTime
{
  short   LevelNo;
  short   SublevelNo;
  short   Samples;
  double  GrowFactor;
  bool SamplesB;
  bool open;
  bool openclosed;
  bool cross;
  bool ExpAvarage;
};

// f�r Exponentialfit: Startwerte (tau's) und Zahl der Zeitkonstanten
const short NMaxTimeConstants = 12;

const short TargetFitMaxChannels = NMaxTimeConstants; // max 12 Kanl�le, 13 Niveaus

typedef double TTau[NMaxTimeConstants+1];

typedef double TAmpli[NMaxTimeConstants+1];

typedef double Tsig[quenched_bins_max+1];

typedef double TTargetWeight[TargetFitMaxChannels+1][quenched_bins_max+1];

struct TDwellTime
{
  short   NTimeConstants;
  TTau     tau;
  TAmpli   Ampli;
  double  Error;    // Fehler f�r Exp. Fit
  double  ErrorCOG; // Fehler f�r COG Fit
  double  ErrorSCO; // Fehler f�r SCO Fit
  double  ErrorFiveState; //Fehler Five State Modell
  int   iter;
  short   Points;
  double  startx;
  bool IsTargetFit;
};

const short closing = 0;
const short opening = 1;
const short openclosing = 2;

typedef TDwellTime TLevel[n_channels_max+1+1][n_extra_channels_max+1+1][3];


// F�r Targetfit: �bergangsraten bei Start des Fits (f�r ExpFitDialog)
struct TRates
{
  double CO;
  double OC;
  double OG;
  double GO;
};

// F�r Targetfit routine func
struct Tpnt
{
  double x;
  double y;
  double s;
};

typedef double Tmama[29][29];
typedef double Tcmatrix[TargetFitMaxChannels+1][TargetFitMaxChannels+1];
typedef Tcmatrix  Tmatrix;

typedef Tpnt Tlevelarray[quenched_bins_max+1];

struct  TLeveldaten //in PASCAL TObject
{
  Tlevelarray    d[TargetFitMaxChannels+1];
// in PASCAL
// Tlevelarray*    d[TargetFitMaxChannels+1];
// TLeveldaten();
// virtual ~TLeveldaten();
// aber konstruktor und destruktor nur fuer speicher 
// des arrays zustaendig 
};

typedef double histotype[quenched_bins_max+1];
// nicht mit vector template verwechseln!
typedef double Vector[TargetFitMaxChannels+1]; //6 Kan�le, 7 Niveaus
typedef short  row_vec[TargetFitMaxChannels+1];


enum TModell{M_COG, M_SCO, M_CGO, M_5State, M_Koop};

// F�r algemeinen TargetFit

const short max_dim_big_matrix = 132;  // 4+9+12+10+1; was wei� ich, wie gro� die Matrix sein mu�!
const short max_dim_unter_matrix = 13; // keine Ahnung was hier stehen mu�
const short b_nr_max = max_dim_big_matrix;
const short n_states_max = 5;

struct TRate_matrix
{
  TRate_matrix();

  double  r[n_states_max][n_states_max];
  bool open[n_states_max];
};

typedef double TRate_m0[n_states_max][n_states_max];
typedef double TRate_m1[n_states_max+1][n_states_max+1];

typedef double TVectorAlg[max_dim_big_matrix+1];
typedef TVectorAlg TMatrixAlg[max_dim_big_matrix+1];

typedef short  Tb[n_states_max+1][b_nr_max+1];
typedef double TSingleSteady[n_states_max+1];
typedef double TMultiSteady[max_dim_big_matrix+1];

typedef double TVectorAlg2[max_dim_big_matrix+2];

typedef TVectorAlg2 TSxS_transpo[max_dim_big_matrix+2];

inline unsigned short  
TSprung::C_1()
{
  return sprung_C_1;
}

inline unsigned short  
TSprung::U_1()
{
  return sprung_U_1; 
}

inline unsigned short  
TSprung::C_2()
{
  return sprung_C_2;
}

inline unsigned short  
TSprung::U_2()
{
  return sprung_U_2;
}

inline int  
TSprung::position()
{
  return sprung_position;
}

inline void    
TSprung::level_1(unsigned short aC_1, unsigned short aU_1)
{ 
  sprung_C_1 = aC_1;
  sprung_U_1 = aU_1;
}

inline void    
TSprung::level_2(unsigned short aC_2, unsigned short aU_2)
{ 
  sprung_C_2 = aC_2;
  sprung_U_2 = aU_2;
}

inline void    
TSprung::position(int aposition)
{ 
  sprung_position = aposition;
}

inline void    
TSprung::setnvs(double anull, double achannel, double aextra_channel)
{ 
  sprung_i_null = anull;
  sprung_i_channel = achannel;
  sprung_i_extra_channel = aextra_channel;
}

inline double  
TSprung::level_1()
{
  return sprung_i_null +
         sprung_i_channel * sprung_C_1 +
         sprung_i_extra_channel * sprung_U_1;
}

inline double  
TSprung::level_2()
{
  return sprung_i_null +
         sprung_i_channel * sprung_C_2 +
         sprung_i_extra_channel * sprung_U_2;
}


#endif /* _DECLARE_H_ */
