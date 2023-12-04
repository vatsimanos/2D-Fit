/*******************************************************************
 *  Kiel-Patch
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest  
 *  CAU KIEL
 *******************************************************************/

#ifndef	_DATEN_H_
#define	_DATEN_H_

#include <vector> //.h weg
#include <random>

//#include "glib.h"
//#include "gtk/gtk.h"
//#include "gdk/gdk.h"

//#include "declare.h"
#include "gtkapp.h"
#include "fit.h"
#include "dwell.h"
#include "SimSettings.h"
#include "matrix2.h"




/* TDaten: Grundlegendes Datenobjekt, ist verantwortlich f�r:                    
 *         keine Vererbung, da nur ein Objekt in day!
 * Methoden:  
 *            -get_data_from_file: Laden der Zeitreihe                           
 *            -min_max: (Dummy d.h. virtual)     
 *            -free: (RAM-) Speicher wieder freigeben                            
 *            -load: Initialisieren von Variablen durch laden eines job-Files  //Konsrtuktor  
 *            -store: Abspeichern von Variabeln auf Platte (job-File)            
 *            -markeddatafit: ...                         (toby)                 
 * TPaintDaten:
 *            -plot: Zeichnet Object (im Speicher)                               
 *            -paint: Kopiert Geplottetes Objekt auf Bildschirm                  
 * THistoDaten:
 *            -ampl_grow: Zeichnet wachsende Amplituden-Histogramme              
 *            -copy_stromsig: ebenfalls, ruft ampl_grow auf                      
 *            -plot_strom_height: kann Linien im Histogram zeichen               
 *            -ampl_min_max: bestimmt min und max aus Histogram der Zeitreihe    
 *            -min_max: ruft ampl_min_max u. histogarm auf, um das Starthistogram
 *                      u. strom_min + strom_max zu berechnen                    
 *            -histogram: begrechnet Histogram aus allen Samples ohne filter                
 *            -paint: Zeichnet wachsende Histogramme durch Aufruf von            
 *                    copy_stomsig, eigentliche Histogramme -> sp�ter (HistoFit) 
 *               hier muss auch entscheden werden, was gezeichnet wird!                                           
 *               oder painthisto
 * TZeitDaten:
 *            -str_filter: Berechnet im Speicher Bitmap von Ausschnitt aus      
 *                          Zeitreihe                                           
 *            
 *            -plot_init: initialisiert Plot-Routine                            
 * TFitDaten:
 *            -plot_levels: Zeichnet Niveaus, ruft f�r jedes Niveau Plot_Strom auf
 *            -plot_strom: Zeichnet ein Niveau
 *            -plot_text: Zeigt Parameter in Statuszeile
 *            -updown: ver�ndert Lage eines Niveaus
 * THistoFitDaten: 
 * Interaktiver Fit im Amplituden-Histogram
 *            -show_ampl_histo: Plotet Histogram im Speicher
 *            -fehlertext: gibt Abweichung Histogram-Gau� aus
 *            -calc-fit: Berechet Amplituden der Gau�-Kurven
 *            -call_simplex: Fit mit Simplex (statt interaktiv)
 *            -simplex: der Simplex Algorithmus, ruft u.a. NR-Routine Amoeba auf
 *            -plot_levels_and_fit: Zeichnet Gau�kurve und Niveaus ein
 * TFilterDaten:
 * Verschiedene Hinkley Filter, die Auf die Zeireihe losgelassen
 * werden k�nnen
 *            -InitDwellTimeHisto: Dwell-Time-Histogramme l�schen (falls vor-
 *                                 handen) und dann neu anlegen und nullsetzen
 *            -hinkley: einen der Detektoren ausw�hlen und aufrufen              
 *            -HOSHD: neuer Typ, experimentell                                   
 *            -HOHD: der Detektor h�herer Ordnung                                
 *            -SHD: Sublevel Detektor                                            
 *            -HOHD_D: Detektor h�herer Ordnung mit wahlfreier Schranke          
 *            -flush_Sprung: Sprung in Liste aufnehmen, f�r alle Detektoren      
 *            -PlotSpruenge: Spr�nge zeichnen, wird von Paint aufgerufen         
 *            -SaveDwellTimeHisto: Dwell-Time-Histogramme Speichern              
 *                      Speichert Dwell-Time-Histogramme in Datei. Es gibt Histo- 
 *                      gramme die open-/und close Zust�nde zusammen enthalten 
 *                      (.hxx-Files ) und welche, die dies Zust�nde getrennt 
 *                      speichern (.dxx-Files ). Die Dateien sind nach 
 *                      aufsteigender Stromreihenfolge numeriert. Zuerst kommen 
 *                      immer ein open Histogramm (gerade Zahlen) und dann ein 
 *                      closed (ungerade Zahlen).              
 *            -SaveDetectedLevelHisto: Zahl der Spr�nge in ein Niveau speichern
 *                      Speichert Histogram �ber die Aufenthaltszeiten in     
 *                      einem Niveau. Ein f�r Allemal:   
 *                      !!!!!!!! DIE Z�HLUNG DER NIVEAUS !!!!!!!      
 *                      1. das nullte Hauptniveau                                                      
 *                      2. alle Sublevels des Nullten Hauptlevels 
 *                                        von 0 - n_extra_Channels             
 *                      3. das erste Hauptlevel                                                        
 *                      4. Sublevel des 1. Hauptlevels                                                 
 *                      5. u.s.w                                                                       
 *                      nach der Formel: 
 *                            out = extra_open + (n_extra_channels+1)*open                 
 *                      Nivaus in aufsteigendes Stromreihenfolge                                       
 *            -CalcStD: Dummy (virtual)                                          
 *            -SaveTransitionMatrix: Zahl der �berg�nge von einem Niveau in ein  
 *                                   anderes Speichern                           
 *            -ClearDetectedLevels: L�scht DetectedLevel-Histo                   
 *            -ClearTransitions: L�scht �berg�nge                                
 *                      L�scht �bergangsmatrix, wichtig, da sonst �berg�nge immer  
 *                      aufaddiert werden bei Mehreren Filter-Durchg�ngen          
 *            -FindLambda: Sucht Schranke f�r False Alarm bei p=1/10000          
 *                      Findet die Schranke Lambda, bei der der 
 *                      HOHD 100 False Alarms feststellt. Bei 1 Mio Samples 
 *                      enspricht dies einer W-keit von 1/10000         
 *                      Ich gehe davon aus, da� die Schranke nicht gr��er als 
 *                      10 000 000 ist (Detektor 3. Ord)                             
 *            -Zeitaufl�sung: berechnet bei vorgegebener Schranke die 
 *                            Zeitaufl�sung nach Draber/Schultze im 
 *                            rauschfreien Fall
 * THistoryDaten: 
 * Statistik �ber erkannte Spr�nge und Standardabweichung der Niveaus
 *            -History: Zeigt Zahl der schon erkannten Spr�nge und die Standard
 *                      Standard-Abweichung in einem Niveau
 * TDwellTimeDaten:
 *            -DetermineMax: ermittelt Maximum im Histogramm (f�r Skalierung)
 *            -GetWeight: virtuelle Dummyfkt.                                
 *            -NoSamples: Berechnet NoSamples aus Grow-Factor                
 *            -ExpEinteil: Berechnet Grow-Factor aus NoSamples               
 *            -mittel: F�hrt Exp. Mittelung durch                            
 *            -axes: Zeichenet Achsenkreuz                                   
 *            -zoomin: vergr��ert Darstellung                                
 *            -zoomout: verkleinert Darstellung                              
 * TExpFitDaten: Expotential Fit im Dwell-Time-Histogram 
 *            -FitExp: Expotential-Fit                                           
 *            -CalcStartIndex: Beginn des Histograms beim gr��ten Wert, da die   
 *                             e-Funktion abf�llt.                               
 *            -Load/Store: Objekt aus Jobfile laden/Jobfile speichern            
 *                         Es wird nur ein Satz an Gewichtsfaktoren gespeichert, 
 *                         da sonst �berm��ig Hauptspeicher/Plattenplatz 
 *                         verschwendet w�rde.        
 *                         Absch�tzung: 1MB f�r Gewichtsfaktoren aller 
 *                         Histogramme ist zu viel.
 *            -GetxyString: p,x,y-Werte einer best. Stelle aus akt. Histogramm   
 *                          als String zur�ckgeben (f�r Combobox im Dialog)                 
 *            -ResetWeight: Gewichtsfaktoren zur�cksetzen                        
 *            -generate_gnuplot: erzeugt Gnuplot Graphik                         
 *            -exec_gnuplot: Ruft gnuplot auf und zeigt Graphik an,              
 *                           R�ckkehr mit Return.                                
 *            -getweight: gibt den aktuellen Satz an Gewichtsfaktoren zur�ck.    
 *                        Dies sind Weightup, Weightdown u. TargetWeight         
 *                        (Gewichte f�r Spr�nge nach oben, unten und oben+unten) 
 * TDwellFitDaten: Rickards Fitprogram COG-Modell                                
 *            -Fit: initialisiert Vertex und startet fit (amoeba)                 
 *            -Read_Data: �bernimmt DwellTimeHistogram aus bisheriger Auswertung 
 *            -NormGewicht: Normierung der aller Gewichte                        
 *            -Gewicht: einen bestimmten Wert st�rke wichten als andere          
 *            -amoeba: Simplex-Fit                                               
 *            -abspeichern: ?                                                    
 * TSCODwellFitDaten: Rickards Fitprogram SCO-Modell                            
 *                    �nderung zu oben:                                          
 *                   -Daten verkehrt herum einlesen (Read_Data),                
 *                    von n_channnel bis 0 herunter                              
 *                   -Zeitkonstanten verkehrt herum abspeichern (f�r Anzeige)   
 * TAlgModellDwellFitDaten: Fitprogram, das allgemeine n-Zustands-Modelle mit    
 *                          m-Kan�len fitten kann                                
 *                          Die Kanalzahl ist auf n_max_channels und die Zust�nde
 *                          auf n_max_states begrenzt                            
 *            -Test_Matrix: Testet Errechnen der b-Vektoren und Aufstellen der     
 *                          Gesamtmatrix                                           
 * TAlgModellDwellFitKoDaten: Fitprogramm mit Kooperativit�t                     
 *
 * neu:
 *            -calc ahisto: berechnet ahisto neu (nach �nderungen in settings!)
 *
 * Variablen: 
 *            -Settings: Eistellungen f�r Settings-Dialog                        
 *            -Wert: Array von werten                     
 *            -Anz_Samples: Zahl der Samples in Datei (im Arrary)                          
 *            -FName: enth�lt nur Dateinamen
 * TPaintDaten:
 *            -schwarzaufweiss: beim Drucken auf wei�en Hintergrund umschalten   
 *            -datatype: Plot noch aktuell? single oder multi?
 *
 * THistoDaten:
 *            -strom_min, strom_max: min, max der Zeitreihe                      
 *            -g_maxx, g_maxy: maximaler Zeichenbereich                          
 *            -graphtype: Histogram, Grow-Histog. oder Zeitreihe oder nichts?    
 *            -ahisto: aktuelles Histogram                            
 *            -bhisto: Starthistogram                                 
 *
 * TZeitDaten:
 *            -Samples_per_Pixel: Zahl der Werte pro Bildpunkt                  
 *            -MemPix: Backing pixmap for drawing area
 * TFitDaten:
 *            -afit: Fit-Parameter, d.h. Lage der Niveaus u.s.w
 *            -command: Legt fest, welcher Parameter sich gerade �ndert
 * THistoFitDaten: 
 *            -max_height: Maximum im Histogram
 *            -Fehler: Abweichung von Gau� - Histogram
 * TFilterDaten:
 *            -.._Spruenge: Sprunglisten, anfangs leer                           
 *            -result: Ergebnis der Sprungdetektion                             
 *            -candidate_restriction: Einschr�nkungen der �berg�nge bei mehreren 
 *                                    Niveaus
 *            -ColorAndOrder: Farbe, Reihenfolge der Darstellung der Spr�nge     
 *            -old_index: Zeitpunkt des letzen erkannten Sprunges (F�r Dwell-Time
 *                        Berechnung wichtig)                                    
 *            -ScTiRe: Scaled-Time-Resolution, wei�es oder blaues Rauschen bei   
 *                     Draber/Schultz                                            
 *            -SNR: Signal-Rausch-Verh�ltnis                                     
 *            -event_too_long: Flag - Soll vor zu langen Events gewarnt werden ? 
 * THistoryDaten: 
 *            -HText: Enth�lt Einstellungen, wird vom Dialog gef�llt
 * TDwellTimeDaten:
 *            -SDwellTime: Welches Dwell-Time Histogram soll angezeigt werden?
 *            -QuenchedDwellTime: Ein Gemitteltes DwellTime                   
 *            -maxx, maxy, xstart, ystart: f�r das Zeichen des Dwell-T-Histog.
 * neu:
 *            -showmarker: sollen Marker angezeigt werden!
 *            -multi_times:   wie oft sollen 
 *            -multi_samples: wie viele samples gezeigt werden
 *			  -Dwell_lb :     fuer log binning dwell time
 */


struct Tsingledetector; // schon mal bekanntmachen, wird unten definiert

struct TDaten
{ 
  // Tobias start
  #include"datentobias.h"
  // Tobias end	 
  
  
  TSettings Settings; 
  long int    Anz_Samples;
  unsigned short   *Wert;             // pointer auf array von werten
  unsigned short   *SimulatWert;      // Tobias: pointer auf array f�r die Zeitreihe der Simulation
  std::vector<TMarkBits> Marker;    // TMarkBitliste
  char      cdir[fsPathName];  // aktuelles verzeichnis merken
  char      FName[fsPathName]; // FileName, z.B. fuer Namensvorschl�ge 

  bool  schwarzaufweiss;   // flag ob schwarz auf weiss darstellung

  Tampl_histo_type  ahisto;
  Tampl_histo_type  bhisto;
  show_graph  graphtype;       // was wird angezeigt (MemPix)
  show_data   datatype;        // was wird gesichert (DataPix)

  int      g_maxx;
  int      g_maxy;
  int      strom_min;
  int      strom_max;
  int      strom_min_old;   // max. u. min. f�r vergr��erte
  int      strom_max_old;   // Darstellung mit enlarge_current

//BEGIN KARSTEN
  int      Maximaler_Strom;
  int      Minimaler_Strom;
//END KARSTEN

//BEGIN KARSTENCUT_PASTE_EDIT
  int      undo_mode;
  int      redo_mode;
  bool    refresh;
  int      copy_paste_Anz_Samples;
  unsigned short     *copy_paste_Wert1;             // pointer auf array von werten
//END KARSTENCUT_PASTE_EDIT

//BEGIN KARSTENNEW
  unsigned short     *crop_Wert1;             // pointer auf array von werten
  unsigned short     *cut_Wert1;             // pointer auf array von werten
  unsigned short     *paste_Wert1;             // pointer auf array von werten

  bool    zoom_50_active;	// Warteabfragen bis Befehl ausgef�hrt
  bool    refresh_active;	// "
  bool    zoom_in_active;	// "
  bool    zoom_out_active;	// "
  bool    cut_active;	// "
  bool    copy_active;	// "
  bool    paste_active;	// "
  bool    crop_active;	// "
//END KARSTENNEW

//BEGIN KARSTENMARKER
  unsigned short   *marker_Wert;             // pointer auf array von werten
//END KARSTENMARKER

//BEGIN KARSTENASCI
  bool import_asci_data;
  double AD_fromwert;
  double AD_towert;
  double Skalierungsvorschlag_min;
  double Skalierungsvorschlag_max;
  bool import_cancel;
//END KARSTENASCI

  double     Samples_per_Pixel;

 //2019 GdkPixmap   *DataPix;        // bitmap for datas
 //2019 GdkPixmap   *MemPix;         // Backing pixmap for drawing area

  Tampl_fit_type a_fit;
  commandtype    command;

  int      max_height;
  double     fehler;

  std::vector<TSprung>  HOHD_Spruenge;
  std::vector<TSprung>  SHD_Spruenge;
  std::vector<TSprung>  HOSHD_Spruenge;  // hier werden die Sprunge gesichert

  // vector von sprunglisten
  // keine Matrix, da vektoren nicht gleich lang sein muessen!
  // Zwischenspeicher fuer gewuerfelte Einzelkanaele
  std::vector< std::vector<TSprung> > SingleChannels;  

  Tresult     result;
  int      candidate_restriction;
  TColorAndOrderRec ColorAndOrder;
  int      old_index;
  double     ScTiRe;
  double     SNR;
  bool    event_too_long;
  double     maximum;
  double     neu_t_res;
  
  THistorytext   HText;

  TSDwellTime    SDwellTime;
  std::vector<TDwell> QuenchedDwellTime;
  double        maxx;              // maxx: maximaler x-Bereich,
  double        maxy;              // wichtig f�r das Zoomen
  double        xstart, ystart;
  double        zoomfactor;

  TLevel         DwellLevel;
  TTargetWeight  Weightup;
  TTargetWeight  Weightdown;
  TTargetWeight  TargetWeight;

  TRates         Rates;
  TModell        modell;

  TRate_matrix   Rate_Matrix;

  double        Kooperativity;

  bool       showmarker;        // sollen Marker angezeigt werden!
  bool       logscalehisto;     // histogramme log(1+n)?
  int         rand_pixel;        // wieviele Randpixel in darstellung 

  TDwell_1d      Dwell_1d_A;        // fuer log binning dwell time
  TDwell_2d      Dwell_2d_A;        // fuer log binning 2d - dwell time
  
  TDwell_2d			 Dwell_2d_MA[11];		// Tobias 2022 for multiple 2D histograms with different level
  int starting_level;
  int level_increment;
  int number_of_levels;

  TDwell_1d      Dwell_1d_B;	
  TDwell_2d      Dwell_2d_B;        // fuer 2. Reihe (Simulation)

  TDwell_2d      Dwell_2d_diff;     // fuer differenz

  TDwell_1d*     Dwell_1d_ptr;
  TDwell_2d*     Dwell_2d_ptr;      // mit welchem wird gearbeitet?

  TDaten();
  ~TDaten();                        // hier wird free() aufgerufen 

  void  freeDaten();                     // dynamischen Speicher freigeben Tobias:Name ge�ndert von free()
  void  reset_level();              // setzt N,C,G,S,U,K zurueck
  void  reset();	        	    // setzt einige Werte zur�ck und ruft free() auf!
  
  void  load(char *Filename);       // jobfile (*.jbf) laden und 
  void  store(char *Filename);      // speichern inkompatibel zu pascal *.job
  int get_data_from_file(char* FileName); // -1: fehler   0: OK
  
  
  
//BEGIN KARSTENMARKER
  void get_marker();
//END KARSTENMARKER
  
  void  markeddatafit();            // (toby)    
 
  void  paint(); 
  void  plot();  

  void  histogram(Tampl_histo_type a);     
  void  min_max();
  void  ampl_min_max(Tampl_histo_type a, int &min, int &max);   
 
 //2019  void  plot_strom_height(double strom, double s_min, double s_max, double height, double max_h, bool aufdruecken, bool logscaling, short   maxx, short maxy,GdkGC *gc = Application.dots_gc);

  //GdkPoint strom_point(double strom,  
  //                     double s_min, double s_max,
  //  		           double height, double max_h,
  //                    bool logscaling, 
  //                     short maxx, short maxy);

  void  copy_stromsig();
  void  amp_grow(Tampl_histo_type a, int max);

  void  str_filter(bool clear = true);
  void  plot_init();
  void  plot_start();
  void  plot_end();

  void  plot_levels();
  void  plot_strom(double strom);
  bool checkampl(Tampl_fit_type ampl);
  void  updown(int w);
  void  plot_text();

  void  show_ampl_histo(Tampl_histo_type a); 
  void  calc_fit(Tampl_fit_type &af, Tampl_histo_type histo, int min, int max);
  void  call_simplex();
  void  simplex(Tampl_fit_type &af, Tampl_histo_type histo, int min, int max);
  void  plot_levels_and_fit(Tampl_fit_type &a, int max_h, double stromrand = 0, double heightrand = 0);
  void  calc_ahisto();

  void  InitDwellTimeHisto(int n_channels, int n_extra_channels, bool neu);
  void  hinkley(detector_type d);
  void  HOSHD();
  void  HOHD();
  void  HOHD_S(double lambda);
  double schranke(int t_res, int ord, double half_jump_magnitude);
  void  cumsum(g_type &g, int ord);
  void  zero(g_type &g, int ord);
  void  SHD();
  void  SHD_AH();
  void  create_neighbors(Tneighbor_type &n,int lastjump, int open, int extra_open);
  bool flush_Sprung(detector_type d,
                        std::vector<TSprung> &Spruenge,
                        int   &index,
                        int   &new_model,
                        bool jump_up,
                        int   n_open_old,
                        int   n_extra_open_old,
                        int   n_open_new,
                        int   n_extra_open_new);
  //void  plotSpruenge(std::vector<TSprung> &Spruenge,
  //                   GdkColor   Farbe,
  //                   bool all);
  
  void  separatechannels();
  void  saveSChannels();
  
  void  saveDwellTimeHisto(char *name);
  void  saveDetectedLevelHisto(char *name);
  double calcStD(int n_channel, int n_extra_channel);
  void  saveTransitionMatrix(char *name);
  void  clearDetectedLevels();
  void  clearTransitions();
  double findLambda();
  double zeitauflsg(double schr, int ord);

  void  History();

  double DetermineMax(int index);
  TTargetWeight* GetWeight();
  int CalcStartIndex();
  void  NoSamples();
  void  ExpEinteil();
  void  mittel();
  void  axes(double mx, double my, double xs, double ys);
  void  zoomin();
  void  zoomout();

  void  ExpFit();
  char* GetxyString(int &index, char* s);
  void  ResetWeight();
  void  generate_gnu_plot(char* name);
  void  exec_gnuplot(char*  name);
  TTargetWeight*  getweight();

  void  tfit(); //ersatz fuer fit, da fit schon bei TFitDaten!
  void  read_data(TTargetFit* TargetFit);
  void  abspeichern(TTargetFit* TargetFit, Vector lambda[], Vector ampli[]);
  void  resetTargetFitWeight();
  void  rates_To_TimeConstants();
  
  void  test_Matrix();
  void  algAbspeichern(VecD lambda[], VecD ampli[], std::vector<TUntermatrix> &Dimension);
  void  Fehler_Bei_Raten();

  void  calc_dwell_log(std::vector<int> &TDwellC, std::vector<int> &TDwellO, double &min, double &max, int bins); 

private:
  void init(); // fuer konstruktor und reset();

  bool ProcessedText(int i); // fuer str_filter
  void SaveDwell(TDwell &Item, std::ofstream &out);    // SaveDwellTimeHisto
  void FindMax(const TDwell &Item, double &max); // DetermineMax
  void PlotItem(TDwell &Item, int &Startindex, int &i, int &k, TTargetWeight *weight, double &Kreuz, double &old); //Paint Dwell

  void StartIndex(TDwell &Dwell, int &index, double &x, double &y); //calcstart
};

extern TDaten Daten;


#endif	/* _DATEN_H_ */
