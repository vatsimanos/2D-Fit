
/*******************************************************************
 *  Kiel-Patch
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <time.h>
//#include <process.h>	//2019
#include <errno.h>
#include <stdio.h>



//#include "gtkapp.h"
//#include "declare.h"
#include "daten.h"
#include "detector.h"
#include "fit.h"
#include "error.h"
#include "round.h"
//#include "gtkfunc.h"
#include "dwell.h"
#include "setfile.h"



TDaten Daten;             //globale Variablen
TApplication Application;


TDaten::TDaten()
{
  // TDaten
  //Tobias start
  #include"datentobias.c"
  //Tobias end
  
   
  Wert = NULL;            
  cdir[0] = 0;    // noch kein verzeichnis    
  FName[0] = 0;   // noch kein Filename  

  // TPaintDaten
  datatype = unknown;
  schwarzaufweiss = false;

  // Backing pixmap for drawing area
//2019  DataPix = NULL;
//2019  MemPix = NULL; 

  g_maxx = 0;
  g_maxy = 0;

  // neu
  showmarker = false;
  logscalehisto = true;
  
  // wie gross ist der rand bei darstellung der histogramme mit achsen
  rand_pixel = 6;

  reset_level();

  // welches log bin dwell soll zuerst genutzt werden?
  Dwell_1d_ptr = &Dwell_1d_A;
  Dwell_2d_ptr = &Dwell_2d_A;  // 1d und 2d muessen konsistent bleiben! 

  // Variablen initialisieren
  init();
}

TDaten::~TDaten()                    // sollte erst beim Progende sein!
{	
  freeDaten();
//2019  if (DataPix) gdk_pixmap_unref(DataPix); 
//2019  DataPix = NULL;                          
//2019 if (MemPix) gdk_pixmap_unref(MemPix); 
//2019  MemPix = NULL;                            
}

void
TDaten::freeDaten()
{
  // TDaten
  delete[] Wert;            // loescht array (nur wenn ptr!=NULL)


  // TMarkBitsequence loeschen
  //Marker.clear();

  // THistoDaten
  // ahsito und bhisto nicht dynamisch! 

  // TZeitDaten
  // aber neuer bildaufbau? (show bzw paint)  
  // oder besser beim oeffnen?

  // TFilterDaten
  HOHD_Spruenge.clear();
  SHD_Spruenge.clear();
  HOSHD_Spruenge.clear();
  InitDwellTimeHisto(0,0,false);
  // TDwellTimeDaten
  //Dispose(QuenchedDwellTime, Done);
}

void
TDaten::init()
{
// hier werden Variablen mit Werten belegt, brauche ich im 
// Konstruktor und f�r reset() (open datafile)
int i,j,k,m;


  // TDaten
  Anz_Samples = 0;       // noch keine werte
  Settings.Zoom = false; // erstmal kein Zoom

  // TPaintDaten
  graphtype = nothing;  
  command = cbaseline;        

  // THistoDaten
  for (i=0;i < 65536;i++) ahisto[i] = 0;
  for (i=0;i < 65536;i++) bhisto[i] = 0;
  strom_min = 0;
  strom_max = 0;
  
  // TZeitDaten
  Samples_per_Pixel = 1;

  // THistoFitDaten
  fehler = 0;

  // TFilterDaten
  result.filter_ord = 4;
  result.n_bins = n_bins_max;
  candidate_restriction = 1;
  old_index = Settings.SamplesToSkip;
  ScTiRe = 32;  // wei�es Rauschen
  clearDetectedLevels();
  clearTransitions();
  event_too_long = true;
  maximum = 1000000;

  // THistoryDaten 
  // besser zu THistorytext 
  HText.level = 0;
  HText.sublevel = 0;
  HText.detected_jumps = 0;
  HText.sigma_of_channel = 0;
  HText.detector = _HOHD;

  // TDwellTimeDaten
  // besser zu TSDwellTime
  SDwellTime.LevelNo = 0;
  SDwellTime.SublevelNo = 0;
  SDwellTime.GrowFactor = 1.0464753417;
  SDwellTime.Samples = 100;
  SDwellTime.SamplesB = true;
  SDwellTime.open = true;
  SDwellTime.openclosed = false;
  SDwellTime.cross = true;
  SDwellTime.ExpAvarage = true;
  zoomfactor = 1;

  // TExpFitDaten 
  SDwellTime.openclosed = true;
  for (i=0;i <= TargetFitMaxChannels;i++)
  {
    SDwellTime.LevelNo = i;
    ResetWeight();
  }
  SDwellTime.openclosed = false;
  SDwellTime.open = true;
  for (i=0;i <= TargetFitMaxChannels;i++)
  {
    SDwellTime.LevelNo = i;
    ResetWeight();
  }
  SDwellTime.open = false;
  for (i=0;i <= TargetFitMaxChannels ;i++)
  {
    SDwellTime.LevelNo = i;
    ResetWeight();
  }
  SDwellTime.LevelNo = 0;
  for (i=0;i <= n_channels_max+1;i++)
    for (j=0;j <= n_extra_channels_max+1;j++)
      for (m=0;m <= 2;m++)
	  {
        DwellLevel[i][j][m].Error = 0;
        DwellLevel[i][j][m].ErrorCOG = 0;
        DwellLevel[i][j][m].ErrorSCO = 0;
        DwellLevel[i][j][m].ErrorFiveState = 0;
        DwellLevel[i][j][m].iter = 0;
        DwellLevel[i][j][m].Points = 0;
        DwellLevel[i][j][m].startx = 0;
        DwellLevel[i][j][m].IsTargetFit = false;
        for (k=1;k <= NMaxTimeConstants;k++)
        {
          if (k==1)
          {
			DwellLevel[i][j][m].tau[k] = 100;
          } else {
            DwellLevel[i][j][m].tau[k] = 0;
		  }
		  DwellLevel[i][j][m].Ampli[k] = 0;
		}
	  }
  // TDwellfitDaten
  Rates.CO = 0.2;
  Rates.OC = 0.1;
  Rates.OG = 0.5;
  Rates.GO = 1.0;
  modell = M_COG;

  //TAlgModellDwellFitDaten
  for (i=0;i <= n_states_max-1;i++)
  {
    for (j=0;j <= n_states_max-1;j++)
	{
      Rate_Matrix.r[i][j] = 0;
	}
	Rate_Matrix.open[i] = false;
  }
  Rate_Matrix.open[0] = true;
  Rate_Matrix.r[0][1] = 0.2;
  Rate_Matrix.r[0][2] = 0.1;
  Rate_Matrix.r[1][0] = 0.5;
  Rate_Matrix.r[2][0] = 1;

  // TAlgModellDwellFitKoDaten
  Kooperativity = 0; // zu Beginn keine Kooperativit�t vorsehen
 //Tobias: Zur�cksetzten der Markertabelle und initialisierung
 Daten.lauf = Daten.Marker.begin();
 Daten.markernumber = 1;
 Daten.nrselected = 0;
 
 Daten.starting_level = 300;
 Daten.level_increment = 100;
 Daten.number_of_levels= 10;
  
}

void
TDaten::reset_level()
{
// hier werden level auf grundeinstellungen gesetzt
// beim programmstart oder auf benutzerwunsch
  a_fit.n_channels = 1;
  a_fit.n_extra_channels = 0;
  a_fit.i_null = (strom_min + strom_max) / 2;
  a_fit.i_channel = (strom_max-strom_min) / 4;
  a_fit.i_extra_channel = 5;
  a_fit.sigma = 15;

}

void
TDaten::reset()
{
  freeDaten();  // gibt dynamischen speicher frei
  init();  // setzt restlichen Variablen zurueck
}


//BEGIN KARSTENMARKER
void TDaten::get_marker()
  {
  int j = 0;
  unsigned short   suchmark;

Marker.clear();

  // erster marker: j=0
 TMarkBits mark(Wert[0]>>12, 0 , false); // 3.Parameter Erweiterung von Tobias
  Marker.push_back(mark);

  // oberste 4bit nach recht schieben, so dass werte von 0 bis 0xF 
  suchmark = Wert[0] >> 12; // fuer vergleich mit neuem marker!

  Wert[0] &= 0x0fff; //marker in obersten 4bit loeschen
  // weitere Datens�tze nach neuen markern durchsuchen

  for (j=1;j < Anz_Samples;j++)
  { 
    // neuer marker, neuer eintrag
    if (Wert[j]>>12 != suchmark)
    {
      TMarkBits mark(Wert[j] >> 12, j, false); 
      Marker.push_back(mark);
      suchmark=Wert[j]>>12;
    }
    // dann markerbits loeschen
    Wert[j] &= 0x0fff; // marker liegt in obersten 4bit
  }

  }
//END KARSTENMARKER

 
 
int
TDaten::get_data_from_file(char* Filename)
{
// hier wird die zeitserie aus einem '*.0*' File ausgelesen
// und eine linked list fuer tobys markerbits wird angelegt

int j = 0;

//BEGIN KARSTENASCI
if (Daten.import_asci_data==true)
  {
  std::cout << ">> open as asci data" << std::endl;
  std::ifstream in(Filename, std::ios::in|std::ios::binary);
  std::cout << "   Filename: " << Filename << std::endl;
  if (!in)
   {
    char infotext[MaxTextLen];
    sprintf(infotext, "Can't open timeseries: %s", Filename);
    warning(infotext);           // fehlermeldung wohin auch immer 
    return -1;                   // keine daten lesen!
   }
  strcpy(FName, Filename);       // Filename sichern!
  reset();                       // falls schon daten im speicher, ruft auch freeDaten() auf

  in.read((char *)&Anz_Samples, 0);  
  in.seekg(0, std::ios::end);                                // Filelaenge in bytes
  Anz_Samples = in.tellg()  / 6;           //Anzahl Samples: so richtig f�r Asci Daten!
  in.seekg(0, std::ios::beg);                                // streampos auf ersten wert
  Settings.SamplesToProcess = Anz_Samples;

Wert= new unsigned short[Anz_Samples+1];

int importWert;
int importWert1[6],importWert2;
int l;
int anzahl_ziffer;
bool abfrage;

AD_fromwert=1000000;
AD_towert=0;

 Application.statustext("calculating min and max of the asci data...");

for (j=0;j<Anz_Samples;j++)
  {
  anzahl_ziffer=0;
  abfrage=false;
  for(l=0;l<100;l++)
      {
       in.read((char *)&importWert, sizeof(importWert));  
       importWert=(importWert^48);
       importWert2=(int)importWert;
       //cout << importWert2 << endl;
       if ((importWert2<=9)&&(importWert2>=0))
           {
           importWert1[anzahl_ziffer]=(int)importWert;
           anzahl_ziffer++;
           abfrage=true;
           }
       if ((importWert2>9)&&(abfrage==true)) break;
       }
   if (anzahl_ziffer==1) Daten.Wert[j]=importWert1[0];
   if (anzahl_ziffer==2) Daten.Wert[j]=importWert1[0]*10+importWert1[1];
   if (anzahl_ziffer==3) Daten.Wert[j]=importWert1[0]*100+importWert1[1]*10+importWert1[2];
   if (anzahl_ziffer==4) Daten.Wert[j]=importWert1[0]*1000+importWert1[1]*100+importWert1[2]*10+importWert1[3];
   if (anzahl_ziffer==5) Daten.Wert[j]=importWert1[0]*10000+importWert1[1]*1000+importWert1[2]*100+importWert1[3]*10+importWert1[4];

   if (Daten.Wert[j]>AD_towert) AD_towert=Daten.Wert[j];
   if (Daten.Wert[j]<AD_fromwert) AD_fromwert=Daten.Wert[j];
   //cout << "Daten.Wert[" << j << "]=" << Daten.Wert[j] << endl;   
   }

if (AD_towert>4095) 
  {
  Skalierungsvorschlag_min=0;
  Skalierungsvorschlag_max=65535;
  }
else 
  {
  Skalierungsvorschlag_min=0;
  Skalierungsvorschlag_max=4095;
  }

  // N, C, U, K, G ausserhalb der grenzen?
  if (!checkampl(a_fit))
  {
    reset_level();
  }
  Settings.SamplesToSkip = 0;
  Settings.SamplesToProcess = Anz_Samples;

  }
else //END KARSTENASCI
  {
   if (!Daten.Geneticfitparameter.mpi_support)
        std::cout << "loading binary data: " <<Filename;

  // TDaten
  std::ifstream in(Filename, std::ifstream::binary);
   if (!in)
    {
      char infotext[MaxTextLen];
      sprintf(infotext, "Can't open timeseries: %s", Filename);

      warning(infotext);           // fehlermeldung wohin auch immer 
      return -1;                   // keine daten lesen!
    }
  strcpy(FName, Filename);       // Filename sichern!
  reset();                       // falls schon daten im speicher, ruft auch freeDaten() auf

/*  in.read((char *)&Anz_Samples, 8);  
  in.seekg(0, std::ios::end);                                // Filelaenge in bytes
  if (Anz_Samples != (((long int)in.tellg()-8)/2))										//Laden auf Anzahl der vorhandenen Bytes beschr�nken und error, wenn Filel�nge nicht zu Datenanzahl passt
   {
   	std::cout<<std::endl;
   	std::cout<<"error: filelength does not match"<<std::endl;
    Anz_Samples= ((long int)in.tellg()-8)/2;
   }		
  
  in.seekg(8);                                		// streampos auf ersten wert
*/  		
	in.seekg(0, std::ios::end);  										//2023 length tag removed from time series
	Anz_Samples = (long int) in.tellg()/2;	
	std::cout<<"length of time series: "<<Anz_Samples<<std::endl;
	in.seekg(0);	
  
  Daten.Wert = new unsigned short[Anz_Samples];          // pointer auf neue daten
  Settings.SamplesToProcess = Anz_Samples;
  in.read((char *)Daten.Wert, Anz_Samples * 2); 
           

/*BEGIN KARSTENMARKER //Tobias 2022 put off
  marker_Wert = new unsigned short[Anz_Samples];          // pointer auf neue daten
//  in.read(marker_Wert, Anz_Samples * 2);           // werte sind 2byte gross 
for(j=0;j<Anz_Samples;j++)
  {	
  Daten.marker_Wert[j]=Daten.Wert[j];
   }	
get_marker();
*/
 
// noch ein paar einstellungen
min_max();

// N, C, U, K, G ausserhalb der grenzen?
  if (!checkampl(a_fit))
  {
    reset_level();
  }
  Settings.SamplesToSkip = 0;
  Settings.SamplesToProcess = Anz_Samples;
  if (!Daten.Geneticfitparameter.mpi_support)
      std::cout<<" done"<<std::endl;
}

 //BEGIN KARSTENASCI
Daten.import_asci_data=false;
//END KARSTENASCI
  //Daten.plot_text();
  return 0;
}

void 
TDaten::load(char *Filename)
{
/*2019	
// jobfile laden

  // TDaten  
  std::ifstream in(Filename,std::ios::in|std::ios::binary);

  if (!in)
  {
    char infotext[MaxTextLen];
    sprintf(infotext,"Can't open file: %s",Filename);

    warning(infotext);           // fehlermeldung in statusleiste 
    return;                      // keine daten lesen!
  }

  // kennung pruefen
  // mit getline koennten noch fehler auftauchen, wenn \n statt \0
  char c;
  if (((in.get(c)) && (c!= '#')) ||
	  ((in.get(c)) && (c!= 'j')) ||
	  ((in.get(c)) && (c!= 'b')) ||
	  ((in.get(c)) && (c!= 'f')) ||
	  ((in.get(c)) && (c!= '#')) ||
	  ((in.get(c)) && (c!= 0)))
  {
    char infotext[MaxTextLen];
    sprintf(infotext,"%s is not a job-file",Filename);

    warning(infotext);           // fehlermeldung in statusleiste 
    return;                      // keine daten lesen!
  }

  // settings lesen
  char datafilename[fsPathName];

  in.getline(datafilename,fsPathName,'\n');
  // settings auch im selben file
  Settings.load(in);                        
  // oder als eiges ascii file
  //in.getline(datafilename,fsPathName,'\n');
  //ifstream stgin(datafilename);
  //stgin >> Settings;
  // fenster einstellungen!
  gdk_window_move_resize (Application.DLG_Main->window, 
    	                  Settings.x, 
					      Settings.y, 
						  Daten.Settings.w, 
						  Daten.Settings.h);

  // daten aber erst im neuen fenster lesen!
  if (get_data_from_file(datafilename) != 0) return; //Fehler schon behandelt!

  //in.read(&Anz_Samples,sizeof(Anz_Samples));  // steht doch in datafile

  // TPaintDaten
  //in.read(&datatype, sizeof(datatype));       // wird sowieso gesetzt

  // THistoDaten
  in.read((char *)&ahisto, sizeof(int)*4096);        
  //in.read(&bhisto, sizeof(int)*4096);      // wird in get_data erstellt
  in.read((char *)&graphtype, sizeof(graphtype));
  //in.read(&strom_min, sizeof(strom_min));
  //in.read(&strom_max, sizeof(strom_max));

  //in.read(&g_maxx, sizeof(g_maxx));           // wird neu berechnet
  //in.read(&g_maxy, sizeof(g_maxy));           // natuerlich auch

  // TZeitDaten
  //in.read(&Samples_per_Pixel, sizeof(Samples_per_Pixel)); // wird bei plot_init berechnet

  // TFitDaten
  in.read((char *)&a_fit, sizeof(a_fit));

  // THistoFitDaten
  //in.read(&fehler, sizeof(fehler));           // wird bei calc_fit berechnet

  // TFilterDaten
  short   i;
  int   HOHD_Size;
  int   SHD_Size;
  int   HOSHD_Size;
  TSprung  Sprung;
 
  in.read((char *)&ColorAndOrder, sizeof(ColorAndOrder));
  in.read((char *)&HOHD_Size, sizeof(HOHD_Size));
  in.read((char *)&SHD_Size, sizeof(SHD_Size));
  in.read((char *)&HOSHD_Size, sizeof(HOSHD_Size));
  for(i=0;i < HOHD_Size;i++)
  {
    Sprung.load(in);	
    HOHD_Spruenge.push_back(Sprung);
  }
  for(i=0;i < SHD_Size;i++)
  {
    Sprung.load(in);	
    SHD_Spruenge.push_back(Sprung);
  }
  for(i=0;i < HOSHD_Size;i++)
  {
    Sprung.load(in);	
    HOSHD_Spruenge.push_back(Sprung);
  }
  in.read((char *)&result.filter_ord, sizeof(result.filter_ord));
  in.read((char *)&result.n_bins, sizeof(result.n_bins));
  in.read((char *)&candidate_restriction, sizeof(candidate_restriction));
  old_index = 0;
  in.read((char *)&ScTiRe, sizeof(ScTiRe));

  InitDwellTimeHisto(a_fit.n_channels, a_fit.n_extra_channels, true);

  size_t x;
  size_t y;
  in.read((char *)&x, sizeof(x));
  in.read((char *)&y, sizeof(y));
  for (size_t i=0;i < x;i++)
  {
    for (size_t j=0;j < y;j++)
	{
      in.read((char *)&result.dwell_time_histo_array[i][j], sizeof(result.dwell_time_histo_array[i][j]));
	}
  }
  //in.read(&result.dwell_time_histo_array, sizeof(result.dwell_time_histo_array));

  in.read((char *)&result.detected_level_histo, sizeof(result.detected_level_histo));
  in.read((char *)&result.transitions, sizeof(result.transitions));

  event_too_long = false;
  maximum = 1000000;

  // THistoryDaten
  in.read((char *)&HText, sizeof(THistorytext));

  // TDwellTimeDaten
  in.read((char *)&SDwellTime, sizeof(TSDwellTime));

  //QuenchedDwellTime:=New(PCollection, Init(n_bins_max, 1));
  in.read((char *)&zoomfactor, sizeof(zoomfactor));

  // TExpFitDaten.Load
  // Zeitkonstanten u. Amplituden laden
  in.read((char *)&DwellLevel, sizeof(TLevel));
  // hier noch Wichtung lesen
  in.read((char *)&Weightup, sizeof(TTargetWeight));
  in.read((char *)&Weightdown, sizeof(TTargetWeight));
  in.read((char *)&TargetWeight, sizeof(TTargetWeight));

  // TDwellfitDaten
  in.read((char *)&Rates, sizeof(Rates));
  in.read((char *)&modell, sizeof(TModell));

  // TAlgModellDwellFitDaten
  in.read((char *)&Rate_Matrix, sizeof(TRate_matrix));

  // TAlgModellDwellFitKoDaten
  in.read((char *)&Kooperativity, sizeof(double));

  // neu
  in.read((char *)&showmarker, sizeof(showmarker));
  in.read((char *)&schwarzaufweiss, sizeof(schwarzaufweiss));
  in.read((char *)&logscalehisto, sizeof(logscalehisto));

  // noch ein paar berechnungen fuer dwell time histogramme
  if (graphtype == dwell_time)
  {
    if (SDwellTime.SamplesB) // entweder Samples oder Grow-Factor Berechnen
	{
      ExpEinteil();
	} else {
      Daten.NoSamples();
	}
    Daten.mittel(); // Exp. Mittelung
  }

  // Dwell lb
  Dwell_1d_A.load(in);
  Dwell_2d_A.load(in);
  Dwell_1d_B.load(in);
  Dwell_2d_B.load(in);
2019*/  
}

void 
TDaten::store(char *Filename)
{
// jobfile speichern 
// anderes format als in PASCAL

  // TDaten
  std::ofstream out(Filename,std::ios::out|std::ios::binary);

  if (!out)
  {
    char infotext[MaxTextLen];
    sprintf(infotext,"Can't open file: %s",Filename);

    warning(infotext);           // fehlermeldung wohin auch immer
    return;                      // keine daten lesen!
  }

  // datafilename schreiben
  out << "#jbf#" << '\0';        // kennung schreiben: #jbf#0!
  out << FName << '\n';		     // beim lesen auf endline testen!

  // settings schreiben
  // binaer
  Settings.store(out);
  // oder ascii file
  //char stgname[MaxTextLen];
  //char* p;
  //strcpy(stgname,Filename);
  //p = strrchr(stgname,'.');
  //if (p != NULL) *p = 0; // endung abschneiden 
  //strcat (stgname, ".ini");
  //out << stgname << '\n';    // beim lesen auf endline testen!
  //ofstream stgout(stgname);
  //stgout << Settings;
  //stgout.close();
  
  //  out.write(&Anz_Samples,sizeof(Anz_Samples)); //warum ? steht doch in datafile

  // TPaintDaten
  //out.write(&datatype, sizeof(datatype));

  // THistoDaten
  out.write((char *)&ahisto, sizeof(int)*4096);
  //out.write(&bhisto, sizeof(int)*4096);   // wird von get_data erzeugt
  out.write((char *)&graphtype, sizeof(graphtype));
  //out.write(&strom_min, sizeof(strom_min));  // wird auch von get_data
  //out.write(&strom_max, sizeof(strom_max));  // erzeugt

  //out.write(&g_maxx, sizeof(g_maxx));        // wird neu berechnet
  //out.write(&g_maxy, sizeof(g_maxy));        // natuerlich auch

  // TZeitDaten
  //out.write(&Samples_per_Pixel, sizeof(Samples_per_Pixel)); // wird bei plot_init berechnet

  // TFitDaten
  out.write((char *)&a_fit, sizeof(a_fit));

  // THistoFitDaten
  //out.write(&fehler, sizeof(fehler));        // wird bei calc_fit berechnet

  // TFilterDaten
  int   SHD_Size;
  int   HOHD_Size;
  int   HOSHD_Size;
  std::vector<TSprung>::iterator iter;

  out.write((char *)&ColorAndOrder, sizeof(ColorAndOrder));
  SHD_Size = SHD_Spruenge.size();
  HOHD_Size = HOHD_Spruenge.size();
  HOSHD_Size = HOSHD_Spruenge.size();
  out.write((char *)&HOHD_Size, sizeof(HOHD_Size));
  out.write((char *)&SHD_Size, sizeof(SHD_Size));
  out.write((char *)&HOSHD_Size, sizeof(HOSHD_Size));
  for(iter=HOHD_Spruenge.begin();iter != HOHD_Spruenge.end();iter++)
    (*iter).store(out);
  for(iter=SHD_Spruenge.begin();iter != SHD_Spruenge.end();iter++)
    (*iter).store(out);
  for(iter=HOSHD_Spruenge.begin();iter != HOSHD_Spruenge.end();iter++)
    (*iter).store(out);
  out.write((char *)&result.filter_ord, sizeof(result.filter_ord));
  out.write((char *)&result.n_bins, sizeof(result.n_bins));
  out.write((char *)&candidate_restriction, sizeof(candidate_restriction));
  out.write((char *)&ScTiRe, sizeof(ScTiRe));

  size_t x = result.dwell_time_histo_array.sizex();
  size_t y = result.dwell_time_histo_array.sizey();
  out.write((char *)&x, sizeof(x));
  out.write((char *)&y, sizeof(y));
  for (size_t i=0;i < x;i++)
  {
    for (size_t j=0;j < y;j++)
	{
      out.write((char *)&result.dwell_time_histo_array[i][j], sizeof(result.dwell_time_histo_array[i][j]));
	}
  }

  out.write((char *)&result.detected_level_histo, sizeof(result.detected_level_histo));
  out.write((char *)&result.transitions, sizeof(result.transitions));

  // THistoryDaten
  out.write((char *)&HText, sizeof(THistorytext));

  // TDwellTimeDaten
  out.write((char *)&SDwellTime, sizeof(TSDwellTime));
  out.write((char *)&zoomfactor, sizeof(zoomfactor));

  // TExpFitDaten.Store
  // Alle Zeitkonstanten u. Amplituden sichern
  out.write((char *)&DwellLevel, sizeof(TLevel));
  // hier noch Wichtung schreiben
  out.write((char *)&Weightup, sizeof(TTargetWeight));
  out.write((char *)&Weightup, sizeof(TTargetWeight));
  out.write((char *)&TargetWeight, sizeof(TTargetWeight));

  // TDwellFitDaten
  out.write((char *)&Rates, sizeof(Rates));
  out.write((char *)&modell, sizeof(TModell));

  // TAlgModellDwellFitDaten
  out.write((char *)&Rate_Matrix, sizeof(TRate_matrix));

  // TAlgModellDwellFitKoDaten
  out.write((char *)&Kooperativity, sizeof(double));

  // neu
  out.write((char *)&showmarker, sizeof(showmarker));
  out.write((char *)&schwarzaufweiss, sizeof(schwarzaufweiss));
  out.write((char *)&logscalehisto, sizeof(logscalehisto));

  // Dwell lb
  Dwell_1d_A.store(out);
  Dwell_2d_A.store(out);
  Dwell_1d_B.store(out);
  Dwell_2d_B.store(out);
  //differenz wird neu berechnet!
}

void 
TDaten::markeddatafit()
{
// still under construction (toby)
}


void
TDaten::PlotItem(TDwell &Item, int &Startindex, int &i, int &k, TTargetWeight *weight, double &Kreuz, double &old)
{
/*2019	
// plotet entweder punkte oder balken in dwell time diagramm

  // Plot Balkendiagramm
  if (!(SDwellTime.cross))
  {
    plot_strom_height(Item.x + xstart,0,maxx, ystart, maxy, false, false, g_maxx, g_maxy, Application.draw_gc);
    plot_strom_height(Item.x + xstart,0,maxx, Item.y + ystart, maxy, true, false, g_maxx, g_maxy, Application.draw_gc);
    plot_strom_height(old,0,maxx, Item.y+ystart, maxy, true, false, g_maxx, g_maxy, Application.draw_gc);
    plot_strom_height(old,0,maxx, ystart, maxy, true, false, g_maxx, g_maxy, Application.draw_gc);
  } else {
  // Plot Kreuze
    i++;

    if (i > Startindex)
    {
      k++;
      if ((*weight)[SDwellTime.LevelNo][k] != 1) 
      {
		// Farbe nach Gewichtung bestimmen, Idee von Maike
        if ((*weight)[SDwellTime.LevelNo][k] > 1)
        {
          Application.drawcol(RGB2Gdk(R1, G1, B1));
        } else {
          Application.drawcol(RGB2Gdk(R2, G2, B2));
        }
      }

    } else {
      Application.drawcol(RGB2Gdk(R3, G3, B3));
    }
    Kreuz = maxx / 200;
    plot_strom_height(Item.x+xstart-Kreuz, 0, maxx, Item.y+ystart, maxy, false, false, g_maxx, g_maxy, Application.draw_gc);
    plot_strom_height(Item.x+xstart+Kreuz, 0, maxx, Item.y+ystart, maxy, true, false, g_maxx, g_maxy, Application.draw_gc);
    
	Kreuz = maxy / 200;
    plot_strom_height(Item.x+xstart, 0, maxx, Item.y+ystart-Kreuz, maxy, false, false, g_maxx, g_maxy, Application.draw_gc);
    plot_strom_height(Item.x+xstart, 0, maxx, Item.y+ystart+Kreuz, maxy, true, false, g_maxx, g_maxy, Application.draw_gc);

    if ((i <= Startindex) || ((*weight)[SDwellTime.LevelNo][k] != 1))
    {
      // urspr�ngliche Farbe wiederherstellen
      Application.drawcol(RGB2Gdk(HOHD_default_Color_R, HOHD_default_Color_G, HOHD_default_Color_B));
    }
  }
  old = Item.x + xstart;
2019*/  
}

void fill_rect_rgb(unsigned char* buffer, int image_width,
				   int x_start, int y_start, 
				   int width, int height,
                   int Red, int Green, int Blue)
{
// einen rechteckigen bereich im puffer fuellen!
// fuer 2d dwell time ausgabe
  for (int y = y_start; y < y_start + height;y++)
  {	
    unsigned char* pos = buffer + (y * image_width + x_start) * 3;
    for (int x = x_start; x < x_start + width;x++)
	{
	  *pos++ = Red;
      *pos++ = Green;
	  *pos++ = Blue;
    }
  }
}


void 
TDaten::paint() 
{
/*2019	
// ausgabe des gewuenschten diagramms in mempix 
// hier kann man noch an client server aufrufen optimieren, wenn man 
// verschiedene gdk_draw befehle in 
// gdk_graw_segments oder gdk_draw_rgb_image aufrufen zusammenfasst

  // falls fenster minimiert ist
  if ((g_maxx <= 1) || (g_maxy <= 1)) return;

  if ((datatype == unknown) ||
	  ((graphtype == time_series) && (datatype != single)) ||
	  ((graphtype == multi_time) && (datatype != multi)))
  {
	// Daten samt filtern neuzeichnen, time_series oder multi_time
    plot(); 
  }

  switch (graphtype)
  {
    case growing_ampl_histo:
      // Hintergrund mit back_gc fuellen
      gdk_draw_rectangle(MemPix, Application.back_gc, true,
                 		 0, 0, g_maxx, g_maxy);
      gdk_draw_line(MemPix, Application.border_gc,
					0, g_maxy-1, g_maxx, g_maxy-1);
      copy_stromsig();
      break;  

    case start_ampl_histo:
      // Hintergrund mit back_gc fuellen
      gdk_draw_rectangle(MemPix, Application.back_gc, true,
                 		 0, 0, g_maxx, g_maxy);

	  show_ampl_histo(bhisto);
	  calc_fit(a_fit, bhisto, strom_min, strom_max);
      plot_levels_and_fit(a_fit, max_height);
	  break;
 
    case ampl_histo:
      // Hintergrund mit back_gc fuellen
      gdk_draw_rectangle(MemPix, Application.back_gc, true,
                 		 0, 0, g_maxx, g_maxy);

      show_ampl_histo(ahisto);
      calc_fit(a_fit, ahisto, strom_min, strom_max);
      plot_levels_and_fit(a_fit, max_height);
      break;

    case time_series:
	{
	short  i;
    short  window_height;
 
    int  j, k, l, x, y;
    double t1, t2;

    short  nwindows;
	
		
	  // Hintergrund mit back_gc fuellen
      // windows, marker, focus loeschen, werden unten sowieso geloescht!
      // gdk_draw_rectangle(MemPix, Application.back_gc, true,
	  //                    0, 0, g_maxx, g_maxy);

	  // daten in Mempix zeichnen
      gdk_draw_pixmap(Daten.MemPix, Application.dots_gc,
  	                  Daten.DataPix, 0,0,0,0,g_maxx,g_maxy);

//BEGIN KARSTENCUT_PASTE_EDIT        
                
if (Daten.refresh==false)
      {
        nwindows = Daten.Settings.StackedWindows();

        // invertieren
        gdk_gc_set_function(Application.dots_gc, GDK_INVERT);

        // im selben stacked window
        if ((Application.yMouse2 - Application.yMouse1) == (Daten.g_maxy / nwindows))
        {
          gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
	  	                     Application.xMouse1, 
                             Application.yMouse1, 
                             Application.xMouse2 - Application.xMouse1, 
                             Application.yMouse2 - Application.yMouse1);
        } else {
		  // von oben nach unten
          if ((Application.yMouse2 - Application.yMouse1) > (Daten.g_maxy / nwindows))
          {
            gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
                               Application.xMouse1, 
                               Application.yMouse1, 
                               Daten.g_maxx - Application.xMouse1, 
                               Daten.g_maxy / nwindows -1);
            gdk_draw_rectangle(Daten.MemPix,
                               Application.dots_gc,
                               true,
                               0, 
                               Application.yMouse2 - (Daten.g_maxy / nwindows),
                               Application.xMouse2,
                               Daten.g_maxy / nwindows);

            if ((Application.yMouse2 - Application.yMouse1) > (2*(Daten.g_maxy / nwindows)))
            {
              gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
                                 0, 
                                 Application.yMouse1 + (Daten.g_maxy / nwindows), 
                                 Daten.g_maxx, 
                                 Application.yMouse2 - (Daten.g_maxy / nwindows)-1 -(Application.yMouse1 + (Daten.g_maxy / nwindows)));
            }
          } else {
            // von unten nach oben
            gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
                               0, 
                               Application.yMouse1, 
                               Application.xMouse1, 
                               Daten.g_maxy / nwindows -1);
            gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
     	    	               Application.xMouse2, 
                               Application.yMouse2 - (Daten.g_maxy / nwindows),
                               Daten.g_maxx - Application.xMouse2,
                               Daten.g_maxy / nwindows);
            if ((Application.yMouse1 - Application.yMouse2) > 0)
            {
              gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
	    	                     0, 
                                 Application.yMouse2, 
                                 Daten.g_maxx, 
                                 Application.yMouse1 - Application.yMouse2);
            }
          }
        }
        // gc zuruecksetzten
        gdk_gc_set_function(Application.dots_gc, GDK_COPY);    
      }                    
//END KARSTENCUT_PASTE_EDIT


      // stackedwindows plotten
      window_height = g_maxy / Settings.StackedWindows();
      for (i=1;i<Settings.StackedWindows()+1;i++)
      {
        gdk_draw_line(Daten.MemPix, Application.border_gc,
                      0, g_maxy - window_height*i,
                      g_maxx, g_maxy - window_height*i );
      }
      gdk_draw_line(Daten.MemPix, Application.border_gc,
                    0, g_maxy -2,
                    g_maxx-1, g_maxy -2);
      gdk_draw_line(Daten.MemPix, Application.border_gc,
                    0, g_maxy - window_height*Settings.StackedWindows(),
                    0, g_maxy);
      gdk_draw_line(Daten.MemPix, Application.border_gc, 
                    g_maxx-1, g_maxy - window_height*Settings.StackedWindows(),
                    g_maxx-1, g_maxy);

      // marker zeichnen
      if (showmarker)
      {
        // Create constant iterator for list
        std::vector<TMarkBits>::const_iterator iter;
	
	//BEGIN KARSTENAVARAGEFILTER
	if (Settings.AvarageFilter()==gleitend)
		{
		t1 = Samples_per_Pixel;		
		}
	else
		{
		t1 = Samples_per_Pixel * Settings.LengthOfAvarageFilter(); // sind dann orig_sam_p_Pix
		}
        //t1 = Samples_per_Pixel * Settings.LengthOfAvarageFilter(); // sind dann orig_sam_p_Pix
        //END KARSTENAVARAGEFILTER
        
        k = g_maxy / Settings.StackedWindows(); // hoehe eines stackedwindows
        j = g_maxy % Settings.StackedWindows(); // restliche freie Pixel ueber stackedwindows

        l = k / 2;  // l: halbe markerhoehe

        for (iter=Marker.begin(); iter != Marker.end(); iter++)
        {
          t2 = (iter->startpos-Settings.SamplesToSkip) / t1;
          x = (int)round(t2) % g_maxx;  // x wert des markers
          i = round(t2) / g_maxx;  // in stackedwindow 0..n-1
    
		  // liegt der marker im darstellbaren bereich
          if ((x<=g_maxx) && (x>=0) &&
              (i<Settings.StackedWindows()) && (i>=0))
          {
            y = round((i+0.5)*k);
            gdk_draw_line(Daten.MemPix,Application.mark_gc,
                          x, y-l+j,
                          x, y+l+j);
          }
        }
      }

      // TFitDaten
      // hier werden die levels gezeichnet
      plot_levels();
      // draw_mouse_focus 
      // button1 innerhalb drawing area pressed 
      if (Application.BtnDownFirst)
      {
        nwindows = Daten.Settings.StackedWindows();

        // invertieren
        gdk_gc_set_function(Application.dots_gc, GDK_INVERT);

        // im selben stacked window
        if ((Application.yMouse2 - Application.yMouse1) == (Daten.g_maxy / nwindows))
        {
          gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
	  	                     Application.xMouse1, 
                             Application.yMouse1, 
                             Application.xMouse2 - Application.xMouse1, 
                             Application.yMouse2 - Application.yMouse1);
        } else {
		  // von oben nach unten
          if ((Application.yMouse2 - Application.yMouse1) > (Daten.g_maxy / nwindows))
          {
            gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
                               Application.xMouse1, 
                               Application.yMouse1, 
                               Daten.g_maxx - Application.xMouse1, 
                               Daten.g_maxy / nwindows -1);
            gdk_draw_rectangle(Daten.MemPix,
                               Application.dots_gc,
                               true,
                               0, 
                               Application.yMouse2 - (Daten.g_maxy / nwindows),
                               Application.xMouse2,
                               Daten.g_maxy / nwindows);

            if ((Application.yMouse2 - Application.yMouse1) > (2*(Daten.g_maxy / nwindows)))
            {
              gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
                                 0, 
                                 Application.yMouse1 + (Daten.g_maxy / nwindows), 
                                 Daten.g_maxx, 
                                 Application.yMouse2 - (Daten.g_maxy / nwindows)-1 -(Application.yMouse1 + (Daten.g_maxy / nwindows)));
            }
          } else {
            // von unten nach oben
            gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
                               0, 
                               Application.yMouse1, 
                               Application.xMouse1, 
                               Daten.g_maxy / nwindows -1);
            gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
     	    	               Application.xMouse2, 
                               Application.yMouse2 - (Daten.g_maxy / nwindows),
                               Daten.g_maxx - Application.xMouse2,
                               Daten.g_maxy / nwindows);
            if ((Application.yMouse1 - Application.yMouse2) > 0)
            {
              gdk_draw_rectangle(Daten.MemPix, Application.dots_gc, true,
	    	                     0, 
                                 Application.yMouse2, 
                                 Daten.g_maxx, 
                                 Application.yMouse1 - Application.yMouse2);
            }
          }
        }
        // gc zuruecksetzten
        gdk_gc_set_function(Application.dots_gc, GDK_COPY);    
      } 
      break;
    }
    case dwell_time:
	{
	// Paint Routine f�r Dwell-Time Histogramme. Achtung: Z�hlung der Histogramme    
    // beachten, die korrekte Reihenfolge ist:                                       
    //   1. Sprung nach oben (up), offene Kan�le                                     
    //   2. Sprung nach unten (down), geschlossene Kan�le                            
    //   a) das Nullniveau                                                           
    //   b) alle Sublevel des Nullniveaus                                            
    //   c) Das 1.Stromniveau                                                        
    //   e) alle Sublevel des ersten Stromnivaus                                     
    //   f) u.s.w 

    double  Kreuz = 0;
    double  old;
    TTargetWeight*  weight;
    short   Startindex;

	short  i;
    int  j, k;

    double  x_d, xx, y_d;   // x -> x_d, y-> y_d: paint exp fit
    short   distance;


	  // hintergrund fuellen
      gdk_draw_rectangle (MemPix, Application.back_gc, true,
	                      0, 0, g_maxx, g_maxy);

	  j = (int)(!(SDwellTime.open)) + 2 * SDwellTime.SublevelNo + 2 * SDwellTime.LevelNo * (a_fit.n_extra_channels + 1);
      if ((j <= 2 * (a_fit.n_channels + 1) * (a_fit.n_extra_channels + 1)) && (!Daten.result.dwell_time_histo_array.empty()))
      {
		maxy = DetermineMax(j);
        ystart = maxy / 20;
        maxy = maxy + ystart;
        maxx = result.n_bins / zoomfactor;
        xstart = maxx / 20;
        old = xstart;
        maxx = maxx + xstart;

        // Achsen zeichnen
        axes(maxx, maxy, xstart, ystart);

        // Farbe f�r Datenpunkt wechseln
    	Application.drawcol(RGB2Gdk(HOHD_default_Color_R, HOHD_default_Color_G, HOHD_default_Color_B));
		
        // Punkte zeichen
        plot_strom_height(xstart,0,maxx,ystart,maxy,false,false, g_maxx, g_maxy, Application.draw_gc);
        i = 0;
        k = 0;
        weight = GetWeight();
        Startindex = CalcStartIndex();

        std::vector<TDwell>::iterator iter;
        for (iter=QuenchedDwellTime.begin(); iter != QuenchedDwellTime.end(); iter++)
		{
		  PlotItem(*iter, Startindex, i, k, weight, Kreuz, old);
		}
   
        // Zeichnet Exp.Fkt in Dwell-Time-Histogram, 
	    // falls schon gefittet
        k = 0;
        if (SDwellTime.open) k = opening;
        if (SDwellTime.openclosed) k = openclosing;
        if (DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].Ampli[1] != 0)
		{
	      //Ist das gew�hlte Histogram schon gefittet?
          // Zeitkonstanten festlegen
          TExpFit E_fit(DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].tau);

          for (i=1;i <= NMaxTimeConstants;i++)
		  {
            E_fit.aa[i] = DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].Ampli[i]; // Amplituden festlegen
		  }
          distance = result.n_bins / 300; // Rasterrung, mit der die e-Fkt gezeichnet wird

          // gerader Strich der den Start markiert (soll warum auch immer verschwinden->sagt Ulf)
		  // y_d = 0;
          // x_d = DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].startx;
          // plot_strom_height(x_d+xstart,0,maxx,y_d+ystart, maxy, false, false, g_maxx, g_maxy, Application.current_gc);
          // y_d = E_fit.summe_expo(x);
          // plot_strom_height(x_d+xstart,0,maxx,y_d+ystart, maxy, true, false, g_maxx, g_maxy, Application.border_gc);

          // fahre zum Anfangswert, kein gerader Strich
		  x_d = DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].startx;
		  y_d = E_fit.summe_expo(x_d);
          plot_strom_height(x_d + xstart, 0, maxx, y_d + ystart, maxy, false, false, g_maxx, g_maxy, Application.current_gc);

		  // zeichnen ab anfangswert, PASCAL zeichnet ab distance
		  for (xx=x_d;xx <= result.n_bins;xx+=distance)
		  {
            y_d = E_fit.summe_expo(xx);  // Summe der Exponentialfunktionen der obigen Amplituden und Zeitkonstanten
            if ((y_d < maxy) && (xx < maxx) && (xx > x_d))
			{
              plot_strom_height(xx + xstart, 0, maxx, y_d + ystart, maxy, true, false, g_maxx, g_maxy, Application.current_gc);
			}
		  }
          // Zeitkonstanten sichern
          E_fit.done(DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].tau);
		}
	  }
	  break;
	}
    case multi_time:
	{
	short window_height;

	  // Hintergrund mit back_gc fuellen
      gdk_draw_rectangle (MemPix, Application.back_gc, true,
		                  0, 0, g_maxx, g_maxy);

      // daten in Mempix zeichnen
      gdk_draw_pixmap(Daten.MemPix, Application.dots_gc,
                      Daten.DataPix, 0,0,0,0,g_maxx,g_maxy);  // ganze area!
      
	  // rahmen zeichnen
      window_height = g_maxy / Settings.StackedWindows();
      gdk_draw_rectangle (MemPix, Application.border_gc, false,
                          0,
                          g_maxy - window_height * Settings.StackedWindows(),
                          g_maxx-1,
                          window_height * Settings.StackedWindows()-1);    
      break;
    }
	case dwell_diff:
	  // dwell 2d mit diff aufrufen und _ptr wieder auf A oder B
      on_CM_log_dwell_diff_activate(NULL, NULL);
	  break;
    case dwell_2d:
	{
    int  j, k, x, y, sizex, sizey;
	//guchar* rgbbuf;
    double  factor;
    unsigned int    R, G, B;
    double  min = 0, max = 0;
    bool grey;
    double  x_d, y_d;   // x -> x_d, y-> y_d: paint exp fit

//	  rgbbuf = new guchar[g_maxx * g_maxy * 3];

	  // wenn schwarzaufweiss, dann grau (1bit macht wohl keinen sinn!)
      grey = Daten.schwarzaufweiss;

      Daten.schwarzaufweiss = false;
	  if (grey)
	  {
        Application.drawcol(RGB2Gdk(255, 255, 255));
//		fill_rect_rgb(rgbbuf, g_maxx,
//			          0, 0, g_maxx, g_maxy,
//                      255, 255, 255);

	  } else {
        Application.drawcol(RGB2Gdk(127, 127, 127));
//		  fill_rect_rgb(rgbbuf, g_maxx,
//         			    0, 0, g_maxx, g_maxy,
//						127, 127, 127);

	  }

	  // Hintergrund mit draw_gc fuellen
	  gdk_draw_rectangle (Daten.MemPix, Application.draw_gc, true,
                          0, 0, g_maxx, g_maxy);
      x_d = g_maxx - 50;
	  y_d = g_maxy - 50;

	  // sind ueberhaupt daten zu zeichnen?
      if (!(Dwell_2d_ptr->CO.empty()))
	  {
        if (((double)Dwell_2d_ptr->CO.sizex() / Dwell_2d_ptr->CO.sizey()) < (x_d / y_d))
		{
          // skalieren nach y
	      factor = y_d / Dwell_2d_ptr->CO.sizey();

	      // skala rechts
	      if (x_d - Dwell_2d_ptr->CO.sizex() * factor > 30)
		  {
		    for (y=5; y < g_maxy-9;y++)
			{
		      ZToRGB(y, g_maxy-10, 5, R, G, B, grey); 
//              fill_rect_rgb(rgbbuf, g_maxx,
//	                        g_maxx - 25, y, 20 , 1,
//                            R, G, B);
              Application.drawcol(RGB2Gdk(R, G, B));
		      gdk_draw_line(Daten.MemPix, Application.draw_gc, 
		   	                g_maxx - 25, y,  g_maxx - 5, y);
			}
		  }	 
		} else {
          // skalieren nach x
	      factor = x_d / Dwell_2d_ptr->CO.sizex();	   

   	      // skala unten
   	      if (y_d - Dwell_2d_ptr->CO.sizey() * factor > 30)
		  {
		    for (x=5; x < g_maxx-9;x++)
			{
		      ZToRGB(x, 5, g_maxx-10, R, G, B, grey); 
//              fill_rect_rgb(rgbbuf, g_maxx,
//	                        x, g_maxy - 25, 1, 20,
//                            R, G, B);

              Application.drawcol(RGB2Gdk(R, G, B));
		      gdk_draw_line(Daten.MemPix, Application.draw_gc,
		   	                x, g_maxy -25, x, g_maxy -5);
			}
		  }
		}
		if (Dwell_1d_ptr != NULL)
		{
          // 1d x von links nach rechts
          // min, max ermitteln 
	      // wenn open und close getrennt, dann besser in TDwell_1d! 
          if (!Dwell_1d_ptr->Close.empty())
		  {
            min = max = Dwell_1d_ptr->Close[0];
		  }

	      for (j=0;j < (int)(Dwell_1d_ptr->Close.size());j++)
		  {
            if (Dwell_1d_ptr->Close[j] < min) min = Dwell_1d_ptr->Close[j]; 
            if (Dwell_1d_ptr->Close[j] > max) max = Dwell_1d_ptr->Close[j]; 
		  }
		  if (min == max) max++; // wenn keine daten, trotzdem bereich von 0 bis 1!

		  for (j=0;j < (int)(Dwell_1d_ptr->Close.size());j++)
		  {
            if (Dwell_2d_ptr->sqrt_e())
			{
			  ZToRGB(sqrt(Dwell_1d_ptr->Close[j]), sqrt(min), sqrt(max), R, G, B, grey); 
			} else {
              ZToRGB(Dwell_1d_ptr->Close[j], min, max, R, G, B, grey); 
			}
            Application.drawcol(RGB2Gdk(R, G, B));

            x = round(Dwell_2d_ptr->CO.sizex() * factor / Dwell_1d_ptr->Close.size() * j);

//            fill_rect_rgb(rgbbuf, g_maxx,
//                          x + 5, (short)(Dwell_2d_ptr->CO.sizey() * factor + 5 + 10),
//                          (short)(Dwell_2d_ptr->CO.sizex() * factor), 30,
//                          R, G, B);
    	    gdk_draw_rectangle  (Daten.MemPix,
                                 Application.draw_gc,
			    				 true,
							     x + 5,
						         (short)(Dwell_2d_ptr->CO.sizey() * factor + 5 + 10),
                                 (short)(Dwell_2d_ptr->CO.sizex() * factor + 5 - x),
                                 30);
		  }

		  // 1d y von unten nach oben 
          //min, max ermitteln
	      // wenn open und close getrennt, dann besser in TDwell_1d! 
          if (!Dwell_1d_ptr->Open.empty())
  	        min = max = Dwell_1d_ptr->Open[0];
	      for (k=0;k < (int)(Dwell_1d_ptr->Open.size());k++)
		  {
            if (Dwell_1d_ptr->Open[k] < min) min = Dwell_1d_ptr->Open[k]; 
            if (Dwell_1d_ptr->Open[k] > max) max = Dwell_1d_ptr->Open[k]; 
		  }
  		  if (min == max) max++; // wenn keine daten, trotzdem bereich von 0 bis 1!

	      for (k=(int)(Dwell_1d_ptr->Open.size()-1);k >= 0;k--)
		  {
            if (Dwell_2d_ptr->sqrt_e())
			{
              ZToRGB(sqrt(Dwell_1d_ptr->Open[k]), sqrt(min), sqrt(max), R, G, B, grey); 
			} else {
              ZToRGB(Dwell_1d_ptr->Open[k], min, max, R, G, B, grey); 
			}
            Application.drawcol(RGB2Gdk(R, G, B));

            y = round(Dwell_2d_ptr->CO.sizey() * factor * (1 - (double)k / (Dwell_1d_ptr->Open.size()-1)));

            gdk_draw_rectangle (Daten.MemPix,
                                Application.draw_gc,
                                true,
		     	    		    (short)(Dwell_2d_ptr->CO.sizex() * factor + 5 + 10),
						        y + 5,
						        30,
                                (short)(Dwell_2d_ptr->CO.sizey() * factor + 5 - y));
		  }
		}

  	    // 2d
		sizex = (int)(Dwell_2d_ptr->CO.sizex());
		sizey = (int)(Dwell_2d_ptr->CO.sizey());
   	    for (j=0;j < sizex;j++)
		{
    	  for (k=sizey-1;k >= 0;k--)
		  {
			if ((Dwell_2d_ptr->sqrt_e()) && (Dwell_2d_ptr != &Dwell_2d_diff))
			{
              min = Dwell_2d_ptr->min_e(false);
              max = Dwell_2d_ptr->max_e(false); // !sym wenn autorange

			  // wenn keine daten, trotzdem bereich von 0 bis 1! damit alles schwarz!
              if (min == max) max++; 
     		  ZToRGB(sqrt(Dwell_2d_ptr->CO[j][k]), sqrt(min), sqrt(max), R, G, B, grey); 
			} else {
			  bool diff = (Dwell_2d_ptr == &Dwell_2d_diff);

			  min = Dwell_2d_ptr->min_e(diff);
              max = Dwell_2d_ptr->max_e(diff);  // sym wenn autorange und diff

              if (!diff) if (min == max) max++; // [0:0] gruen fuer diff, schwarz fuer A oder B

              ZToRGB(Dwell_2d_ptr->CO[j][k], min, max, R, G, B, grey); 
			}
            Application.drawcol(RGB2Gdk(R, G, B));

            x = round(factor * j);
            y = round((sizey-1 - k) * factor);
 
//    		fill_rect_rgb(rgbbuf, g_maxx,
//                          x + 5, y + 5, 
//                          (short)factor + 1, (short)factor + 1,
//						  R, G, B);

       	    gdk_draw_rectangle (Daten.MemPix,
                                Application.draw_gc,
                                true,
                                x + 5, 
				    		    y + 5, 
                                (short)factor + 1,
                                (short)factor + 1);
		  }
		}
	  }
      Daten.schwarzaufweiss = grey;  // alten zustand wieder herstellen

	  // rgb ausgeben!
//	  gdk_draw_rgb_image(Application.drawingarea->window,
//					     Application.dots_gc,
//                         0, 0, g_maxx, g_maxy,
//                         GDK_RGB_DITHER_MAX, rgbbuf, g_maxx * 3);
//	  gdk_draw_rgb_image(MemPix,
//					     Application.dots_gc,
//                         0, 0, g_maxx, g_maxy,
//                         GDK_RGB_DITHER_MAX, rgbbuf, g_maxx * 3);
	  // puffer wieder freigeben
//	  delete[] rgbbuf; 
      break;
	}
    case nothing:
		// Hintergrund loeschen und linie ueber statusbar zeichnen 
		gdk_draw_rectangle(MemPix, Application.back_gc, true,
		                   0, 0, g_maxx, g_maxy);
		gdk_draw_line(MemPix, Application.border_gc,
					  0, g_maxy-1, g_maxx, g_maxy-1);
	  break;
    default:
      break;
  }
  // MemPix in drawingarea ausgeben
  gdk_draw_pixmap(Application.drawingarea->window, Application.dots_gc,
                  Daten.MemPix, 0,0,0,0,g_maxx,g_maxy);  // ganze area!
  // und statusbar aktualisieren
  plot_text();
  Daten.paintinaction = false; // Tobias hebt die Tastatursperre wieder auf

//BEGIN KARSTENCUT_PASTE_EDIT	
Daten.refresh=true; 
//END KARSTENCUT_PASTE_EDIT  
2019*/
}

void 
TDaten::plot()
{


//BEGIN KARSTENNEW
if (Settings.enlarge_current==true)
  {
  strom_min=Minimaler_Strom;   //�nderung hier!
  strom_max=Maximaler_Strom;   //�nderung hier!
  }
else
  {
  strom_min=Settings.CurrentFrom;
  strom_max=Settings.CurrentTo;
  }
//END KARSTENNEW


// Daten in datapix zeichnen (mit str_filter)
TSettings Settings_back;
int multi_times;
/*
  switch (graphtype)
  {
    case time_series:
	  datatype = single;
      // daten ploten
      str_filter();
      // TFilterDaten
      // erkannte Spruenge zeichnen
      switch (ColorAndOrder.DrawingOrder)
	  {
        case 0:
          if ((!HOHD_Spruenge.empty()) && (ColorAndOrder.HOHD_Showflag))
            plotSpruenge(HOHD_Spruenge, RGB2Gdk(ColorAndOrder.HOHD_R,ColorAndOrder.HOHD_G,ColorAndOrder.HOHD_B), true);
          if ((!SHD_Spruenge.empty()) && (ColorAndOrder.SHD_Showflag))
            plotSpruenge(SHD_Spruenge, RGB2Gdk(ColorAndOrder.SHD_R,ColorAndOrder.SHD_G,ColorAndOrder.SHD_B), true);
          if ((!HOSHD_Spruenge.empty()) && (ColorAndOrder.HOSHD_Showflag))
            plotSpruenge(HOSHD_Spruenge, RGB2Gdk(ColorAndOrder.HOSHD_R,ColorAndOrder.HOSHD_G,ColorAndOrder.HOSHD_B), true);
          break;
        case 1:
          if ((!HOHD_Spruenge.empty()) && (ColorAndOrder.HOHD_Showflag))
            plotSpruenge(HOHD_Spruenge, RGB2Gdk(ColorAndOrder.HOHD_R,ColorAndOrder.HOHD_G,ColorAndOrder.HOHD_B), true);
          if ((!HOSHD_Spruenge.empty()) && (ColorAndOrder.HOSHD_Showflag))
            plotSpruenge(HOSHD_Spruenge, RGB2Gdk(ColorAndOrder.HOSHD_R,ColorAndOrder.HOSHD_G,ColorAndOrder.HOSHD_B), true);
          if ((!SHD_Spruenge.empty()) && (ColorAndOrder.SHD_Showflag))
            plotSpruenge(SHD_Spruenge, RGB2Gdk(ColorAndOrder.SHD_R,ColorAndOrder.SHD_G,ColorAndOrder.SHD_B), true);
          break;
        case 2:
          if ((!SHD_Spruenge.empty()) && (ColorAndOrder.SHD_Showflag))
            plotSpruenge(SHD_Spruenge, RGB2Gdk(ColorAndOrder.SHD_R,ColorAndOrder.SHD_G,ColorAndOrder.SHD_B), true);
          if ((!HOHD_Spruenge.empty()) && (ColorAndOrder.HOHD_Showflag))
            plotSpruenge(HOHD_Spruenge, RGB2Gdk(ColorAndOrder.HOHD_R,ColorAndOrder.HOHD_G,ColorAndOrder.HOHD_B), true);
          if ((!HOSHD_Spruenge.empty()) && (ColorAndOrder.HOSHD_Showflag))
            plotSpruenge(HOSHD_Spruenge, RGB2Gdk(ColorAndOrder.HOSHD_R,ColorAndOrder.HOSHD_G,ColorAndOrder.HOSHD_B), true);
          break;
        case 3:
          if ((!SHD_Spruenge.empty()) && (ColorAndOrder.SHD_Showflag))
            plotSpruenge(SHD_Spruenge, RGB2Gdk(ColorAndOrder.SHD_R,ColorAndOrder.SHD_G,ColorAndOrder.SHD_B), true);
          if ((!HOSHD_Spruenge.empty()) && (ColorAndOrder.HOSHD_Showflag))
            plotSpruenge(HOSHD_Spruenge, RGB2Gdk(ColorAndOrder.HOSHD_R,ColorAndOrder.HOSHD_G,ColorAndOrder.HOSHD_B), true);
          if ((!HOHD_Spruenge.empty()) && (ColorAndOrder.HOHD_Showflag))
            plotSpruenge(HOHD_Spruenge, RGB2Gdk(ColorAndOrder.HOHD_R,ColorAndOrder.HOHD_G,ColorAndOrder.HOHD_B), true);
          break;
        case 4:
          if ((!HOSHD_Spruenge.empty()) && (ColorAndOrder.HOSHD_Showflag))
            plotSpruenge(HOSHD_Spruenge, RGB2Gdk(ColorAndOrder.HOSHD_R,ColorAndOrder.HOSHD_G,ColorAndOrder.HOSHD_B), true);
          if ((!HOHD_Spruenge.empty()) && (ColorAndOrder.HOHD_Showflag))
            plotSpruenge(HOHD_Spruenge, RGB2Gdk(ColorAndOrder.HOHD_R,ColorAndOrder.HOHD_G,ColorAndOrder.HOHD_B), true);
          if ((!SHD_Spruenge.empty()) && (ColorAndOrder.SHD_Showflag))
            plotSpruenge(SHD_Spruenge, RGB2Gdk(ColorAndOrder.SHD_R,ColorAndOrder.SHD_G,ColorAndOrder.SHD_B), true);
          break;
        case 5:
          if ((!HOSHD_Spruenge.empty()) && (ColorAndOrder.HOSHD_Showflag))
            plotSpruenge(HOSHD_Spruenge, RGB2Gdk(ColorAndOrder.HOSHD_R,ColorAndOrder.HOSHD_G,ColorAndOrder.HOSHD_B), true);
          if ((!SHD_Spruenge.empty()) && (ColorAndOrder.SHD_Showflag))
            plotSpruenge(SHD_Spruenge, RGB2Gdk(ColorAndOrder.SHD_R,ColorAndOrder.SHD_G,ColorAndOrder.SHD_B), true);
          if ((!HOHD_Spruenge.empty()) && (ColorAndOrder.HOHD_Showflag))
            plotSpruenge(HOHD_Spruenge, RGB2Gdk(ColorAndOrder.HOHD_R,ColorAndOrder.HOHD_G,ColorAndOrder.HOHD_B), true);
          break;  
	  }
      break;

    case multi_time:
	  datatype = multi;
      short i;
      bool clear;
        
      // daten sichern
      Settings_back = Settings;
      Tampl_histo_type  ahisto_back;
      for (i=0;i < 4096;i++) ahisto_back[i] = ahisto[i];
     
      // zeichnen
      multi_times = Settings.SamplesToProcess / Settings.multi_samples + 1;
      Settings.Zoom = true;
      Settings.enlarge_current = false;
      Settings.SamplesToProcess = Settings.multi_samples;  
      clear = true;
      for (i=0;i < multi_times;i++)
      {
        str_filter(clear); 
        Settings.SamplesToSkip += Settings.multi_samples;
        if (Settings.SamplesToSkip + Settings.multi_samples > Anz_Samples) break;
        clear = false;
      }

      // daten restaurieren
      Settings = Settings_back;
      for (i=0;i < 4096;i++) ahisto[i] = ahisto_back[i];
      break;

    case nothing:
    case growing_ampl_histo:
    case dwell_time:
    case ampl_histo:
    case start_ampl_histo:
	case dwell_2d:
	case dwell_diff:
    default:
      break;
  }
  */
}

void
TDaten::plot_init()
{
double nenner;

  // TZeitDaten
  nenner = g_maxx * Settings.StackedWindows();
  if (Settings.AvarageFilter() == normal) nenner *= Settings.LengthOfAvarageFilter();
  // bei gleitend werden es nicht weniger pixel!
  Samples_per_Pixel = Settings.SamplesToProcess / nenner;
}

void
TDaten::plot_start()
{
int  i;
int  y;

  // TZeitDaten
  if ((Settings.enlarge_current) && (Settings.Zoom))
  {  
    strom_min_old = strom_min;
    strom_max_old = strom_max;
    strom_min = 4095;
    strom_max = 0;
    for (i=Settings.SamplesToSkip;i < Settings.SamplesToProcess+Settings.SamplesToSkip+1;i++)
    {
      if (i>= Anz_Samples) break;
      y=Wert[i];
      if (y > strom_max) strom_max=y;
      if (y < strom_min) strom_min=y;
    }
  }
}

void
TDaten::plot_end()
{
int  a;
int  b;

  // TZeitDaten
  if ((Settings.enlarge_current) && (Settings.Zoom))
  {
    a=strom_min_old;
    b=strom_max_old;
    strom_min_old = strom_min;
    strom_max_old = strom_max;
    strom_min = a; 
    strom_max = b;
	// jetzt stehen strom_min, strom_max fuer enlarge also in old!
	// warum auch immer! strom_min ist jetzt wieder = settings.currentfrom...
  }
}


bool
TDaten::ProcessedText(int i)
{
// drawing text in statusbar ausgeben
bool back = false;

  i = 100 * i / (Settings.StackedWindows() * g_maxx);
  if (i != 0)
  {
    char buf[MaxTextLen];
    sprintf(buf,"drawing... %i%%", i);
	Application.statustext(buf); 
  }
  return back;
}

void
TDaten::str_filter(bool clear)
{
/*2019	
// hier werden punkte in datapix gezeichnet
int  i,j,k,l,x;

const int MAXPT=4096;  // kompromiss zeit gegen kB     
GdkPoint point[MAXPT];
int     pointnr=0; 
bool *nopixel=NULL;	// auf array (stackedwindow spalte) - pixel schon gemalt?

double awert = 0;
double rel_height;
unsigned short y_pixel;
short  window_height;
//
 short  rhythm;
double t1, t2;
int  t1r, t2r;

short vor;
short nach;
short laenge;

  plot_init();
  plot_start();

  if (clear)
  {
    // background loeschen
    gdk_draw_rectangle(DataPix, Application.back_gc, true,
                        0, 0, g_maxx, g_maxy);
  }
  for (i=0;i < 4096;i++) ahisto[i]=0;

  window_height = g_maxy / Settings.StackedWindows();

  rhythm = g_maxx / 8;
  x = Settings.SamplesToSkip;

  nopixel = new bool[window_height+1];   // array fuer stackedwindows spalte erzeugen
 
  // wie weit schaut der gleitende filter vor und zurueck?
  vor = Settings.LengthOfAvarageFilter()/2;
  if (Settings.LengthOfAvarageFilter() % 2 == 0)
  {
    nach = vor;
  } else {
    nach = vor + 1;
  }

  for (i=Settings.StackedWindows()-1;i >= 0;i--)
  {
    for (j=0;j < g_maxx;j++)
    {
      for (k=0;k < window_height;k++) nopixel[k]=true; // noch kein pixel in dieser spalte
 
      if (j % rhythm == 0)
	  {
        if (ProcessedText(j + (Settings.StackedWindows()-i-1)*g_maxx)) 
            return;
           
	  }
      t1 = j * Samples_per_Pixel + (Settings.StackedWindows()-i-1) * Samples_per_Pixel * g_maxx;
      t2 = t1 + Samples_per_Pixel;  // *1 Sample
      t1r = round(t1);
      t2r = round(t2) - 1;

      for (k=t1r;k < t2r+1;k++)
      {
        if (x >= (Anz_Samples-Settings.LengthOfAvarageFilter())) break;

		switch (Settings.AvarageFilter())
		{
          case normal:
		    awert = 0;
		    for (l=0;l < Settings.LengthOfAvarageFilter();l++)
			{
              awert += Wert[x];
              x++;
			}
            awert /= Settings.LengthOfAvarageFilter();
            break;
          case keine:
            awert = Wert[x];
			x++;
            break;
          case gleitend:      
     	    awert = 0;
            for (l=x-vor;l < x+nach;l++)
			{
              if ((l >= 0) && (l < Anz_Samples)) awert += Wert[l];
			}
		    laenge = Settings.LengthOfAvarageFilter();
            if (x - vor < 0) laenge += (x - vor); // x-vor ist negativ!
            if (x + nach >= Daten.Anz_Samples) laenge -= (x+nach - Daten.Anz_Samples +1);
		    awert /= laenge;
			x++;
		}	
		// hier kann es zu Laufzeitfehlern kommen (PASCAL), wenn bereich zu klein
        // wennn Division durch 1 statt durch 0, fehler spaeter!
        rel_height=(awert-strom_min)/(strom_max-strom_min);

        y_pixel=(short)(rel_height*window_height) & 0x03FF;

        // nur pixel, die noch nicht in liste sind aufnehmen    
        if (y_pixel <= window_height) //unsigned int Y-pixel immder >=0
        {
          if (nopixel[y_pixel])
          {
            point[pointnr].x = j;
            point[pointnr].y = g_maxy - y_pixel - (i*window_height);
            nopixel[y_pixel] = false;
            pointnr++;
            if (pointnr >= MAXPT)
            {
              gdk_draw_points(DataPix, Application.dots_gc,point,MAXPT);
              pointnr = 0;
            }
          }
        }
        ahisto[(int)round(awert)]++;
      }
    }
  }
  if (pointnr!=0) gdk_draw_points(DataPix, Application.dots_gc,point,pointnr);
  delete[] nopixel; // dynamisch erzeugtes array loeschen!
  plot_end();
2019*/  
}

void
TDaten::histogram(Tampl_histo_type a)
{
// start amplituden histogramm wird berechnet
// (alle datenpunkte ohne filter)
int  i;
int awert;
int  rhythm;

  // THistoDaten
  rhythm = Anz_Samples / 100;
  
  for (i=0;i < 65536;i++) a[i]=0;

  for (i=0;i < Anz_Samples;i++)
  {
    awert = Wert[i];
    if (i % rhythm == 0)
	{
      char buf[MaxTextLen];
      sprintf(buf,"samples: %i", i);
	  Application.statustext(buf);
	}
	/*Tobias 2022	
    if ((awert & 0xF000) != 0)
    {
    // kann eigentlich nicht mehr sein!, 
	// marker wurden aus oberen bits extrahiert und diese dann zurueckgesetzt
    warning((char *)"corrupt data");
    std::cout<<"i: "<<i<<"Wert: "<<awert<<std::endl;
    }
    awert &= 0x0FFF;
      */ 
    Wert[i] = awert;
    a[awert]++;
   }  
}


void
TDaten::min_max()
{
  // THistoDaten
  histogram(bhisto);
  ampl_min_max(bhisto, strom_min, strom_max);
}

void
TDaten::ampl_min_max(Tampl_histo_type a, int &min, int &max)          
{
// minimaler und maximaler strom werden berechnet
int i;
int breite;
int one_percent;
int sum;

  // THistoDaten
  min = 0;
  while (a[min]==0)
  {
    min++;
    if (min==65536) break;
  }
  max = 65535;
  while (a[max]==0)
  {
    max--;
    if (max==-1) break;
  };
  if (min >= max)
  {
    min=0;
    max=0;
    return;
  }
  sum=0;
  for (i=min;i<max+1;i++) sum+=a[i];
  one_percent = sum / 100;
  sum=0;
  while (sum < one_percent)
  {
    if (min >=65536) break;
    sum+=a[min];
    min++;
  }
  min--;
  sum=0;
  while (sum < one_percent)
  {
    if (max <= 0) break;
    sum+=a[max];
    max--;
  }
  max++;
  while ((min>0) && (a[min] != 0)) min--;
  while ((max<65536) && (a[max] != 0)) max++;
  if (max < min + 50)
    if (max < 2000) {max+=50;} else {min-=50;}
  breite= max-min+1;
  min = min - breite / 4;
  max = max + breite / 4;
  if (min < 0) min =0;
  if (max > 65535) max = 65535;
  Settings.CurrentFrom = min;
  Settings.CurrentTo = max;

//BEGIN KARSTEN
Minimaler_Strom=min;
Maximaler_Strom=max;
//END KARSTEN

}
/*2019
void
TDaten::plot_strom_height(double strom,
                          double s_min, double s_max,
                          double height, double max_h,
                          bool aufdruecken, bool logscaling, 
                          short   maxx, short maxy,
                          GdkGC *gc)																	//Tobias 2017 weg: = Application.dots_gc
{
// ausgabe einer linie vom letzten punkt zu strom/height
// wenn aufdruecken == false nur setzten der alten koordinaten
short  x_pixel;
short  y_pixel;
short  xrand = 0;
short  yrand = 0;

static short x_old;
static short y_old;
double rel_height;


  switch (graphtype)
  {
    case ampl_histo:
    case start_ampl_histo:
      // fuer achsen 
	  //maxx -= rand_pixel;
      //maxy -= rand_pixel; // keine referenz, also nur hier!
      //xrand = rand_pixel;

	  // fuer box auch oben rand
	  maxx -= 2*rand_pixel;
      maxy -= 2*rand_pixel;
      xrand = yrand = rand_pixel;
      break;
	default:
	  break;
  }

  // THistoDaten
  rel_height = (strom - s_min)/(s_max - s_min);
  x_pixel = (short)(rel_height * maxx) + xrand;
  
  if (logscaling)
  {
	if (height > 0)
	{
	  y_pixel = (short)(log10(height) / log10(max_h) * maxy);
	} else {
	  y_pixel = 0;
	}
  } else {
    y_pixel = (short)(height / max_h * maxy);
  }
  if (y_pixel < -rand_pixel / 3) y_pixel = 0; 
  if (y_pixel > maxy) y_pixel = maxy;

  y_pixel -= yrand; //oberen rahmen bei box beruecksichtigen

  if (aufdruecken)
  {
    gdk_draw_line(Daten.MemPix, gc,
                  x_old, maxy - y_old,
                  x_pixel, maxy - y_pixel);
  }
  x_old = x_pixel;
  y_old = y_pixel;
}
2019*/

/*GdkPoint
TDaten::strom_point(double strom,  
                    double s_min, double s_max,
    			    double height, double max_h,
                    bool logscaling, 
                    short maxx, short maxy)
{
// liefert position des zu malenden Punktes strom/height
// fuer gdk_draw_points, ..._lines, ...-segments
GdkPoint point;
short  xrand = 0;
short  yrand = 0;

double rel_height;

  switch (graphtype)
  {
    case ampl_histo:
    case start_ampl_histo:
      // fuer achsen 
	  //maxx -= rand_pixel;
      //maxy -= rand_pixel; // keine referenz, also nur hier!
      //xrand = rand_pixel;

	  // fuer box auch oben rand
	  maxx -= 2*rand_pixel;
      maxy -= 2*rand_pixel;
      xrand = yrand = rand_pixel;
      break;
	default:
	  break;
  }

  // THistoDaten
  rel_height = (strom - s_min)/(s_max - s_min);
  point.x = (short)(rel_height * maxx) + xrand;
  
  if (logscaling)
  {
	if (height > 0)
	{
	  point.y = (short)(log10(height) / log10(max_h) * maxy);
	} else {
	  point.y = 0;
	}
  } else {
    point.y = (short)(height / max_h * maxy);
  }
  if (point.y < -rand_pixel / 3) point.y = 0; 
  if (point.y > maxy) point.y = maxy;          

  point.y -= yrand; //oberen rahmen bei box beruecksichtigen

  point.y = maxy - point.y;

  return point;
}*/

void
TDaten::copy_stromsig()
{
/*2019	
int  i;
int  grow_rhythm;
int  maxheight;
unsigned short awert;
Tampl_histo_type  growhisto;

  // THistoDaten
  maxheight=1;
  for (i=strom_min;i < strom_max+1;i++)
  {
    if (maxheight < bhisto[i]) maxheight = bhisto[i];
  }
  grow_rhythm = Settings.SamplesToProcess / 15;
  for (i=0;i < 4096;i++) growhisto[i] = 0;
  for (i=0;i < Settings.SamplesToProcess;i++)
  {
    if (i+Settings.SamplesToSkip > Anz_Samples) break;
    awert=Wert[i+Settings.SamplesToSkip];
    if ((awert & 0xF000) != 0) 
      {
       warning((char *)"corrupt data!");
       std::cout<<"i: "<<i<<"Wert: "<<awert<<std::endl;
      } 
      
    awert &= 0x0FFF;
    Wert[i+Settings.SamplesToSkip] = awert;
    growhisto[awert]++;
    if (i % grow_rhythm == 0) amp_grow(growhisto, maxheight);
  }
  amp_grow(growhisto, maxheight);

  // levels zeichnen (nicht logscale, immer von oben nach unten)
  // sollte eigentlich ausgelagert werden, zusammen mit dem teil aus
  // plot_levels_and_fit.
  // dann nach plot levels -> mit flag, ob zeitreihe oder histogramm
  for (int level=0;level < a_fit.n_channels+1;level++)
  {
    for (int sublevel=0;sublevel < a_fit.n_extra_channels+1;sublevel++)
    {
      plot_strom_height(a_fit.i_null + level*a_fit.i_channel + sublevel*a_fit.i_extra_channel,
                        strom_min,strom_max,
                        maxheight, maxheight, false, false, g_maxx, g_maxy);
      plot_strom_height(a_fit.i_null + level*a_fit.i_channel + sublevel*a_fit.i_extra_channel,
                        strom_min,strom_max,
                        0, maxheight, true, false, g_maxx, g_maxy,
                        Application.current_gc);
    }
  }

2019*/
}

void
TDaten::amp_grow(Tampl_histo_type a, int max)
{
/*2019	
short  i;
short  j;
double aa;
// zeichnet growing amplituden histogramm

  // THistoDaten
  plot_strom_height(0, 0, 4096, 0, max, false, false, g_maxx, g_maxy);
  aa = 0;
  j = 0;
  for (i=strom_min;i < strom_max+1;i++)
  {
    aa += a[i];
    j++;
    if (j % 20 == 0)
    {
      plot_strom_height(i, strom_min, strom_max, aa/20, max, true, logscalehisto, g_maxx, g_maxy, Application.hisfit_gc);
      aa = 0;
    }
  }
2019*/  
}

void
TDaten::plot_levels()
{
int   i;
int   j;
// zeichnet levels in zeitreihe
  for (i=0;i<a_fit.n_channels+1;i++)
    for (j=0;j<a_fit.n_extra_channels+1;j++)
      plot_strom(a_fit.i_null + i*a_fit.i_channel + j*a_fit.i_extra_channel);
}

void
TDaten::plot_strom(double strom)
{
/*2019	
short  i;
short  window_height;
short  y_pixel;
double rel_height;

  // solange keine 1bit pixmaps fuer spruenge
  // level mit XOR �ber die spruenge malen
  if(!schwarzaufweiss) gdk_gc_set_function (Application.current_gc, GDK_XOR);

  window_height = g_maxy / Settings.StackedWindows();
  if ((Settings.enlarge_current) && (Settings.Zoom))
  {
     rel_height = (double)(strom-strom_min_old)/(strom_max_old-strom_min_old);
  } else {
     rel_height = (double)(strom-strom_min)/(strom_max-strom_min);
  }
  for (i=0;i<Settings.StackedWindows();i++)
  {
    y_pixel=(short)(rel_height * window_height);
    gdk_draw_line(MemPix ,Application.current_gc,
                  0, g_maxy - y_pixel -(i * window_height),
                  g_maxx, g_maxy - y_pixel -(i * window_height));
  }
  gdk_gc_set_function (Application.current_gc,GDK_COPY);
2019*/  
}

bool
TDaten::checkampl(Tampl_fit_type ampl)
{
// prueft, ob alle werte innerhalb des erlaubten bereichs liegen
// besser bei Tampl_fit_type aufgehoben
int min, max;

  if ((Settings.enlarge_current) && (Settings.Zoom))
  {  
    min = strom_min_old;
    max = strom_max_old;
  } else {
    min = Settings.CurrentFrom;
	max = Settings.CurrentTo;
  }

  return
     // baseline innerhalb des stackedwindows
	((ampl.i_null > min) && 
     (ampl.i_null < max) &&
     // great channels innerhalb
     (ampl.i_null + ampl.n_channels * ampl.i_channel > min) &&
     (ampl.i_null + ampl.n_channels * ampl.i_channel < max) &&
     // small channels innerhalb
     (ampl.i_null + ampl.n_channels * ampl.i_channel + ampl.n_extra_channels * ampl.i_extra_channel > min) &&
     (ampl.i_null + ampl.n_channels * ampl.i_channel + ampl.n_extra_channels * ampl.i_extra_channel < max) &&
     (ampl.i_null + ampl.n_extra_channels * ampl.i_extra_channel) > min &&
     (ampl.i_null + ampl.n_extra_channels * ampl.i_extra_channel) < max &&
     // anzahl channels gueltig
     (ampl.n_channels >= 0) && (ampl.n_channels <= n_channels_max) && 
     (ampl.n_extra_channels >= 0) && 
	 (ampl.n_extra_channels <= n_extra_channels_max));
}

void
TDaten::updown(int w)
{
// andern eines wertes von a_fit_
Tampl_fit_type ampl = a_fit; // so kann ich testen, bevor ich aendere

  switch (command)
  {
  case cnone:
    return;
    break;
  case cbaseline: 
    ampl.i_null += w;
    if(checkampl(ampl)) a_fit.i_null += w;
    break;  
  case cccurrent:
    ampl.i_channel += w;
    if(checkampl(ampl)) a_fit.i_channel += w;
    break;
  case cucurrent: 
    ampl.i_extra_channel += w;
    if(checkampl(ampl)) a_fit.i_extra_channel += w;
    break;
  case cgreatchannels: 
    ampl.n_channels += w;
    if(checkampl(ampl)) a_fit.n_channels += w;
    break;
  case csmallchannels:
    ampl.n_extra_channels += w;
    if(checkampl(ampl)) a_fit.n_extra_channels += w;
    break;
  case csigma: 
    if ((a_fit.sigma + w >= 4.0) && (a_fit.sigma + w <= 4095.0)) 
	{  
	  a_fit.sigma = a_fit.sigma + w;
	}
    break;  
  }
}

void
TDaten::plot_text()
{
// text in statusbar ausgeben
char textbuf[MaxTextLen];
char fehlerbuf[MaxTextLen];
char filterbuf[MaxTextLen];

  // fehlertext in statusbar?
  switch (graphtype)
  {
    case start_ampl_histo:
    case ampl_histo:
      sprintf(fehlerbuf, "     E: %.2f", fehler);
      break;
    default:
      fehlerbuf[0] = 0; 
      break;
  }

  // filterlaengetext fuer statusbar
  switch (graphtype)
  {
    case start_ampl_histo:
      sprintf(filterbuf, "FL: 1");
      break;
    case growing_ampl_histo:
    case ampl_histo:
    case time_series:
    case dwell_time:
    case multi_time:
  	  if (Settings.AvarageFilter() == gleitend)
	  {
        sprintf(filterbuf, "mFL: %i", Settings.LengthOfAvarageFilter());
	  } else {
        sprintf(filterbuf, "FL: %i", Settings.LengthOfAvarageFilter());
	  }
      break;
    case dwell_2d:
    case nothing:
    default:
	  filterbuf[0] = 0;
	  break;
  }   

  // statusbar basteln 
  switch (graphtype)
  {
    case growing_ampl_histo:
    case start_ampl_histo:
    case ampl_histo:
    case time_series:
    case dwell_time:
    case multi_time:
   								
/*		sprintf(textbuf,"N: %i     C: %i    %.2fpA     U: %d   %.2fpA     G: %d  K: %d     S: %d   %.2fpA    %s%s",
              (int)a_fit.i_null,                // Nullstrom in A/D-Units ausgeben
              (int)a_fit.i_channel,             // Kanalstrom in A/D-Units ausgeben
              a_fit.i_channel*Settings.VFactor/Settings.Gain, // Kanalstrom in pA ausgeben
              (int)a_fit.i_extra_channel,       // Sub-Level-Strom in A/D-Units ausgeben
              a_fit.i_extra_channel*Settings.VFactor/Settings.Gain, // Sub-Level-Strom in pA
              (int)a_fit.n_channels,            // Zahl der Niveaus ausgeben
              (int)a_fit.n_extra_channels,      // Zahl der Sub-Niveaus ausgeben
              (int)a_fit.sigma,                 // Standardabweichung in A/D-Units ausgeben
              a_fit.sigma*Settings.VFactor/Settings.Gain, // Standardabweichung in pA   
              filterbuf,                          // Filterl�nge ausgeben
              fehlerbuf);                         // fehler oder nichts (s.o.)    
*/    
	  break;
    case dwell_2d:
      if (Dwell_2d_ptr->sqrt_e())
	  {
        // wurzel aus events bringt nichtlineare stauchung
        sprintf(textbuf,"%s: C: %i - %i     O: %i - %i     sE: %.1f - %.1f     S: %.1f     BpL: %i",
              Dwell_2d_ptr->name(),
 		      Dwell_2d_ptr->min_close_log(), 
		      Dwell_2d_ptr->max_close_log(),
			  Dwell_2d_ptr->min_open_log(),
			  Dwell_2d_ptr->max_open_log(),
              sqrt(Dwell_2d_ptr->min_e(false)),
              sqrt(Dwell_2d_ptr->max_e(false)),
              Dwell_2d_ptr->events(),
	  		  Dwell_2d_ptr->bins_per_log());
	  } else {
        sprintf(textbuf,"%s: C: %i - %i     O: %i - %i     E: %.1f - %.1f     S: %.1f     BpL: %i",
              Dwell_2d_ptr->name(),
		      Dwell_2d_ptr->min_close_log(), 
		      Dwell_2d_ptr->max_close_log(),
			  Dwell_2d_ptr->min_open_log(),
			  Dwell_2d_ptr->max_open_log(),
              Dwell_2d_ptr->min_e(false),
              Dwell_2d_ptr->max_e(false),
              Dwell_2d_ptr->events(),
			  Dwell_2d_ptr->bins_per_log());
	  }
      break;
	case dwell_diff:
      // wurzel hier nicht moeglich, negative werte!
	  // und loglikelihood ausgeben
      sprintf(textbuf,"%s: C: %i - %i     O: %i - %i     D: %.1f - %.1f     BpL: %i     QA: %.1f (%.1f%%)",
              Dwell_2d_diff.name(),
		      Dwell_2d_diff.min_close_log(),
		      Dwell_2d_diff.max_close_log(),
			  Dwell_2d_diff.min_open_log(),
			  Dwell_2d_diff.max_open_log(),
              Dwell_2d_diff.min_e(true),
              Dwell_2d_diff.max_e(true),
			  Dwell_2d_diff.bins_per_log(),
              Dwell_2d_A.lnlikelihood(Dwell_2d_B),
			  Dwell_2d_A.checkrange(Dwell_2d_B) / Dwell_2d_A.events() * 100); // % in bereich der simulierten daten
	  break;
    case nothing:
    default:
   	  textbuf[0] = 0;
      break;
  }
  Application.statustext(textbuf);
}

void 
TDaten::show_ampl_histo(Tampl_histo_type a) 
{
/*2019	
// hier werden rahmen und daten (amplituden histogramm) gezeichnet
int i;

  max_height = 1;
  for (i = strom_min; i < strom_max+1;i++)
  {
    if (max_height < a[i]) max_height = a[i];
  }

  // speicher fuer liste der punkte reservieren
  GdkPoint* points = new GdkPoint[strom_max + 1 - strom_min];

  for (i=strom_min; i < strom_max+1;i++)
  {
    points[i - strom_min] = strom_point(i, 
	                                    strom_min,strom_max, 
			                            a[i], max_height, 
			                            logscalehisto, 
	     	                            g_maxx, g_maxy);
  }
  gdk_draw_lines(MemPix, Application.dots_gc, points, strom_max + 1 - strom_min);

  //speicher wieder freigeben
  delete[] points;

  int j;
  short stromrand = (short)((double)(strom_max - strom_min) / (g_maxx - rand_pixel) * rand_pixel);// / 3);
  short yticslen = stromrand;
  //short heightrand = (short)((double)max_height / (g_maxy - rand_pixel) * rand_pixel / 3);

//
  // Achsen mit Pfeilen
//  gdk_draw_line(Daten.MemPix, Application.border_gc,
//                rand_pixel, g_maxy - rand_pixel,
//                rand_pixel, 0);
//  gdk_draw_line(Daten.MemPix, Application.border_gc,
//                rand_pixel, 0,
//	            rand_pixel * 2/3, rand_pixel / 2);
//  gdk_draw_line(Daten.MemPix, Application.border_gc,
//  	            rand_pixel, 0,
//				rand_pixel * 4/3, rand_pixel / 2);
// 
//  gdk_draw_line(Daten.MemPix, Application.border_gc,
//                rand_pixel, g_maxy - rand_pixel,
//                g_maxx, g_maxy - rand_pixel);
//  gdk_draw_line(Daten.MemPix, Application.border_gc,
//                g_maxx, g_maxy - rand_pixel,
//                g_maxx - rand_pixel / 2, g_maxy - rand_pixel * 2/3);
//  gdk_draw_line(Daten.MemPix, Application.border_gc,
//                g_maxx, g_maxy - rand_pixel,
//                g_maxx - rand_pixel / 2, g_maxy - rand_pixel * 4/3);

  // box
  gdk_draw_line(Daten.MemPix, Application.border_gc,
                rand_pixel, g_maxy - rand_pixel,
                rand_pixel, rand_pixel);
  gdk_draw_line(Daten.MemPix, Application.border_gc,
                rand_pixel, g_maxy - rand_pixel,
                g_maxx - rand_pixel, g_maxy - rand_pixel);
  gdk_draw_line(Daten.MemPix, Application.border_gc,
                g_maxx - rand_pixel, g_maxy - rand_pixel,
                g_maxx - rand_pixel, rand_pixel);
  gdk_draw_line(Daten.MemPix, Application.border_gc,
                rand_pixel, rand_pixel,
                g_maxx - rand_pixel, rand_pixel);


  // y-achsen marker 
  if (logscalehisto)
  {
    i = 1; j = 1;
    while ((i*j <= max_height) && (i*j > 0))
	{
      for (i=1;i <= 9;i++)
	  {
		if (i*j <= max_height)
		{  
		  if (i==1)
		  { 
		    yticslen = 2 * stromrand;
		  } else {
		    yticslen = stromrand;
		  }
		  // fuer achsen
          //plot_strom_height(strom_min, strom_min, strom_max, i * j, max_height, false, logscalehisto, g_maxx, g_maxy);  
          //plot_strom_height(strom_min - yticslen, strom_min, strom_max, i * j, max_height, true, logscalehisto, g_maxx, g_maxy, Application.border_gc);  
		  // fuer box
		  plot_strom_height(strom_min, strom_min, strom_max, i * j, max_height, false, logscalehisto, g_maxx, g_maxy);  
		  plot_strom_height(strom_min + yticslen, strom_min, strom_max, i * j, max_height, true, logscalehisto, g_maxx, g_maxy, Application.border_gc);
          plot_strom_height(strom_max, strom_min, strom_max, i * j, max_height, false, logscalehisto, g_maxx, g_maxy, Application.border_gc);
          plot_strom_height(strom_max - yticslen, strom_min, strom_max, i * j, max_height, true, logscalehisto, g_maxx, g_maxy, Application.border_gc);
		};
	  }
      i=1;j *= 10;
	}
  } else {
    i = 1;
    for (j=1;j <= (int)log10(max_height);j++) i *= 10;
	if (i*2 > max_height) i /= 10;
    j = i;
    while (j <= max_height)
	{
      // fuer achsen
	  //plot_strom_height(strom_min, strom_min, strom_max, j, max_height, false, logscalehisto, g_maxx, g_maxy);  
      //plot_strom_height(strom_min - stromrand, strom_min, strom_max, j, max_height, true, logscalehisto, g_maxx, g_maxy, Application.border_gc);  
	  // fuer box
      plot_strom_height(strom_min, strom_min, strom_max, j, max_height, false, logscalehisto, g_maxx, g_maxy);  
	  plot_strom_height(strom_min + stromrand, strom_min, strom_max, j, max_height, true, logscalehisto, g_maxx, g_maxy, Application.border_gc);
      plot_strom_height(strom_max, strom_min, strom_max, j, max_height, false, logscalehisto, g_maxx, g_maxy, Application.border_gc);
      plot_strom_height(strom_max - stromrand, strom_min, strom_max, j, max_height, true, logscalehisto, g_maxx, g_maxy, Application.border_gc);
      j += i;
	}
  }
//
//  // x-achsen marker
//  i = 1;
//  for (j=1;j <= (int)log10(strom_max - strom_min);j++) i *= 10;
//  if (i*2 > (strom_max - strom_min)) i /= 10; 
//  j = i;
//  while (j <= strom_max)
//  {
//    if (j >= strom_min)
//	{
//      plot_strom_height(j, strom_min, strom_max, -heightrand, max_height, false, false, g_maxx, g_maxy);  
//      plot_strom_height(j, strom_min, strom_max, 0, max_height, true, false, g_maxx, g_maxy, Application.border_gc);  
//	}
//    j += i;
//  }
//

//
//  // systemfont ermitteln
//  GdkGCValues values;
//  gdk_gc_get_values((Application.DLG_Main)->style->bg_gc[GTK_STATE_NORMAL],&values);
//
//  // Ausgabe der Texte
//  gdk_draw_string(Daten.MemPix,
//			      values.font,
//			      Application.border_gc,
//			      g_maxx - 40,
//			      g_maxy - 2,
//			      "I/pA");
//  gdk_draw_string(Daten.MemPix,
//			      values.font,
//			      Application.border_gc,
//			      2,
//			      22,
//			      "N");
2019*/
}
  
void 
TDaten::calc_fit(Tampl_fit_type &af, Tampl_histo_type histo, int min, int max)
{
int    i,k;
TGaussFit gauss(af);      // Neues Object f�r Linear Least Squares

//double   akt_niveau = 0;    
  
  //cout << "calc_fit\n";
  //cout << gauss.ax[1] << gauss.ay[1] << gauss.asig[1]
  //     << gauss.NumPoints << gauss.aa[1]
  //     << gauss.ama << gauss.alista[1]
  //     << gauss.amfit << gauss.acovar[1][1] << gauss.achisq << endl;

  gauss.NumPoints = strom_max - strom_min;    // Zahl der Datenpunkte initialisieren

  for(i=1;i < gauss.NumPoints+1;i++)
  {
    gauss.ax[i] = i + strom_min;                  // x,
    gauss.ay[i] = histo[i-1+strom_min];           // y,

    if (gauss.ay[i] != 0)
    {
      gauss.asig[i] = sqrt(fabs(gauss.ay[i]));    // und Sigma setzen
    } else {
      gauss.asig[i] = 1;
    }
  }
  k = gauss.determine_base();
  // PASCAL kommentar
  // for (i=1;i < k+1;i++) gauss.alista[i]=i;
  // gauss.ama=k;
  // gauss.amfit=k;

  gauss.lfit(gauss.ax, gauss.ay, gauss.asig, gauss.NumPoints, 
             gauss.aa, gauss.ama, gauss.alista, gauss.amfit, 
             gauss.acovar, gauss.achisq);  // Gauss-Routine aufrufen

  // Koeffizienten der Gauss-Kurven sichen
  for (i=1;i < k+1;i++) af.a[i] = gauss.aa[i];

  fehler = gauss.func(gauss.ap);
}

void 
TDaten::call_simplex()
{
  if (graphtype == ampl_histo)
  {
    simplex(a_fit, ahisto, strom_min, strom_max);
  } else { 
    if (graphtype == start_ampl_histo)
    {
      simplex(a_fit, bhisto, strom_min, strom_max);
    }
  }
}

void 
TDaten::simplex(Tampl_fit_type &af, Tampl_histo_type histo, int min, int max)
{
int     i;//,j,k;
TGaussFit  simpl(af);  // Neues Object f�r Simplex
short     bla;
double    e;

  simpl.NumPoints = strom_max - strom_min;      // Zahl der Datenpunkte initialisieren
  for (i=1;i < simpl.NumPoints+1;i++)
  {
    simpl.ax[i] = i + strom_min;                // x,
    simpl.ay[i] = histo[i-1+strom_min];         // y,
    if (simpl.ay[i] != 0)
    {
      simpl.asig[i] = sqrt(fabs(simpl.ay[i]));  // und Sigma setzen        
    } else {
      simpl.asig[i] = 1;
    }
  }
  e = simpl.fit(0.0000001,bla);    // fitten
  simpl.done(af);
}

void 
TDaten::plot_levels_and_fit(Tampl_fit_type &a, int max_h, double stromrand, double heightrand) //Tobias 2017 weg: =0
{
/*2019	
short level = 0, sublevel = 0, int_strom = 0, i = 0;
short anz_pixel = 0;
double   height;
TGaussFit gauss(a);

  // zeichenposition zuruecksetzten
  plot_strom_height(strom_min, strom_min,strom_max, 0, max_h, false, logscalehisto, g_maxx, g_maxy);

  // levels zeichnen (nicht logscale, immer von oben nach unten)
  // sollte eigentlich ausgelagert werden, zusammen mit dem teil aus
  // copy_stromsig.
  // dann nach plot levels -> mit flag, ob zeitreihe oder histogramm
  for (level=0;level < a.n_channels+1;level++)
  {
    for (sublevel=0;sublevel < a.n_extra_channels+1;sublevel++)
    {
      plot_strom_height(a.i_null + level*a.i_channel + sublevel*a.i_extra_channel,
                        strom_min - stromrand,strom_max,
                        max_h, max_h, false, false, g_maxx, g_maxy);
      plot_strom_height(a.i_null + level*a.i_channel + sublevel*a.i_extra_channel,
                        strom_min - stromrand,strom_max,
                        0, max_h, true, false, g_maxx, g_maxy,
                        Application.current_gc);
    }
  }

  for (i=1;i < map+1;i++) gauss.aa[i] = a.a[i];

  // fuer summe und jede glocke
  std::vector<GdkPoint*> points;
  points.resize((gauss.niveaus.size()+1));
  
  // brauche ich dynamischen speicher zum ablegen der stuetzpunkte
  anz_pixel = (strom_max + 1 - strom_min);
  for (unsigned short ii=0;ii < gauss.niveaus.size()+1;ii++)
  {
   points[ii] =  new GdkPoint[anz_pixel];
  }

  RealArrayMA  afunc;
  for (int_strom=strom_min;int_strom < strom_max+1;int_strom++)
  {
    gauss.funcs(int_strom, afunc);

    // summe der gauss glocken
    height = gauss.summe_gauss_glocken(afunc);
    points[0][int_strom - strom_min] = strom_point(int_strom, 
	                                               strom_min,strom_max, 
				                                   height, max_h, 
				                                   logscalehisto, 
				                                   g_maxx, g_maxy);

	// einzelne gauss glocken
    for (unsigned short glocke=1;glocke < gauss.niveaus.size()+1;glocke++)
	{
      height = gauss.gauss_glocke(glocke, afunc);
      points[glocke][int_strom - strom_min] = strom_point(int_strom, 
                                                   strom_min,strom_max, 
	                                               height, max_h, 
				                                   logscalehisto, 
				                                   g_maxx, g_maxy);
	}
  }
  gdk_draw_lines(MemPix, Application.border_gc, points[0], anz_pixel);

  // zeichnen und dynamischen speicher wieder freigeben
  for (unsigned short ii=0;ii < gauss.niveaus.size()+1;ii++)
  {
    gdk_draw_lines(MemPix, Application.border_gc, points[ii], anz_pixel);
    delete[] points[ii];
  }
  
//for (unsigned short ii=1;ii < gauss.niveaus.size()+1;ii++)
//{
//  cout << gauss.aa[ii] << endl;
//}
2019*/
}

void 
TDaten::calc_ahisto()
{
// wird auch schon in str_filter (PASCAL) berechnet, 
// aber auch durch settings �nderung n�tig!
int  x,k,l;
double awert;
int  vor, nach, laenge;

  // wie weit schaut der gleitende filter vor und zurueck?
  vor = Settings.LengthOfAvarageFilter()/2;
  if (Settings.LengthOfAvarageFilter() % 2 == 0)
  {
    nach = vor;
  } else {
    nach = vor + 1;
  }

  // histo zuruecksetzen
  for (k=0;k < 4096;k++) ahisto[k]=0;

  x = Settings.SamplesToSkip;

  switch (Settings.AvarageFilter())
  {
    case normal:
      for (k=x;k <= x + Settings.SamplesToProcess;k++)
      {
        if (k >= (Anz_Samples-Settings.LengthOfAvarageFilter())) break;
        awert = 0;
        for (l=0;l < Settings.LengthOfAvarageFilter();l++)
        {
          awert += Wert[k];
          k++;
        }
        awert = awert / Settings.LengthOfAvarageFilter();
             
        ahisto[(int)round(awert)]++;
      }
      break;
   
    case keine:
      for (k=x;k <= x + Settings.SamplesToProcess;k++)
      {
        if (k >= Anz_Samples) break;
        ahisto[Wert[k]]++;
      }
      break;
    case gleitend:
      for (k=x;k <= x + Settings.SamplesToProcess;k++)
      {
        if (k >= (Anz_Samples-vor)) break;
  	    awert = 0;

		for (l=x-vor;l < x+nach;l++)
		{
          if ((l >= 0) && (l < Anz_Samples)) awert += Wert[l];
		}

	    laenge = Settings.LengthOfAvarageFilter();
        if (x - vor < 0) laenge += (x - vor); // x-vor ist negativ!
        if (x + nach >= Daten.Anz_Samples) laenge -= (x+nach - Daten.Anz_Samples +1);
		awert /= laenge;

        ahisto[(int)round(awert)]++;

		x++;
	  }
      break;
  }
}
