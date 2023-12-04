/*******************************************************************
 *  Kiel-Patch
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
// #include <process.h>    //2019 
#include <errno.h>
#include <stdio.h>

#include <algorithm>


//#include "gtkapp.h"
//#include "declare.h"
#include "daten.h"
#include "detector.h"
#include "fit.h"
#include "error.h"
#include "round.h"
//#include "gtkfunc.h"
#include "matu.h"


#include "testout.h"

void
TDaten::InitDwellTimeHisto(int n_channels, int n_extra_channels, bool neu)
{
  // TDaten
  result.dwell_time_histo_array.clear();
  if (neu)
  {
	result.dwell_time_histo_array.resize(2*(n_extra_channels+1)*(n_channels+1), n_bins_max+1, 0);
  }
}

void
TDaten::hinkley(detector_type d)
{
//  PASCAL testet auf nil nicht auf count=0, nur integritaet?
//  if ((!HOHD_Spruenge.empty()) && (!SHD_Spruenge.empty()) && (!HOSHD_Spruenge.empty()))
  {
    clearTransitions();
    clearDetectedLevels();
    switch (d)
    {
      case _HOHD: 
        HOHD();
      break;
      
      case _SHD:
        SHD();
      break;
      
      case _SHD_AH:
        SHD_AH();
      break;
      
      case _HOSHD: 
        HOSHD();
      break;
    }
  }
}

void
TDaten::HOSHD()
{
std::cout<<"----"<<std::endl;		
int  i;
int  n_open, n_extra_open;
int  Samples;
Thohdetector  *detector;

  InitDwellTimeHisto(a_fit.n_channels, a_fit.n_extra_channels, true);
  Samples = Settings.SamplesToProcess + Settings.SamplesToSkip;
  i = Settings.SamplesToSkip;
  HOSHD_Spruenge.clear();
  n_open = 0;
  n_extra_open = 0;
  while (i < Samples)
  {
    if (i > Anz_Samples) break;

    detector = new Thohdetector(&Daten,
                                i,
                                n_open,n_extra_open,
                                a_fit.i_null,a_fit.i_channel,a_fit.i_extra_channel,
                                Samples);
//    detector^.done;
  }
}

void
TDaten::HOHD()
{	
double  lambda;

  SNR = a_fit.i_channel / a_fit.sigma;
  SNR = SNR * SNR;
  lambda = schranke(round(0.5 + ScTiRe / (SNR)), result.filter_ord, 0.5);
  HOHD_S(lambda);
}

void 
TDaten::HOHD_S(double lambda)
{
//cout<<"test2012"<<endl;  					//Tobias 2012	
//cout<<"a_fit.n_channels:"<<a_fit.n_channels<<endl;
//cout<<"a_fit.n_extra_channels:"<<a_fit.n_extra_channels<<endl;
//cout<<" a_fit.i_null:"<<a_fit.i_null<<endl;
//cout<<"a_fit.i_channel:"<<a_fit.i_channel<<endl;

g_type  g_up, g_down;
int  i,ii,Samples;
int  z;
int  lastmin, lastmax;
int  u_act;
int  n_open_bla = 0;

  InitDwellTimeHisto(a_fit.n_channels, a_fit.n_extra_channels, true);
  u_act = 0;
  Samples = Settings.SamplesToProcess + Settings.SamplesToSkip;
  i = Settings.SamplesToSkip;
  zero(g_up, result.filter_ord);
  zero(g_down, result.filter_ord); 
  lastmin = i;
  lastmax = i;
  HOHD_Spruenge.clear();
  if (Daten.SimSettings.actionpotentials == true) //Tobias: f�gt den Beginn der Zeitreihe (jeweiliger Abschnitt) in die Sprungliste ein
     {	
      ii=Daten.Settings.SamplesToSkip;
      // F�gt Sprung nach oben ein (nach unten geht nicht!), der vom Detektor mit sprung nach unten beantwortet wird
      // =>Pseudosprung am Anfang der Marker-Zeitreihe, der aber nicht! im 2D-Dwell-Time aufaucht 
      flush_Sprung(_HOHD, HOHD_Spruenge, ii, u_act, true, n_open_bla, n_open_bla, n_open_bla, n_open_bla);
      flush_Sprung(_HOHD, HOHD_Spruenge, ii, u_act, false, n_open_bla, n_open_bla, n_open_bla, n_open_bla);
     }

  while (i < Samples)
  {
    if (i > Anz_Samples) break;
    z = Wert[i];
    if (u_act == a_fit.n_channels)
    {
      lastmin = i;
    } else {
      g_up[0] = (z - a_fit.i_null) / a_fit.i_channel - 0.5 - u_act;
      cumsum(g_up, result.filter_ord);
      if (g_up[1] <= 0)
      {
        zero(g_up, result.filter_ord);
        lastmin = i;
      }
    }
    if (u_act == 0)
    {
      lastmax = i;
    } else {
      g_down[0] = (z - a_fit.i_null) / a_fit.i_channel + 0.5 - u_act;
      cumsum(g_down, result.filter_ord);
      if (g_down[1] >= 0)
      {
        zero(g_down, result.filter_ord);
        lastmax = i;
      }
    }
    if ((g_up[result.filter_ord] > lambda) && (u_act != a_fit.n_channels))
    {
      i = lastmin;
      lastmax = i;
      flush_Sprung(_HOHD, HOHD_Spruenge, i, u_act, true, n_open_bla, n_open_bla, n_open_bla, n_open_bla);
      zero(g_up, result.filter_ord);
      zero(g_down,  result.filter_ord);
    } else {
	  if ((g_down[result.filter_ord] <-lambda) && (u_act != 0))
	  {
        i = lastmax;
        lastmin = i;
        flush_Sprung(_HOHD, HOHD_Spruenge, i, u_act, false, n_open_bla, n_open_bla, n_open_bla, n_open_bla);
        zero(g_up, result.filter_ord);
        zero(g_down, result.filter_ord);
      }
	}
    i++;
  }
  old_index = Settings.SamplesToSkip;
}

double 
TDaten::schranke(int t_res, int ord, double half_jump_magnitude)
{
g_type  g;
int  i, t;
double lambda;

  for (i=0;i <= ord;i++) g[i] = 0.0;
  g[0] = half_jump_magnitude;           // half jump magnitude
  for (t=1;t <= t_res;t++) cumsum(g,ord);
  lambda = g[ord];
  g[0] = -half_jump_magnitude;          // nach dem event steigt g[ordd] noch weiter !!!
  for (t=1;t <= t_res;t++)              // solange wird g_up nicht nullgesetzt
  {
    cumsum(g,ord);
    if (g[ord] > lambda) lambda = g[ord];
  }
  neu_t_res = t_res * Settings.Samplininterval;
  return lambda;
}

inline void 
TDaten::cumsum(g_type &g, int ord)
{
int  i;

  for (i=1;i <= ord;i++) g[i]=g[i]+g[i-1];
}

inline void 
TDaten::zero(g_type &g, int ord)
{
int  i;

  for (i=0;i <= ord;i++)  g[i]=0;
}


//Indra&Tobias
void 
TDaten::SHD_AH()
{

int    Samples, i;
int    j, diagnosis_counter = 0, new_model;
double   epsilon,hmax, e;
bool  diagnosis_mode, positive_h, up;
unsigned short   z;
Tneighbor_type  neighbor;
int    n_open_old, n_extra_open_old;
int    n_open_new, n_extra_open_new;
//Indra
std::ofstream logouti("levelhist.log",std::ios::app);
int    level ,sublevel ,strom;
int    max_level, max_sublevel;
int   iold;
int   ii;
int    ***Level_Histogramm;

max_level = a_fit.n_channels;
max_sublevel = a_fit.n_extra_channels;

//erzeugen eines 3-dim Arrays f�r die Level_Histogramme, nach Anzahl der level/sublevel


Level_Histogramm = new int**[max_level+1]; 

for (level=0; level < (max_level+1); level++)
   Level_Histogramm[level] = new int*[max_sublevel+1];
 
  
for (level=0; level < (max_level+1); level++)
  for (sublevel=0; sublevel < (max_sublevel+1); sublevel++)
      { 
       Level_Histogramm[level][sublevel] = new int[4096];
      } 

for (level=0; level < (max_level+1); level++)
  for (sublevel=0; sublevel < (max_sublevel+1); sublevel++)
    for (strom=0; strom < 4096; strom++)	   
       Level_Histogramm[level][sublevel][strom] = 0;

   

  InitDwellTimeHisto(a_fit.n_channels,a_fit.n_extra_channels,true);
  SHD_Spruenge.clear();
  Samples = Settings.SamplesToProcess + Settings.SamplesToSkip;
  i = Settings.SamplesToSkip;
  
  //indra
  iold = Settings.SamplesToSkip;
  
  epsilon = 8.0 * a_fit.sigma * a_fit.sigma;
  create_neighbors(neighbor, i, 0, 0);
  diagnosis_mode = false;
  while (i < Samples)
  {
    if (i >= Anz_Samples) break;
    z = Wert[i];
    e = z - neighbor.m[0];
    for (j=1;j <= neighbor.n_neighbors;j++)
    {
      neighbor.h[j] = neighbor.h[j] + neighbor.p[j] * (e - neighbor.p[j]);
      if ((neighbor.h[j] <= 0) && (!diagnosis_mode))
      {
        neighbor.lastmin[j] = i;
        neighbor.h[j] = 0;
      }
    }
    if (!diagnosis_mode)
    {
      hmax = epsilon;
      for (j=1;j <= neighbor.n_neighbors;j++)
      {
        if (neighbor.h[j] > hmax)
        {
          hmax = neighbor.h[j];
          diagnosis_mode = true;
          diagnosis_counter = 0;
        }
      }
    }
    hmax = epsilon;
    new_model = 0;
    positive_h = false;
    if (diagnosis_mode)
    {
      diagnosis_counter++;
      for (j=1;j <= neighbor.n_neighbors;j++)
      {
        if (neighbor.h[j] > neighbor.greatest_h[j])
        {
          neighbor.greatest_h[j] = neighbor.h[j];
          neighbor.lastmax[j] = i;  // Zu welchem Zeitpunkt maximal
        }
        if (neighbor.h[j] > 0) positive_h = true;
      }
      if ((diagnosis_counter >= 7) || (positive_h = false)) // now a decision has to be made
      {      
        for (j=1;j <= neighbor.n_neighbors;j++)
        {
          if (neighbor.greatest_h[j] > hmax)
          {
            hmax = neighbor.greatest_h[j];
            new_model = j;
          }
        }
      }
    }
    if (new_model != 0)
    {
      if ((neighbor.lastmax[new_model]-neighbor.lastmin[new_model]) < 7)
      {
        i = neighbor.lastmax[new_model];
        diagnosis_mode = false;
        diagnosis_counter = 0;
        create_neighbors(neighbor, i, neighbor.n_open[0], neighbor.n_extra_open[0]);
      } else {
        i=neighbor.lastmin[new_model];  // geschaetzter Sprungzeitpunkt
        
        if (neighbor.n_open[new_model] > neighbor.n_open[0])  // Sprung nach oben ?
        {
          up = true;
        } else {
          if (neighbor.n_open[new_model] < neighbor.n_open[0])
          {
            up = false;
          } else { 
            if (neighbor.n_extra_open[new_model] > neighbor.n_extra_open[0])
            {
              up = true;
            } else {
              up = false;
            }
          }
        }             
        n_open_old = neighbor.n_open[0];
        n_extra_open_old = neighbor.n_extra_open[0];
        n_open_new = neighbor.n_open[new_model];
        n_extra_open_new = neighbor.n_extra_open[new_model];
        
        flush_Sprung(_SHD, SHD_Spruenge,i, new_model, up,n_open_old,n_extra_open_old,n_open_new,n_extra_open_new);

        if (iold < i)
          {     		
           for (ii = iold; ii < i; ii++)
             {
             Level_Histogramm[n_open_old][n_extra_open_old][Daten.Wert[ii]]++;
             iold = i;
             }
          } 
        diagnosis_mode = false;
        diagnosis_counter = 0;
        create_neighbors(neighbor, i, neighbor.n_open[new_model], neighbor.n_extra_open[new_model]);       
      }

    }
    i++;
  }
  old_index = Settings.SamplesToSkip;
  
  
  //indra  
  for (strom=0; strom < 4096; strom++)
    {
     logouti<<strom<<" ";	
     for (level=0; level < max_level+1; level++)
        for (sublevel=0; sublevel < max_sublevel+1; sublevel++)
       	   logouti<<Level_Histogramm[level][sublevel][strom]<<" ";
     logouti<<std::endl;	       	       	       	       	       	       	
    }
    
for (level=0; level < max_level+1; level++)
  for (sublevel=0; sublevel < max_sublevel+1; sublevel++)
    delete[] Level_Histogramm[level][sublevel];

for (level=0; level < (max_level+1); level++)
   delete[] *Level_Histogramm[level];   
   
delete[] **Level_Histogramm;  
     
}	


void 
TDaten::SHD()
{
int    Samples, i;
int    j, diagnosis_counter = 0, new_model;
double   epsilon,hmax, e;
bool  diagnosis_mode, positive_h, up;
unsigned short   z;
Tneighbor_type  neighbor;
int    n_open_old, n_extra_open_old;
int   n_open_new, n_extra_open_new;



  InitDwellTimeHisto(a_fit.n_channels,a_fit.n_extra_channels,true);
  SHD_Spruenge.clear();
  Samples = Settings.SamplesToProcess + Settings.SamplesToSkip;
  i = Settings.SamplesToSkip;
  
  
  epsilon = 8.0 * a_fit.sigma * a_fit.sigma;  //original von olli:
//epsilon = 8.0 * sqrt(a_fit.sigma);  // in Pascal aber Wurzel aus!!
  create_neighbors(neighbor, i, 0, 0);
  diagnosis_mode = false;
  while (i < Samples)
  {
    if (i >= Anz_Samples) break;
    z = Wert[i];
    e = z - neighbor.m[0];
    for (j=1;j <= neighbor.n_neighbors;j++)
    {
      neighbor.h[j] = neighbor.h[j] + neighbor.p[j] * (e - neighbor.p[j]);
      if ((neighbor.h[j] <= 0) && (!diagnosis_mode))
      {
        neighbor.lastmin[j] = i;
        neighbor.h[j] = 0;
      }
    }
    if (!diagnosis_mode)
    {
      hmax = epsilon;
      for (j=1;j <= neighbor.n_neighbors;j++)
      {
        if (neighbor.h[j] > hmax)
        {
          hmax = neighbor.h[j];
          diagnosis_mode = true;
          diagnosis_counter = 0;
        }
      }
    }
    hmax = epsilon;
    new_model = 0;
    positive_h = false;
    if (diagnosis_mode)
    {
      diagnosis_counter++;
      for (j=1;j <= neighbor.n_neighbors;j++)
      {
        if (neighbor.h[j] > neighbor.greatest_h[j])
        {
          neighbor.greatest_h[j] = neighbor.h[j];
          neighbor.lastmax[j] = i;  // Zu welchem Zeitpunkt maximal
        }
        if (neighbor.h[j] > 0) positive_h = true;
      }
      if ((diagnosis_counter >= 7) || (positive_h = false)) // now a decision has to be made
      {      
        for (j=1;j <= neighbor.n_neighbors;j++)
        {
          if (neighbor.greatest_h[j] > hmax)
          {
            hmax = neighbor.greatest_h[j];
            new_model = j;
          }
        }
      }
    }
    if (new_model != 0)
    {
      if ((neighbor.lastmax[new_model]-neighbor.lastmin[new_model]) < 7)
      {
        i = neighbor.lastmax[new_model];
        diagnosis_mode = false;
        diagnosis_counter = 0;
        create_neighbors(neighbor, i, neighbor.n_open[0], neighbor.n_extra_open[0]);
      } else {
        i=neighbor.lastmin[new_model];  // geschaetzter Sprungzeitpunkt
        
        if (neighbor.n_open[new_model] > neighbor.n_open[0])  // Sprung nach oben ?
        {
          up = true;
        } else {
          if (neighbor.n_open[new_model] < neighbor.n_open[0])
          {
            up = false;
          } else { 
            if (neighbor.n_extra_open[new_model] > neighbor.n_extra_open[0])
            {
              up = true;
            } else {
              up = false;
            }
          }
        }             
        n_open_old = neighbor.n_open[0];
        n_extra_open_old = neighbor.n_extra_open[0];
        n_open_new = neighbor.n_open[new_model];
        n_extra_open_new = neighbor.n_extra_open[new_model];
        
        flush_Sprung(_SHD, SHD_Spruenge,i, new_model, up,n_open_old,n_extra_open_old,n_open_new,n_extra_open_new);
        diagnosis_mode = false;
        diagnosis_counter = 0;
        create_neighbors(neighbor, i, neighbor.n_open[new_model], neighbor.n_extra_open[new_model]);
        
      }

    }
    i++;
  }
  old_index = Settings.SamplesToSkip;
      	       	       	       	       	       	       	
}	  

 
void 
TDaten::create_neighbors(Tneighbor_type &n,int lastjump, int open, int extra_open)
{
int   candi, i, k, g;
bool allowed = false;

  n.m[0] = a_fit.i_null + a_fit.i_channel * open + a_fit.i_extra_channel * extra_open;
  n.n_open[0] = open;
  n.n_extra_open[0] = extra_open;
  n.lastmin[0] = lastjump;

  for (i=1;i <= map;i++)
  {
    n.m[i] = n.m[0];
    n.p[i] = 0.0;
    n.h[i] = 0.0;
    n.greatest_h[i] = 0.0;
    n.lastmin[i] = lastjump;
  }
  candi = 0;
  for (k=0;k <= a_fit.n_extra_channels;k++)
  {
    for (g=0;g <= a_fit.n_channels;g++)
    {
      if ((k==extra_open) && (g==open))
      {
        allowed = false;  // aktuelles Modell[0]
      } else {
        switch (candidate_restriction)
        {
        	
          case 0:
            if ((abs((g+k) - (open+extra_open)) == 1) &&
                (abs(g - open ) <= 1 ) &&
                (abs( k - extra_open ) <= 1 ))
            {
               allowed = true;
            } else {
               allowed = false;
            }
            break;
          case 1: 
            allowed = true;
            break;
          case 2:
            if ((k+g <= a_fit.n_channels) &&
                (abs((g+k) - (open+extra_open)) <= 1) &&
                (abs(g - open) <= 1) &&
                (abs(k - extra_open) <= 1))
            {    
              allowed = true;
            } else {
              allowed = false;
            }
          case 3:
            if ((abs((g+k) - (open+extra_open)) <= 1) &&
                (abs(g - open) <= 1) &&
                (abs(k - extra_open) <= 1))
            {
              allowed = true;
            } else { 
              allowed = false;
            }
            break;
          case 4 :
            if ((abs((g+k) - (open+extra_open)) <= 2) &&
                (abs(g - open) <= 2) &&
                (abs(k - extra_open) <=2))
            {
              allowed = true;
            } else {
              allowed = false;
            }
            break;
        }
      }
      if (allowed)
      {
        candi++;
        n.m[candi] = a_fit.i_null + a_fit.i_channel * g + a_fit.i_extra_channel * k;
        n.p[candi] = (n.m[candi]-n.m[0]) / 2 ;
        n.n_open[candi] = g;
        n.n_extra_open[candi] = k;
      }
    }
  }
  n.n_neighbors = candi;
  if (a_fit.n_extra_channels == 0)
  {
    SNR = a_fit.i_channel/a_fit.sigma;
  } else {
    SNR = a_fit.i_extra_channel/a_fit.sigma;
  }
  SNR = SNR*SNR;
  for (i=1;i <= n.n_neighbors;i++) // lambda f�r jedes niveau berechnen, p[]: half jump magnitude
  {
    n.lambda[i] = schranke(round(0.5 + ScTiRe/(SNR)),//result.filter_ord
                                                     1,fabs(n.p[i]));
    // n.lambda[i]:=8*sqr(a_fit.sigma); PASCAL kommentar
  }
}


bool 
TDaten::flush_Sprung(detector_type d,
                     std::vector<TSprung> &Spruenge,
                     int   &index,
                     int   &new_model,
                     bool jump_up,
                     int   n_open_old,
                     int   n_extra_open_old,
                     int   n_open_new,
                     int   n_extra_open_new)
{
TSprung Sprung;
int  out = 0, new_out = 0;
int  event_length;
unsigned int  Spruenge_size = Spruenge.size();
 
 
  if (Spruenge_size < Spruenge.max_size())
  {
    Sprung.setnvs(a_fit.i_null, a_fit.i_channel, a_fit.i_extra_channel);  

    if ((d == _SHD) || (d == _HOSHD)) // Berechnung der Niveaus f�r SHD
    {
      // Strom ausrechnen: Level->Niveau, von dem der Sprung ausging
      Sprung.level_1(n_open_old, n_extra_open_old);
      Sprung.level_2(n_open_new, n_extra_open_new);
      // Position vermerken
      Sprung.position(index);
      out = n_extra_open_old +                       // erst alle kleinen Kan�le,
            (a_fit.n_extra_channels+1) * n_open_old; //  dann alle gro�en
      new_out = n_extra_open_new +                   // ->Reihenfolge nach Stromh�he geordnet
                (a_fit.n_extra_channels+1) * n_open_new;
    }
    if (d == _HOHD)  // Berechnung der Niveaus f�r HOHD
    {
      if (Spruenge.empty()) old_index = Settings.SamplesToSkip;  // Startwert wird auf Ausschnittsstartpunkt gesetzt
      out = new_model;  // dem Sub_level-Detector anpassen
      Sprung.level_1(new_model, 0);
      if (jump_up)
	  {
        new_model++;
	  } else {
        new_model--;
	  }
      new_out = new_model;
      Sprung.level_2(new_model, 0);
      Sprung.position(index);
    }
    Spruenge.push_back(Sprung); // Sprung in Liste einf�gen->f�r Darstellung �ber der Zeitreihe
    // Flush_Sprung := false  PASCAL return value;
    
    event_length = index - old_index; // L�nge des Impulses berechnen

    // PASCAL kommentar
    //if event_length >= result.n_bins then
    //  begin
    //    if not(event_too_long) then
    //      begin
    //        if (not(autom_steuer)) then
    //          begin
    //            if MessageBox(Parent^.HWindow,'event too long','abort Detection?', mb_OKCancel) = id_OK then
    //              index:=Anz_Samples;
    //          end;
    //      end;
    //    event_too_long:=true;
    //  end;

    if (event_length < 0)
    {
      event_length = 0;
      //if MessageBox(Parent^.HWindow,'event length < 0','abort Detection?', mb_OKCancel) = id_OK then index:=Anz_Samples;
      warning((char *)"event length < 0\n");
      // index = Anz_Samples; // Abbruch OK!
    }
    result.transitions[out][new_out]++; // Eintrag in �bergangsmatrix vornehmen
    result.detected_level_histo[out] += event_length; // Histogram �ber die Zeit in einem Niveau um Event-L�nge erh�hen
    if (event_length < result.n_bins)  // Sicherheitsabfrage: Events > n_bins (2000) werden vernachl��igt
    {
      result.dwell_time_histo_array[2 * out + (int)(!(jump_up))][event_length]++;  // erst: up dann: down f�r out=0
    }                                                                              // dann: up, down f�r out=1 usw
    old_index = index;
 
    // Neuzeichnen der Spr�nge w�hrend der Erkennung
    // PASCAL kommentar
    // DC:=GetDC(parent^.HWindow);
    // with ColorAndOrder do
    //   PlotSpruenge(Spruenge,DC,RGB(HOHD_R,HOHD_G,HOHD_B), True);
    // ReleaseDC(parent^.HWindow,DC);
    // PASCAL kommentar
 if(!TwoD_fitinaction)
    if (Spruenge_size % 512 == 0)
    {
      char buf[MaxTextLen];
      sprintf(buf, "# events: %i", Spruenge_size);
      Application.statustext(buf);
    }
	return false;
  } else {
    index = Anz_Samples;
    // MessageBox(Parent^.HWindow,'Collection Overflow','Detection aborted',mb_IconStop);
    warning((char *)"Collection Overflow - Detection aborted");
    return true;
  }
} 

/*void 
TDaten::plotSpruenge(std::vector<TSprung> &Spruenge, 
                     GdkColor   Farbe,
                     bool all)
{
std::vector<TSprung>::iterator iter;

  Application.drawcol(Farbe);
  if ((Settings.enlarge_current) && (Settings.Zoom))
  {
    //Spruenge.plotSpruenge(Settings, strom_min_old, strom_max_old, g_maxx, g_maxy,all);
    for(iter=Spruenge.begin();iter != Spruenge.end();iter++) (*iter).plot(strom_min_old, strom_max_old); 
  } else {
    for(iter=Spruenge.begin();iter != Spruenge.end();iter++) (*iter).plot(strom_min, strom_max);
  }
  TSprung::reset();
}*/

void
TDaten::separatechannels()
{
// vectoren von zeigern auf open und closed kanaele
std::vector< std::vector<TSprung>* > open_channels;  
std::vector< std::vector<TSprung>* > close_channels; 

// iterator ueber spruenge in sprunglisten
std::vector<TSprung>::iterator  sprung_iter;   

// zeiger auf sprunglisten 
std::vector<TSprung>*           sprung_vec_ptr;

  if (Daten.HOHD_Spruenge.empty()) return;
 
  // bei x kanaelen, geht G von 0 bis x
  // dito bei unterkanaelen  
  unsigned short channels = 0;
  unsigned short subchannels = 0;
  for(sprung_iter=Daten.HOHD_Spruenge.begin();sprung_iter != Daten.HOHD_Spruenge.end();sprung_iter++) 
  {
    if ((*sprung_iter).C_2() > channels) channels = (*sprung_iter).C_2();
    if ((*sprung_iter).U_2() > subchannels) subchannels = (*sprung_iter).U_2();
  }
  //vektor mit x=channel vektoren 
  SingleChannels.resize(channels);  

  // HOHD legt immer nur einen sprung ueber ein niveau ab 
  // und faengt immer bei 0 an!
  // wenn ich im 2 niveau starte hab' ich spruenge von 0 -> 1 und von 1 -> 2
  for (unsigned short i=0;i < channels;i++) close_channels.push_back(&SingleChannels[i]);
  open_channels.clear();

  // hier nun alle sprunge in verschiedene listen stecken
  for(sprung_iter=Daten.HOHD_Spruenge.begin();sprung_iter != Daten.HOHD_Spruenge.end();sprung_iter++) 
  {
    TSprung kopie(*sprung_iter);

	unsigned short c1 = kopie.C_1();
	unsigned short c2 = kopie.C_2();
    if (c1 < c2)
	{
   	  // jump up
      kopie.level_1(0, 0);
 	  kopie.level_2(1, 0); // werte anpassen

	  // liste zufaellig permutieren (zufaellige position schneller?)
      random_shuffle(close_channels.begin(), close_channels.end());

      sprung_vec_ptr = close_channels.back();   // diesen sprung im letzten
      sprung_vec_ptr->push_back(kopie);         // einzelkanal speichern

      close_channels.pop_back();                // der ist jetzt nicht mehr closed
      open_channels.push_back(sprung_vec_ptr);  // sonder open!
	} else {
      // jump down
      kopie.level_1(1, 0);
 	  kopie.level_2(0, 0); // werte anpassen

	  // liste zufaellig permutieren (zufaellige position schneller?)
      random_shuffle(open_channels.begin(), open_channels.end());

      sprung_vec_ptr = open_channels.back();    // diesen sprung im lezten
      sprung_vec_ptr->push_back(kopie);         // einzelkanal speichern

      open_channels.pop_back();                 // der ist jetzt nicht mehr open
      close_channels.push_back(sprung_vec_ptr); // sonder closed!
	}
  }
  return;
}

void TDaten::saveSChannels()
{
/* 2019
TDwell_2d dwell_2d;
TDwell_1d dwell_1d;  // fuer einzelne singechannels
char      buffer[MaxTextLen];

std::ofstream out;

  for(unsigned short i=0;i < Daten.SingleChannels.size();i++)
  {
	sprintf(buffer, "SChannel %i", i);
    dwell_1d.init(buffer);
    dwell_2d.init(buffer);

	// nicht addjumps, um SingleChannels nicht zu veraendern!
    dwell_2d.calc(&Daten.SingleChannels[i]);
	dwell_1d.calc(&Daten.SingleChannels[i]); 

	strcat(buffer,".dlb");
	out.open(buffer);

	if (!out)
	{
      char infotext[MaxTextLen];
      sprintf(infotext,"Can't open file: %s",buffer);

      warning(infotext);           // fehlermeldung wohin auch immer
      return;                      // keine daten lesen!
	}
    out << dwell_1d << dwell_2d;
	out.close();
  }
  return;
 */ 
}

void
TDaten::SaveDwell(TDwell &Item, std::ofstream &out)
{

  out << std::setw(10) << std::setprecision(0) << Item.x << " " 
	  << std::setw(10) << std::setprecision(3) << Item.y << std::endl;
}

//void 
//TDaten::SaveAddedDwell(TDwell &Item, ofstream &out)
//{
//  writeln(TextFile2, Item^.x:1:0, '   ',Item^.y:1:3);
//}

void 
TDaten::saveDwellTimeHisto(char *name)
{
// SaveDwellTimeHisto - Speichern der Dwell-Time-Histogramme in gemittelter Form 
// Wie man sie ungemittelt speichert, s.u. in den TFilterDaten-Routinen          
// Die 'cxx'-Files enthalten open+closed-Histogramme, die 'ixx'-Files open und   
// closed in verschiedenen Dateien. Gerade->open, Ungerade->Closed               
/*2019
char      FileName[fsPathName]; 	
short    i,j,k,l,m,n,o;

  //if MessageBox(Parent^.HWindow,'Quetch data?','Save dwell-time-histogram', mb_YesNo)= IDYES then
  // Die Frage, ob Gequetched oder nicht hat sich er�brigt, da immer alles gespeichert wird
  l = 0;
  o = 0;
  for (i=0;i <= a_fit.n_channels;i++)
  {
    for (j=0;j <= a_fit.n_extra_channels;j++)
	{
      sprintf(FileName, "%sc%i", name, l);

	  std::ofstream Textfile2(FileName,std::ios::in|std::ios::app);
      if (!Textfile2)
	  {
        char infotext[MaxTextLen];
        sprintf(infotext,"Can't open file: %s",FileName);

        warning(infotext);           // fehlermeldung 
        return;                      // keine daten lesen!
	  }
	  Textfile2 << std::endl << result.n_bins << std::endl
		        << "jump up + jump down, n_Channel: " << i 
				<< " n_Extra_Channel: " << j << std::endl;
       	
	  if (!Daten.result.dwell_time_histo_array.empty())
	  {
        SDwellTime.LevelNo = i;
        SDwellTime.SublevelNo = j;
        SDwellTime.openclosed = true;
        mittel();

        std::vector<TDwell>::iterator iter;
        for (iter=QuenchedDwellTime.begin(); iter != QuenchedDwellTime.end(); iter++)
		{
          SaveDwell(*iter, Textfile2);
		}
	  }
      Textfile2.close();
      if (((i > 0) || (j > 0) || (k = 0)) &&
          !((i = a_fit.n_channels) && (j = a_fit.n_extra_channels) && (k=0)))
      {
        // for k:=1 to 0 do ??????
        for (k=1;k<=0;k++)
		//for(k=1;k >= 0;k--)
		{  
          if (l >= 100) break;
          sprintf(FileName,"%si%i", name, o);
		  std::ofstream Textfile(FileName);
          if (!Textfile)
		  {
            char infotext[MaxTextLen];
            sprintf(infotext,"Can't open file: %s",FileName);

            warning(infotext);           // fehlermeldung 
            return;                      // keine daten lesen!
		  }

          Textfile << std::endl << result.n_bins << std::endl; 

          if (k == 0)
		  {
            Textfile << "jump up, n_channel: " << i << " n_extra_channel: " << j << std::endl;     
          } else {
            Textfile << "jump down, n_channel: " << i << " n_extra_channel: " << j << std::endl;     
          }
		  if (!Daten.result.dwell_time_histo_array.empty())
          {
            SDwellTime.LevelNo = i;
            SDwellTime.SublevelNo = j;
            if (k == 0) 
			{
			  SDwellTime.open = true;
			} else {
			  SDwellTime.open = false;
			}
            SDwellTime.openclosed = false;
            mittel();

            std::vector<TDwell>::iterator iter;
            for (iter=QuenchedDwellTime.begin(); iter != QuenchedDwellTime.end(); iter++)
			{
              SaveDwell(*iter, Textfile);
			}
		  }
          Textfile.close();
          o++;
		}
        l++;
	  }
	}
  }
  // else:
  // Immer alles speichern -> else Zweig immer ausf�hren, immer auch nicht gequetcht speichern
  // hier aus TFilterDaten
  l=0;
  o=0;
  for (i=0;i <= a_fit.n_channels;i++)
  {
    for (j=0;j <= a_fit.n_extra_channels;j++)
    {
      n=j+(a_fit.n_extra_channels+1)*i;
      n=2*n;
      sprintf(FileName,"%sh%i",name,l/2);
      std::ofstream Textfile2(FileName);
      Textfile2 << std::endl;
      Textfile2 << result.n_bins << std::endl;
      Textfile2 << "jump up + jump down, n_Channel: " << i << " n_Extra_Channel: " << j << std::endl;
      if (!Daten.result.dwell_time_histo_array.empty())
      {
        for (m=0;m <= result.n_bins-1;m++) 
        {
          Textfile2 << m << "   " 
                    << result.dwell_time_histo_array[n][m] + 
                       result.dwell_time_histo_array[n+1][m] 
                    << std::endl;
        }
      }
      Textfile2.close();
      for (k=1;k >= 0;k--)
      {

		if (l >= 100) break;
        if (
            ((i > 0) || (j > 0) || (k == 0)) &&
            !((i == a_fit.n_channels) && (j == a_fit.n_extra_channels) && (k==0))
           )
        { 
          sprintf(FileName,"%sd%i",name,o);
          std::ofstream Textfile(FileName);
          Textfile << std::endl;
          Textfile << result.n_bins << std::endl;
          if (k == 0)
          {
            Textfile << "jump up, n_Channel: " << i << " n_Extra_Channel: " << j << std::endl;
          } else {
            Textfile << "jump down, n_Channel: " << i << " n_Extra_Channel: " << j << std::endl;
          }
          if (!Daten.result.dwell_time_histo_array.empty())
          {
            for (m=0;m <= result.n_bins-1;m++) 
            {
              Textfile << m << "   " << result.dwell_time_histo_array[n+k][m] << std::endl;
            }
          }
          Textfile.close();
          o++;
        }
        l++;
      }
    }
  }
  2019*/
}

void 
TDaten::saveDetectedLevelHisto(char *name)
{
int  sum;
int  i, j, k,l;

double   popen0 = 0, popenl = 0;
std::ofstream Textfile (name);

  Textfile << std::endl;
  sum=0;
  for (i=0;i <= (a_fit.n_channels + 1) * (a_fit.n_extra_channels + 1) - 1;i++)
  {
    sum += result.detected_level_histo[i];
  }
  Textfile << "Level  Current/pA   # Level  P(Level)  Sigma of Level" << std::endl;
  Textfile << "-----------------------------------------------------" << std::endl;
  k=0;
  for (i=0;i <= a_fit.n_channels;i++)
  {
    for (j=0;j <= a_fit.n_extra_channels;j++)
    {
      l = j + (a_fit.n_extra_channels + 1) * i;

	  Textfile.setf(std::ios::fixed,std::ios::floatfield);
      Textfile << std::setw(4) << k
               << std::setw(12) << std::setprecision(2) << (i * a_fit.i_channel + j * a_fit.i_extra_channel) * Settings.VFactor/Settings.Gain;
      if (l == 0)
	  {
        Textfile << std::setw(10) << Settings.SamplesToProcess - sum + result.detected_level_histo[l];
	  } else {
        Textfile << std::setw(10) << result.detected_level_histo[l];
	  }
      if (sum != 0)
	  {
        if (l == 0)
		{
          Textfile << std::setw(10) << ((double)(Daten.Settings.SamplesToProcess - sum + Daten.result.detected_level_histo[l]) / 
			                       (double)Daten.Settings.SamplesToProcess);
		} else {
          Textfile << std::setw(10) << ((double)Daten.result.detected_level_histo[l] / (double)Daten.Settings.SamplesToProcess);
		}
	  }

	  // offenwahrscheinlichkeiten speichern
      if (l == 0)
	  {
		popen0 = 1.0 - exp((1.0 / (double)Daten.a_fit.n_channels) *
                 log((Daten.Settings.SamplesToProcess - sum + 
				 (double)Daten.result.detected_level_histo[l]) /
				 (double)Daten.Settings.SamplesToProcess));
	  }
	  if (l == (Daten.a_fit.n_extra_channels+Daten.a_fit.n_channels))
	  {
		popenl = exp((1.0 / Daten.a_fit.n_channels) * 
			     log(Daten.result.detected_level_histo[l] /
                 (double)Daten.Settings.SamplesToProcess));
	  }

	  Textfile << std::setw(16) << std::setprecision(3) << calcStD(i,j)
               << std::endl;
	  Textfile.setf(std::ios::floatfield); //Tobias 2017 weg: 0,
      k++;
    }
  }
  Textfile << "-----------------------------------------------------" << std::endl;

  Textfile << std::endl << "Offenwahrscheinlichkeiten:" << std::endl;
  Textfile << "berechnet aus Level 0: ";
  Textfile << std::setprecision(3) << popen0 << std::endl;
  Textfile << "berechnet aus Level " << Daten.a_fit.n_extra_channels+Daten.a_fit.n_channels << ": "; 
  Textfile << std::setprecision(3) << popenl << std::endl;
}

void 
TDaten::saveTransitionMatrix(char *name)
{
int  i,j;

  std::ofstream TextFile (name);
  TextFile << std::endl;
  TextFile << "Transitions:" << std::endl;
  TextFile << std::endl;
  for (i=0;i <= (a_fit.n_channels+1)*(a_fit.n_extra_channels+1)-1;i++)
  {
    for (j=0;j <= (a_fit.n_channels+1)*(a_fit.n_extra_channels+1)-1;j++)
    {
      if (i==j)
      {
        TextFile << std::setw(5) << 0 << "  ";
      } else {
        TextFile << std::setw(5) << result.transitions[i][j] << "  ";
      }
    }
    TextFile << std::endl;
  }
}

void  
TDaten::clearDetectedLevels()
{
// waere  wohl besser bei Tresult aufgehoben
int i;
  for (i=0; i <= map;i++) result.detected_level_histo[i] = 0;
}

void
TDaten::clearTransitions()
{
// auch besser bei Tresult
int i,j;

  for (i=0;i <= map;i++)
    for (j=0;j <= map;j++)
      result.transitions[i][j] = 0;
}

double 
TDaten::findLambda()
{
double  la,sa;
double  lb;//,sb;
double  lc,sc;
  sa=0;
  la=maximum;
  do
  {
    HOHD_S(la);
    sa=HOHD_Spruenge.size();
    la=la/2;
  } while (sa <= 100);

  la=2*la;
  lb=2*la;
  do
  {
    lc=(la+lb)/2;
    HOHD_S(lc);
    sc=HOHD_Spruenge.size();
    if (sc > 100)
    {
      la=lc;
    } else {
      lb=lc;
    }
  } while (fabs(sc-100) > 2);

  return lc;
}

double 
TDaten::zeitauflsg(double schr, int ord)
{
int   t1, t2, t3, t4;  // 4 Zeitaufl�sungen
double  c1, c2, c3;      // 3 Schranken

  t1=1;     // kleinste Aufl�sung: 1 Sample
  t2=30000; // gr��te Auflk�sung: 1000 Samples
  t3=0;
  c1=schranke(t1,ord,0.5);
  //writeln(c1);
  c2=schranke(t2,ord,0.5);
  //writeln(c2);
  //Zeitauflsg:=0; (PASCAL return value)
  if ((schr > c2) || (schr < c1))
  {
    // MessageBox(parent^.HWindow, 'can`t find tres', 'tres > 30000 or < 1', mb_IconStop);
    warning ((char *)"can`t find tres\ntres > 30000 or < 1");
  } else {
    do
    {
      t4=t3;
      t3=(t1+t2) / 2;
      c3=schranke(t3,ord,0.5);
      if (schr > c3)
      {
        t1=t3;
        c1=c3;
      } else {
        t2=t3;
        c2=c3;
      }
    } while (!((fabs(schr - c3) < schr/20) || (t3 == t4)));
  }
  return t3;
}

void
TDaten::History()
{
  switch (HText.detector)
  {
    case _HOHD: 
      HText.detected_jumps = HOHD_Spruenge.size();
      break;
    case _SHD:  
      HText.detected_jumps = SHD_Spruenge.size();
      break;
    case _HOSHD: 
      HText.detected_jumps = HOSHD_Spruenge.size();
      break;
    case _SHD_AH:
      break;					//Tobias 2013 
  }
  HText.sigma_of_channel = 0;
}

double 
TDaten::calcStD(int n_channel, int n_extra_channel)
{
double  level;

  HText.level = n_channel;
  HText.sublevel = n_extra_channel;
  //Spruenge:=nil;
  level = a_fit.i_null
        + a_fit.i_channel * HText.level
        + a_fit.i_extra_channel * HText.sublevel;

  switch (HText.detector)
  { 
    case _HOHD:
      level = a_fit.i_null + a_fit.i_channel * HText.level;
      if (!HOHD_Spruenge.empty()) HText.sigma_of_channel = TSprung::calcSD(level,HOHD_Spruenge);
      return HText.sigma_of_channel;
      break;
    case _SHD:
      if (!SHD_Spruenge.empty()) HText.sigma_of_channel = TSprung::calcSD(level,SHD_Spruenge);
      return HText.sigma_of_channel;
      break;
    case _HOSHD:
      if (!HOSHD_Spruenge.empty()) HText.sigma_of_channel = TSprung::calcSD(level,HOSHD_Spruenge);
      return HText.sigma_of_channel;
      break;
    case _SHD_AH:
      break;					//Tobias 2013   
  }
  // spruenge empty!
  HText.sigma_of_channel = 0.0;
  return 0;
}

inline void
TDaten::FindMax(const TDwell &Item, double &max)
{
  if (Item.y > max) max = (Item.y)+1;
}

double 
TDaten::DetermineMax(int index)
{
// Maximum des DwellTimeHistograms bestimmen -> f�r Normierung auf Bildschirm

double max = 1;

  if (!Daten.result.dwell_time_histo_array.empty())
  {
    std::vector<TDwell>::const_iterator iter;
    for (iter=QuenchedDwellTime.begin(); iter != QuenchedDwellTime.end(); iter++)
	{
	  FindMax(*iter, max);
	}
  } else {
    max = 0;
  }
  return max;
}

TTargetWeight* 
TDaten::GetWeight()
{
TTargetWeight*  weight;

  weight = &Weightdown;
  if (SDwellTime.open)
  {
    weight = &Weightup;
  }
  if (SDwellTime.openclosed)
  {
    weight = &TargetWeight;
  }
  return weight;
}

void
TDaten::StartIndex(TDwell &Dwell, int &index, double &x, double &y)
{
int i=0;

  if (Dwell.y > y)
  {
    y = Dwell.y;
    x = Dwell.x;

	while (&QuenchedDwellTime[i] != &Dwell && i < (int)QuenchedDwellTime.size()) i++; 
    index = i; //QuenchedDwellTime.IndexOf(Dwell);
  }
}

int 
TDaten::CalcStartIndex()
{
int  index;
double x,y;
std::vector<TDwell>::iterator iter;

  y = 0;
  index = -1;

  for(iter=QuenchedDwellTime.begin();iter != QuenchedDwellTime.end();iter++) 
  {
	 StartIndex(*iter,index,x,y);
  }
  return index;
}
  
void 
TDaten::NoSamples()
{
// NoSamples: Aus Grow-Factor wird die Zahl der Samples nach dem Mittelungs-
// prozess errechnet. Gegenst�ck zu ExpEinteil

double  w,s,g;
int   n;

  if (!Daten.result.dwell_time_histo_array.empty())
  {
    w = SDwellTime.GrowFactor;
    s = 1;
    n = 1;
    g = result.n_bins;
    do
    {
      s = s + w;
      w = w * SDwellTime.GrowFactor;
      n = n + 1;
    } while (s <= g - 20);
    SDwellTime.Samples = n;
  }
}

void 
TDaten::ExpEinteil()
{
// ExpEinteil: Expotentielle Mittelung im DwellTimeHistogram vornehmen
// Der Wachstumsfaktor w wird mit der Bisektionsmethode bestimmt man gibt die
// Zahl der Samples vor 

double n,w1,w2,w3,test1,test2,test3;

  if (!Daten.result.dwell_time_histo_array.empty())
  {
    n = SDwellTime.Samples;
    w1 = 1.01;
    w2 = 5;
    test1 = (1-exp(log(w1)*n))/(1-w1)-result.n_bins;  // Geometrische Reihe
    test2 = (1-exp(log(w2)*n))/(1-w2)-result.n_bins;
    if (((test1 > 0) && (test2 > 0)) || ((test1 < 0) && (test2 < 0)))
    {
      std::cout<<"Please choose other start-values\nCan`t find root"<<std::endl;
    }
    if ((test1 > 0) && (test2 < 0))
    {
      test3 = test1;
      test1 = test2;
      test2 = test3;
      w3 = w1;
      w1 = w2;
      w2 = w3;
    }
    if ((test1 < 0) && (test2 > 0))
    {
      do
      {
        w3 = (w1 + w2) / 2;
        test3 = (1-exp(log(w3)*n))/(1-w3)-result.n_bins;
        if (test3 > 0)
        {
          w2 = w3;
          test2 = test3; 
        } else {
          w1 = w3;
          test1 = test3;
        }
      } while (!((fabs(test1) < 0.000001) && (fabs(test2) < 0.000001)));
      SDwellTime.GrowFactor = w3;
    }
  }
}
  
void 
TDaten::mittel()
{
// Exp. Mittelung f�r ein bestimmtes Dwell-Time-Histogram
// modifiziert, so da� auch open/closed Histogramme aufaddiert werden k�nnen

int  i=0,j=0,k=0,l=0;
double a=0,b=0,w=0,s=0,t=0,m=0,mx=0,g=0;

  if (SDwellTime.ExpAvarage)
  {
    if (!Daten.result.dwell_time_histo_array.empty())
    {
      QuenchedDwellTime.clear();
      if (SDwellTime.openclosed)
      {
        k = 2 * SDwellTime.SublevelNo + 2 * SDwellTime.LevelNo * (a_fit.n_extra_channels + 1);
      } else {
        k = (int)(!(SDwellTime.open)) + 2 * SDwellTime.SublevelNo + 2 * SDwellTime.LevelNo * (a_fit.n_extra_channels + 1);
	  }
      if (k <= (2 * (a_fit.n_channels + 1) * (a_fit.n_extra_channels + 1)))
      {
        w = 1;
        s = 1;
        j = 0;
        i = 0;
        m = 0;
        mx = 0;
        t = 0;
        g = result.n_bins;
        do
		{
          i = i+1;
          j = j+1;
          if (i <= result.n_bins)
          {
            if (SDwellTime.openclosed)
			{
              a = result.dwell_time_histo_array[k][i-1];
              b = result.dwell_time_histo_array[k+1][i-1];
              m = m + a + b;
			} else {
              m = m + result.dwell_time_histo_array[k][i-1];
			}
            //if (a <> 0) or (b <> 0) then PASCAL kommentar
            {
              mx = mx+i;
			}
            if (j >= w)
			{
		      //Einfach mal nicht durch Zahl der Samples teilen		
              //QuenchedDwellTime^.Insert(New(PDwell, Init(mx/j,m))); // PASCAL kommentar

              // nur abspeichern, wenn Wert <> 0, da 0-Werte nicht gefittet werden.}
              if (m != 0)
			  {
				TDwell dwell(mx/j,m/j);
                QuenchedDwellTime.push_back(dwell);
			  }
              w = w * SDwellTime.GrowFactor;
              t = s;
              s = s + w;
              m = 0;
              mx = 0;
              j = 0;
			}
		  }
		   //until t > g - 20; //PASCAL kommentar
		} while (i < result.n_bins);
	  }
	}
  } else {
    if (!Daten.result.dwell_time_histo_array.empty())
    {
      QuenchedDwellTime.clear();
      if (SDwellTime.openclosed)
      {
        k = 2 * SDwellTime.SublevelNo + 2 * SDwellTime.LevelNo * (a_fit.n_extra_channels + 1);
	  } else {
        k = (int)(!(SDwellTime.open)) + 2 * SDwellTime.SublevelNo + 2 * SDwellTime.LevelNo * (a_fit.n_extra_channels + 1);
	  }
      if (k <= (2 * (a_fit.n_channels + 1) * (a_fit.n_extra_channels + 1)))
      {
        j = 0;
        m = 0;  // Mittelwert y
        l = 0;  // z�hlt Samples eines Intervals
        mx = 0; // Mittelwert x
        for (i=0;i <= result.n_bins-1;i++)
        {
          j = j + 1;
          l = l + 1;
          if (SDwellTime.openclosed)
          {
            a = result.dwell_time_histo_array[k][i];
            b = result.dwell_time_histo_array[k+1][i];
            m = m + a + b;
		  } else {
            m = m + result.dwell_time_histo_array[k][i];
          }
          mx = mx + j;
          if  ((i % (2000 / SDwellTime.Samples)) == 0)
		  {
			TDwell dwell(mx/l,m/l);
            QuenchedDwellTime.push_back(dwell);
            m = 0;
            l = 0;
            mx = 0;
		  }
		}
	  }
	}
  }
}
  
void 
TDaten::axes(double mx, double my, double xs, double ys)
{
/*2019	
// Achsenkreuz zeichnen

short   i;
double  x,y;
char     textbuf[MaxTextLen];

  plot_strom_height(mx,0,mx,ys,my,false,false, g_maxx, g_maxy, Application.border_gc);
  plot_strom_height(xs,0,mx,ys,my,true,false, g_maxx, g_maxy, Application.border_gc);
  plot_strom_height(xs,0,mx,my, my, true,false, g_maxx, g_maxy, Application.border_gc);
  x = xs;
  for (i=1;i <= 10;i++)
  {
    x = x + result.n_bins/(10);    // x-Achse in 1/10 n_bins einteilen
    plot_strom_height(x,0,mx,ys,my,false,false, g_maxx, g_maxy, Application.border_gc);
    plot_strom_height(x,0,mx,2.0/3*ys,my,true,false, g_maxx, g_maxy, Application.border_gc);
  }
  y = ys;
  for (i=1;i <= (short)my;i++)
  {
    y = y + 1;                      // y-Achse in ganze Ereignisse teilen
    plot_strom_height(xs,0,mx,y,my,false,false, g_maxx, g_maxy, Application.border_gc);
    plot_strom_height(2.0/3*xs,0,mx,y,my,true,false, g_maxx, g_maxy, Application.border_gc);
  }

  // Ausgabe der Achsenbeschriftung

  x = result.n_bins * Settings.Samplininterval/10;
  sprintf(textbuf,"x: %.0f Microseconds / Unit",x);

  // systemfont ermitteln
  GdkGCValues values;
  gdk_gc_get_values((Application.DLG_Main)->style->bg_gc[GTK_STATE_NORMAL],&values);

  // Ausgabe der Texte
  gdk_draw_string(Daten.MemPix,
			      values.font,
			      Application.dots_gc,
			      g_maxx / 2,
			      2 + gdk_string_height(values.font,"I"),
			      textbuf);
  gdk_draw_string(Daten.MemPix,
			      values.font,
			      Application.dots_gc,
			      g_maxx / 2,
			      2 * (2 + gdk_string_height(values.font,"I")),
			      "y: 1 Event / Unit");
2019*/			      
}

void 
TDaten::zoomin()
{
// Darstellung der Dwell-Time-Histogramme vergr��ern

  zoomfactor = zoomfactor + 1;
}
  
void 
TDaten::zoomout()
{
// Darstellung der Dwell-Time-Histogramme verkleinern

  if (zoomfactor > 1) zoomfactor = zoomfactor - 1;
}

void
TDaten::ExpFit()
{
// Der Expotentialfit

short        i,j,k,index,iter;
double       e,x,y;
TTargetWeight* weight;

  k = closing;
  weight = GetWeight();
  if (SDwellTime.open)
  {
    k = opening;
  }
  if (SDwellTime.openclosed)
  {
    k = openclosing; // k auf closing, opening oder open + closing setzen
  }

  // Es gib verschiedene Zeitkonstanten f�r diese Histogramme
  // Neues Fit-Objekt mit gew�hlten Start-Zeitkonstanten anlegen
  TExpFit  E_fit(DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].tau);

  //index = Zeitpunkt der gr��ten Amplitude im Histogram, Werte vorher werden nicht beachtet
  index = CalcStartIndex();
  //Zeitpunkt abspeichern (f�r Darstellung: Paint-Routine zeichnet e-Fkt. erst ab diesem Zeitpkt.)
  if (index >= 0)
  {
    DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].startx = QuenchedDwellTime[index].x;
    j = 0;
    iter = 0;
    e = 0;

    for (i=index;i <= (int)QuenchedDwellTime.size()-1;i++)
	{
      y = QuenchedDwellTime[i].y;
      if (y != 0)
	  {
        j++;
        x = QuenchedDwellTime[i].x;
        E_fit.ax[j] = x;              // x, y -Werte des Fit-Objekts mit Werten des Gequetschten Histogramms initialisieren
        E_fit.ay[j] = y;
        if ((*weight)[SDwellTime.LevelNo][j] != 0) // wir erlauben keine 0-Wichtung, wer 0-wichtet ist selber schuld!
        {
          E_fit.asig[j] = sqrt(fabs(y)) / (*weight)[SDwellTime.LevelNo][j];
          // Sigma-Werte des Fit-Objekts mit Gewichten initialisieren
		} else {
          E_fit.asig[j] = sqrt(fabs(y));
        }
	  }
	}
    E_fit.NumPoints = j; // Fitroutine Zahl der Datenpunkte �bergeben
    DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].Points = E_fit.NumPoints; // Dialog Zahl der Datenpunkte �bergeben


    e = E_fit.fit(0.0000001, iter); // Fitten

    // Ergebnisse kopieren, damit sie angezeigt werden k�nnen
    DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].Error = e;  // Dialog Fehler �bergeben
    DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].iter += iter; // im Dialog Zahl der Iterationen erh�hen
    for (i=1;i <= NMaxTimeConstants;i++) // Amplituden extra sichern, da Dispose nur Zeitkonstanten liefert
    {
	  DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].Ampli[i] = E_fit.aa[i];
	}
    DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].IsTargetFit = false;
  }
  // Zeitkonstanten sichern
  E_fit.done(DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].tau);
}

char*  
TDaten::GetxyString(int &index, char* s)
{
// GetxyString liefert einen Punkt, ZeitDauer und Zahl der Werte dieser Zeitdauer aus dem aktuellen,
// gequetschten Dwell-Histogram
// s: zeiger auf textbuffer in dem string liegt

int  i,startindex;
double  x,y;

  startindex = CalcStartIndex();
  i = startindex + index;  // Index beginnt bei Null
  if ((i < (int)QuenchedDwellTime.size()) && (i >=0))
  {
    x = (QuenchedDwellTime[i]).x;
    x = x * Settings.Samplininterval;
    y = (QuenchedDwellTime[i]).y;
    
	sprintf(s,"%3i %.1f %.3f",index, x, y);
    return s;
  }	else {
    return NULL;
  }
}
  
void   
TDaten::ResetWeight()
{
// ResetWeight setzt alle Gewichtsfaktoren auf 1

int  i;
TTargetWeight*  weight;

  weight = GetWeight();
  if (SDwellTime.LevelNo <= TargetFitMaxChannels)
  for (i=1;i <= quenched_bins_max;i++)
  {
    (*weight)[SDwellTime.LevelNo][i] = 1;
  }
}

  
void   
TDaten::generate_gnu_plot(char* name)
{
// gnuplot Graphik erzeugen und speichern

const int hoch = 39;

char*   p;
char    FileName[MaxTextLen];
char    DatenFileName[MaxTextLen];
char    Datename[MaxTextLen];
char    Datename2[MaxTextLen];
char    Title[MaxTextLen];
char    Zahl[MaxTextLen];
char    Zahl2[MaxTextLen];
char    Expofunk[MaxTextLen];
int  i,startindex,k;
double a,b,x,y;

  strcpy(FileName,name);
  strcpy(DatenFileName,name);
  p = strrchr(DatenFileName,'.');
  if (p != NULL) *p = 0; // endung abschneiden 
  strcat (DatenFileName, ".dat");

/*  sprintf (Title, "%s: Level %i sublevel %i jump down", FName, SDwellTime.LevelNo, SDwellTime.SublevelNo);					//2019 tOBIAS							
  if (SDwellTime.open) sprintf (Title, "%s: Level %i sublevel %i jump up", FName, SDwellTime.LevelNo, SDwellTime.SublevelNo);
  if (SDwellTime.openclosed) sprintf (Title, "%s: Level %i sublevel %i jump up + down", FName, SDwellTime.LevelNo, SDwellTime.SublevelNo);
*/
  // Plt-File f�r gnuplot erzeugen
  std::ofstream TextFile(FileName);
  TextFile << "set terminal windows color \"Arial\" 18" << std::endl
           << "set output" << std::endl
           << "set noclip points" << std::endl
           << "set clip one" << std::endl
           << "set noclip two" << std::endl
           << "set border" << std::endl
           << "set boxwidth" << std::endl
           << "set dummy x,y" << std::endl
           << "set format x \"%g\"" << std::endl
           << "set format y \"%g\"" << std::endl
           << "set nogrid" << std::endl
           << "set key" << std::endl
           << "set nolabel" << std::endl
           << "set logscale x 10" << std::endl
           << "set offsets 0, 0, 0, 0" << std::endl
           << "set nopolar" << std::endl
           << "set noparametric" << std::endl
           << "set view 60, 30, 1, 1" << std::endl
           << "set samples 100,100" << std::endl
           << "set isosamples 10,10" << std::endl
           << "set surface" << std::endl
           << "set nocontour" << std::endl
           << "set clabel" << std::endl
           << "set nohidden3d" << std::endl
           << "set size 1,1" << std::endl
           << "set data style points" << std::endl
           << "set function style lines" << std::endl
           << "set xzeroaxis" << std::endl
           << "set yzeroaxis" << std::endl
           << "set tics in" << std::endl
           << "set ticslevel 0.5" << std::endl
           << "set xtics" << std::endl
           << "set ytics" << std::endl
           << "set title \"" << Title << "\" 0,0" << std::endl
           << "set ylabel \"Number of events \" 0,0" << std::endl
           << "set xlabel \"Time / �s\" 0,0" << std::endl
           << "set autoscale xy" << std::endl
           << "set zero 1e-08" << std::endl;

  //Tobias 2019 sprintf(Datename,"%c%s%c title %cDwell- time histogram%c",hoch, DatenFileName, hoch, hoch, hoch);
  // gefittete exp. fkt. mit in gnuplot Graphik einbauen

  k = 0;
  if (SDwellTime.open) k = opening;
  if (SDwellTime.openclosed) k = openclosing;

  strcpy(Expofunk," ");

  // aus Zeitkonstanten und Faktoren die Exp.Funk. zusammenbasteln
  for (i=1;i <= NMaxTimeConstants;i++)
  {
    a = DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].tau[i] * Settings.Samplininterval; // in us belassen
    if (a == 0) break;
    b = DwellLevel[SDwellTime.LevelNo][SDwellTime.SublevelNo][k].Ampli[i];
    if (b == 0) break;

	sprintf(Zahl,"%f", b);
    sprintf(Zahl2, "%f", a);
    if (i == 1) strcat(Expofunk," , ");
    if (i > 1)  strcat(Expofunk," + ");
    strcat (Expofunk, Zahl);
	strcat (Expofunk, "*exp( -(x)/");
	strcat (Expofunk, Zahl2);
	strcat (Expofunk, ")");
  }
  if (strlen(Expofunk) > 1)
  {
	// titel nur anhaengen, wenn fit ex.
	sprintf(Datename2," title %cFit%c",hoch,hoch);
    strcat (Expofunk, Datename2);
  }

  TextFile << "plot " << Datename << Expofunk << std::endl
           << "pause -1" << std::endl;
  TextFile.close();
  // Datenfile mit Dwell-Time Histogram erzeugen
  TextFile.open(DatenFileName);

  startindex = CalcStartIndex(); // Startindex hat gr��tem Wert
  if (startindex >=0)
  {
    for (i=startindex;i <= (int)QuenchedDwellTime.size()-1;i++) // alle Datenpunkte ab gr��ten Wert durchgehen
	{		
      y = QuenchedDwellTime[i].y;
      if (y != 0)
	  {
        x = QuenchedDwellTime[i].x * Settings.Samplininterval; // in us belassen
        TextFile << x << "   " << y << std::endl; //4:3
	  }
	}
  }
}

void   
TDaten::exec_gnuplot(char*  name)
{
/*
// gnuplot ausf�hren und Graphik anzeigen lassen

short  returnvalue;
char cmd[] = "wgnuplot";
char arg[MaxTextLen];

  //brauche "" um filename
  strcpy(arg,"\"");
  strcat(arg,name);
  strcat(arg,"\"");

  generate_gnu_plot(name);

  errno=0;
  returnvalue = _spawnlp (_P_NOWAIT, cmd, cmd, arg, NULL); 

  if (returnvalue == -1) 
  {
    switch (errno)
	{
	  case ENOENT:
	    warning((char *)"file not found\nerror executing gnuplot");
        break;
	  case ENOMEM:
        warning((char *)"not enough memory\nerror executing gnuplot");
		break;
	  default:
        warning((char *)"error executing gnuplot");
		break;
	}
  }
*/
}



void 
TDaten::algAbspeichern(VecD lambda[], VecD ampli[], std::vector<TUntermatrix> &Dimension)
{
int  i,j;
int  start, stop, dim;
double l;

  for (i=0;i <= a_fit.n_channels;i++)
  {
    start = Dimension[i].start_index;
    stop = Dimension[i].stop_index;
    dim = Dimension[i].dimension;
    for (j=1;j <= dim;j++)
	{
      l = lambda[i][j];
      if (l == 0) l = 1e-10;
      DwellLevel[i][0][openclosing].tau[j] = -(1000 / Settings.Samplininterval) / l;  // auf Samples umrechnen // = -200/l;
      DwellLevel[i][0][openclosing].Ampli[j] = ampli[i][j];
	} 
    for (j=dim+1;j <= NMaxTimeConstants;j++) 
    {                                       
      // etwas problematisch: 
      // was ist wenn unter-Matrix gr��ere Dimension
	  // als NMaxTimeConstants?
      DwellLevel[i][0][openclosing].Ampli[j] = 0;
      DwellLevel[i][0][openclosing].tau[j] = 0;
	}
  }
}
