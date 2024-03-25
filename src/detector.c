/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <math.h>
#include <vector>

//#include "glib.h"
//#include "gtk/gtk.h"
//#include "gdk/gdk.h"

#include "declare.h"
#include "daten.h"
#include "detector.h"
#include "fit.h"
#include "error.h"
#include "round.h"


Tsingledetector::Tsingledetector(int open, int extra_open,double m0, int alastmin, TDaten *aparent)
{
double SNR;
std::cout<<"----"<<std::endl;	
  parent=aparent;
  n_open=open;
  n_extra_open=extra_open;
  jump=false;
  zero();
  greatest_h = 0;
  
  m = parent->a_fit.i_null + 
      parent->a_fit.i_channel * n_open + 
      parent->a_fit.i_extra_channel * n_extra_open;
    
  lastmin=alastmin;
  p = (m-m0) / 2;

  //std::cout<<"m: "<<m<<std::endl;	
  SNR= ((2*p)/(parent->a_fit.sigma)) * ((2*p)/(parent->a_fit.sigma));
  lambda=parent->schranke(round(0.5 + parent->ScTiRe/ SNR),filter_ord, fabs(p));
}

inline void 
Tsingledetector::zero()
{
  for (int i=0;i <= filter_ord;i++) h[i]=0.0;
}

inline void 
Tsingledetector::cumsum()
{
  for (int i=1;i <= filter_ord;i++) h[i]=h[i] + h[i-1];
}
 
Thohdetector::Thohdetector(TDaten  *aparent,
               int  &i, // Startpunkt
               int  &an_open, int &an_extra_open, // aktuelles Niveau
               double ai_null, double ai_channel, double ai_extra_channel, // Strom der Niveaus
               int  aSamples) // bis max. hierher
{
std::cout<<"----"<<std::endl;		
int   Samples;
double  i_null, i_channel, i_extra_channel;
unsigned int  z;
double  e;
bool up;
std::vector<Tsingledetector>::iterator iter;

  parent = aparent;
  detector = NULL;
  detectors.clear();
  n_open = an_open;
  n_extra_open = an_extra_open;
  i_null = ai_null;
  i_channel = ai_channel;
  i_extra_channel = ai_extra_channel;
  
  Samples = aSamples;
  m = i_null + i_channel * n_open + i_extra_channel * n_extra_open;
  jump = false;
  Create_N(i);
  while (i < Samples)
  {
    if ((i >= parent->Anz_Samples) || jump) break;
    z = parent->Wert[i];
    e = z - m;
    
    // PASCAL Kommentar
    // CalcHinkley(e, i);  // Testwerte der Detektoren berechnen
    // CheckJump(i);
    
    for (iter=detectors.begin(); iter != detectors.end(); iter++)
	{
      CalcH(*iter, i, e);
	}
    for (iter=detectors.begin(); iter != detectors.end(); iter++)
	{
      CheckJ(*iter);
	}

    if (jump)
    {
	  detector = &(detectors[new_model]); 	

	  // PASCAL kommentar
      //writeln('Sprung erkannt bei: ', i, 'lastmin: ', detector^.lastmin);
      //writeln('in das Niveau: ', detector^.m);
    
      i=detector->lastmin;               // beim letzen min. weitermachen
      if (detector->n_open > n_open)     // Sprung nach oben ?
      {
      	up = true;
      } else { 
        if (detector->n_open < n_open)
        {
          up = false;
        } else { 
          if (detector->n_extra_open > n_extra_open)
          {
            up = true;
          } else {
            up = false;
          }
        }
      }
      an_open = detector->n_open;                // nivea-nr. von level
      an_extra_open = detector->n_extra_open;    // und sub-level merken
      parent->flush_Sprung(_HOSHD,
                           parent->HOSHD_Spruenge,
                           i,
                           new_model,
                           up,
                           n_open, n_extra_open,
                           an_open, an_extra_open);
      detectors.clear(); // Detektoren l�schen, Speicher Freigeben
    }
    i++;
  }
}
/*
inline int 
operator==(const Tsingledetector& x, const Tsingledetector& y)
{
  return ((x.parent == y.parent) &&
	      (x.h[0] == y.h[0]) &&
	      (x.h[1] == y.h[1]) &&
	      (x.h[2] == y.h[2]) &&
	      (x.h[3] == y.h[3]) &&
	      (x.h[4] == y.h[4]) &&
	      (x.h[5] == y.h[5]) &&
	      (x.h[6] == y.h[6]) &&
	      (x.h[7] == y.h[7]) &&
	      (x.h[8] == y.h[8]) &&
          (x.greatest_h == y.greatest_h) &&
          (x.p == y.p) &&
          (x.m == y.m) &&
          (x.n_open == y.n_open) &&
          (x.n_extra_open == y.n_extra_open) &&
          (x.lastmin == y.lastmin) &&
          (x.lastmax == y.lastmax) &&
          (x.lambda == y.lambda) &&
		  (x.sigma == y.sigma) &&
          (x.jump == y.jump));
}
*/
int 
index(Tsingledetector &item, std::vector<Tsingledetector> &vec)
{
std::cout<<"----"<<std::endl;		
std::vector<Tsingledetector>::iterator iter;
int back = -1;
bool schleife=false;

    while (schleife == false)					//Tobias 2017 Schleife neu, original unten
      {
				back++;		
        if (&vec.at(back) == &item) schleife=true;
      }
  //for (iter=vec.begin();iter != vec.end();iter++)
  //{
  //
  //  if (&item == iter) break;    // item muss in liste sein // Tobias 2017: &item == iter
    // if (item == *iter) break; // item hat gleichen inhalt (braucht operator==)
  //  back++;
  //}
  return back;
}

inline int 
Thohdetector::Vorzeichen(double x)
{	
  //std::cout<<"item.lambda: "<<std::endl;	
  if (x >= 0) return 1; else return -1;
}

void 
Thohdetector::CalcH(Tsingledetector &item, int i, double e)
{
	std::cout<<"----"<<std::endl;	
  item.h[0] = Vorzeichen(item.p)*(e-item.p); // Hinkley Formel, p: HalfJumpMag.
  item.cumsum();
  if (item.h[1] <= 0) // Testwert Null setzen, wenn negativ
  {
    if (!item.jump) item.lastmin=i;
    item.zero();
  }
}

void
Thohdetector::CheckJ(Tsingledetector &item)
{
	std::cout<<"----"<<std::endl;	
  if ((item.h[filter_ord] > item.lambda) && !(item.jump))
  {

    //std::cout<<"item.lambda: "<<item.lambda<<std::endl;	
    new_model = index(item, detectors);

	jump = true;
    item.jump = true;
  }
}

void
Thohdetector::CalcHinkley(double e, int &i)
{
std::vector<Tsingledetector>::iterator iter;
std::cout<<"----"<<std::endl;		

  // Testwerte f�r alle Detektoren berechnen
  for (iter=detectors.begin(); iter != detectors.end(); iter++)
  {
    CalcH(*iter, i, e);
  }
}

void
Thohdetector::CheckJump(int &i)
{
std::vector<Tsingledetector>::iterator iter;
std::cout<<"----"<<std::endl;	
  for (iter=detectors.begin(); iter != detectors.end(); iter++)
  {
    CheckJ(*iter);
  }
}

void 
Thohdetector::Create_N(int &i)
{
  detectors.clear();

  if (n_open < parent->a_fit.n_channels)
  {
	Tsingledetector sd(n_open+1,0,m,i,parent);
    detectors.push_back(sd);
  }  // Detektor der Spr�nge nach oben erkennt
  if (n_open > 0)
  {
    Tsingledetector sd(n_open-1,0,m,i,parent);
	detectors.push_back(sd);
  }  // Detektor der Spr�nge nach unten erkennt
}

Tshdetector::Tshdetector(TDaten  *aparent,
                         int  &i,
                         int  &an_open, int &an_extra_open,
                         double ai_null, double ai_channel, double ai_extra_channel,
                         int  aSamples)
            : Thohdetector(aparent, i, an_open, an_extra_open,
                           ai_null, ai_channel, ai_extra_channel, aSamples)
{
// ruft nur Thohdetector Konstruktor auf 
}

void
Tshdetector::Create_N(int &i)
{
int   k, g;
bool allowed = false;
//Tsingledetector *pd; auch nicht in PASCAL

  diagnosis_mode = false;
  epsilon = 8*sqrt(parent->a_fit.sigma); // Schranke f�r erster Ordung, W-Rauschen, Bessel TP
  detectors.clear();

  for (k=0;k <= parent->a_fit.n_extra_channels;k++)
  {
    for (g=0;g <= parent->a_fit.n_channels;g++)
    {
      if ((k==n_extra_open) && (g==n_open))
      {
        allowed = false; //aktuelles Modell nicht erlaubt
      } else {
        switch (parent->candidate_restriction)
        {
          case 0:
            if ((abs((g+k) - (n_open+n_extra_open)) == 1 ) &&
                (abs(g - n_open ) <= 1) &&
                (abs(k - n_extra_open) <= 1))
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
            if ((k+g <= parent->a_fit.n_channels) &&
                (abs((g+k) - (n_open+n_extra_open)) <= 1) &&
                (abs(g - n_open) <= 1) &&
                (abs(k - n_extra_open) <= 1))
            {    
              allowed = true;
            } else {
              allowed = false;
            }  
            break;
          case 3:
            if ((abs((g+k) - (n_open+n_extra_open)) <= 1) &&
                (abs(g - n_open) <= 1) &&
                (abs(k - n_extra_open) <= 1))
            { 
              allowed = true;
            } else { 
              allowed = false;
            }  
            break;
          case 4:
            if ((abs((g+k) - (n_open+n_extra_open)) <= 2) &&
                (abs(g - n_open) <= 2) &&
                (abs(k - n_extra_open) <=2))
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
		Tsingledetector sd(g,k,m,i,parent);
        detectors.push_back(sd);
      }
    } 
  }
}

void
Tshdetector::CalcH(Tsingledetector &item, int i, double e)
{
  item.h[0] = (item.p)*(e-item.p);         // Hinkley Formel SHD, p: HalfJumpMag.
  item.h[1] = item.h[1] + item.h[0];
  if ((item.h[1] <=0) && !(diagnosis_mode))  // Testwert Null setzen, wenn negativ
  {
    item.lastmin=i;  // Position speichern
    item.zero();     // Null Setzen
  }
}

void
Tshdetector::CalcHinkley(double e, int &i)
{
std::vector<Tsingledetector>::iterator iter;

  // Testwerte f�r alle Detektoren berechnen
  for (iter=detectors.begin(); iter != detectors.end(); iter++)
  {
    CalcH(*iter, i, e);
  }
}

void 
Tshdetector::CheckJ(Tsingledetector &item)
{
  if (item.h[1] > epsilon)
  {
    diagnosis_mode = true;
    diagnosis_counter =0;
  }
}

void 
Tshdetector::Greatest_H(Tsingledetector &item, int &i, bool &positiv_h)
{
  if (item.h[1] > item.greatest_h)
  {
    item.greatest_h = item.h[1];
    item.lastmax = i;
  }
  if (item.h[1] > 0) positiv_h = true;
}

void 
Tshdetector::Decide(Tsingledetector &item, double hmax)
{
  if (item.greatest_h > hmax)
  {
    hmax = item.greatest_h;

    new_model = index(item, detectors);
  }
}

void
Tshdetector::CheckJump(int &i)
{
double  hmax;
bool positiv_h;
std::vector<Tsingledetector>::iterator iter;

  if (!diagnosis_mode)
  {
    for (iter=detectors.begin(); iter != detectors.end(); iter++)
    {
      CheckJ(*iter);
    }
  }
  new_model = -1;
  positiv_h = false;
  hmax = epsilon;
  if (diagnosis_mode)
  {
    diagnosis_counter++;
    for (iter=detectors.begin(); iter != detectors.end(); iter++)
    {
      Greatest_H(*iter, i, positiv_h);
    }
    if ((diagnosis_counter >= 7) || (positiv_h = false))
    {
      for (iter=detectors.begin(); iter != detectors.end(); iter++)
      {
        Decide(*iter, hmax);
      }
    }
    if (new_model != -1)
    {
      detector = &(detectors[new_model]);
      diagnosis_mode = false;  // Wichtig diagnosis_mode zur�cksetzen, sonst Endlosschleife
      diagnosis_counter = 0;
      if ((detector->lastmax - detector->lastmin) < 7)
      {
      // Dieser Spung z�hlt nicht
        i = detector->lastmax;
        Create_N(i);
      } else {
      // Dieser Sprung z�hlt
        jump = true;
      }
    }
  }
}

void
Tdetector::calc_sigma(Tsingledetector &item, int &t_res)
{
// Berechnung der s-Abweichung auch einfache Differenz denkbar
// Beginn der Berechung beim letzten Min. des jeweiligen Niv.

int  l;                                         
double m2;                                          

  m2 = 0;
  item.sigma = 999999999;
  if (item.jump)
  {
    for (l=item.lastmin+1;l <= item.lastmin+1+t_res;l++)
    {
      m2 = m2 + (item.m - parent->Wert[l]) * (item.m - parent->Wert[l]);
    }
    item.sigma = sqrt(m2/t_res+1);
    // writeln('Mittelwert: ', item^.m,' Sigma: ',item^.Sigma);
    // kein PASCAL kommentar
  }
}

void
Tdetector::FindMinSigma(Tsingledetector &item)
{
// Minimales Sigma bei gesetztem jump-flag

  if ((item.jump) && (sigma > item.sigma))
  {
    sigma = item.sigma;
    detector = &item;
  }
}

void 
Tdetector::FirstJump(Tsingledetector &item)
{
// Ersten erkannten Sprung ermitteln

  if ((item.jump) && (item.lastmin < detector->lastmin))
  {
    detector = &item;
  }
}

void
Tdetector::CheckTresTooLong(Tsingledetector &item, int &t_res)  
{
// Ueberpruefen, Ob Sprung zu weit entfernt

  if ((item.jump) && ((detector->lastmin-item.lastmin) > t_res))
  {
    item.jump = false;
  }
}

void
Tdetector::CalcHinkley(Tsingledetector &item, double e, int &i)
{
  item.h[0] = Vorzeichen(item.p)*(e-item.p); // Hinkley Formel, p: HalfJumpMag.
  item.cumsum();
  if (item.h[1] <= 0)                          // Testwert Null setzen, wenn negativ
  {
    if (!item.jump) item.lastmin = i;         // nur Speichern, falls noch kein Sprung erkannt
    item.zero();
  }
}

Tdetector::Tdetector(TDaten *aparent,
                     int &i,
                     int &ajump_pos,
                     int &an_open, int &an_extra_open)
{
std::cout<<"----"<<std::endl;		
int   Samples;
unsigned short  z;
double  e;
bool jump;
int   jump_zaehler;
bool up;
double  SNR;
int   t_res;
std::vector<Tsingledetector>::iterator iter;

  parent = aparent;
  // detectors := nil; // sollte eigentlich schon leer sein
  detector = NULL;
  n_open = an_open;
  n_extra_open = an_extra_open;
  jump_pos = ajump_pos;
  Samples = aparent->Settings.SamplesToProcess+aparent->Settings.SamplesToSkip;

  m = parent->a_fit.i_null + parent->a_fit.i_channel * n_open + parent->a_fit.i_extra_channel * n_extra_open;

  Create_Neighbors();
  jump = false;
  jump_zaehler = 0;
  if (parent->a_fit.n_extra_channels == 0)
  {
    SNR = parent->a_fit.i_channel/parent->a_fit.sigma;
  } else {
    SNR = parent->a_fit.i_extra_channel/parent->a_fit.sigma;
  }
  SNR = SNR*SNR;
  t_res = round(0.5 + parent->ScTiRe/SNR);
  while (i < Samples)
  {
    if ((i >= parent->Anz_Samples) || (jump_zaehler > (t_res))) break;
    z = parent->Wert[i];
    e = z - m;
    for (iter=detectors.begin(); iter != detectors.end(); iter++)
    {
      CalcHinkley(*iter, e, i);
    }
    
    // Alle Detektoren durchgehen, ob Sprung vorhanden
    int count=0;
    for (iter=detectors.begin(); iter != detectors.end(); iter++)
    {
      detector = &(*iter);
      if ((detector->h[filter_ord] > detector->lambda) && !(detector->jump))
      { 
        // Falls Testwert > Schranke dann Sigmaberechnung starten
        new_model = count; // das Niveau, da� Schranke �bersteigt merken
        jump = true;
        detector->jump = true;
      }
      count++;
    }
    if (jump) jump_zaehler++;
    if (jump_zaehler > (t_res))
    {
      detector = &(detectors[new_model]);     // irgendeine Annahme f�r detector
      for (iter=detectors.begin(); iter != detectors.end(); iter++)
      {
        FirstJump(*iter);                     // ersten Sprung in detector festhalten
      }
      for (iter=detectors.begin(); iter != detectors.end(); iter++)
      {
        CheckTresTooLong(*iter, t_res);       // Spr�nge verwerfen, die weiter als t_res vom ersten Sprung entfernt}
      }
      for (iter=detectors.begin(); iter != detectors.end(); iter++)
      {
        calc_sigma(*iter, t_res);             // Sigmaberechung
      }
      sigma = detector->sigma;
      for (iter=detectors.begin(); iter != detectors.end(); iter++)
      {
        FindMinSigma(*iter);                  // Das Niveau mit min. Sigma ist es
      }
      //PASCAL kommentar
      //detectors^.ForEach(@calc_Sigma);                           //Standartabweichung f�r alle Niveaus bereichen
      //sigma:=psingledetector((detectors^.At(new_model)))^.sigma; //erst mal irgendeine Annahme f�r Sigma machen
      //detector:=psingledetector(detectors^.At(new_model));       //kleinste Standartabweichung ermitteln
      //detectors^.ForEach(@FindMinSigma);                         //und das niveau in detector abspeichern
      //detectors^.ForEach(@CheckTresTooLong);  }

      // PASCAL mit index of new model
      new_model = index(*detector, detectors);

      // kein PASCAL kommentar
      // writeln('Sprung erkannt bei: ', i, 'lastmin: ', detector^.lastmin);
      // writeln('in das Niveau: ', detector^.m);
      i = detector->lastmin;                       // beim letzen min. weitermachen
      jump_pos = i;
      ajump_pos = jump_pos;                        // Jump Pos merken
      if (detector->n_open > n_open)               // Sprung nach oben ?
      {
      	up = true;
      } else {
      	if (detector->n_open < n_open)
        { 
          up = false;
        } else {
	  if (detector->n_extra_open > n_extra_open)
          {
            up = true;
          } else {
            up = false;
          }
        }
      }
      an_open = detector->n_open;                // nivea-nr. von level
      an_extra_open = detector->n_extra_open;    // und sub-level merken
      parent->flush_Sprung(_HOSHD,
                           parent->HOSHD_Spruenge,
                           i,
                           new_model,
                           up,
                           n_open, n_extra_open,
                           an_open, an_extra_open);
      detectors.clear();
    }
    i++;
  }
}


inline int 
Tdetector::Vorzeichen(double x)
{
  if (x >= 0) return 1; else return -1;
}

void
Tdetector::Create_Neighbors()
{
int  g,k;

  detectors.clear();

  //with parent^.a_fit
  for (k=0;k <= parent->a_fit.n_extra_channels;k++)
  {
    for (g=0;g <= parent->a_fit.n_channels;g++)
    {
      if (!((k==n_extra_open) && (g=n_open)))
      {
      	// PASCAL kommentar
        // if abs(g-n_open) <=1 then            nur Spr�nge um ein Niveau zugelassen, zum testen
        
		Tsingledetector sd(g,k,m,jump_pos,parent);
		detectors.push_back(sd);
      }
    }
  }
}
