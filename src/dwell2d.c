/*******************************************************************
 *  Patch Clamp 2d-dwell_time
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <fstream>
#include <math.h>
//#include <process.h>	2019
#include <errno.h>
#include <stdio.h>

//#include "glib.h"

#include "daten.h"
#include "gtkfunc.h"
#include "gtkapp.h"
#include "dwell.h"
#include "error.h"

void 
ZToRGB(const double z, const double zmin, const double zmax, 
	   unsigned int &R, unsigned int &G, unsigned int &B, bool grey = false)
{
// hier wird einen z-wert in einen intervall eine farbe zugeordnet (RGB tripel)
double Z, Red, Green, Blue;

  // zuerst intervall auf [0:1] abbilden
  if (zmax - zmin != 0) 
  {  
    Z = (z - zmin) / (zmax - zmin); 
  } else {
    Z = 0.5;
  }

  if (grey)
  {
    // grau zuordnen
    Red = Z;
    Green = Z;
    Blue = Z;
  } else {
  // farbe zuordnen
  if (Z < 0)
  {
     Red = 0;
	 Green = 0;
	 Blue = 0;
  } else {
    if (Z < 0.2)
    {
      Red = 0;
      Green = 0;
      Blue = 1 + (Z - 0.2) * 5;
    } else {
      if (Z < 0.4)
	  {
        Red = 0;
        Green = 1 + (Z - 0.4) * 5;
        Blue = 1;
	  } else {
        if (Z < 0.6)
		{
          Red = 0;
          Green = 1;
          Blue = -(Z - 0.6) * 5;
		} else {
          if (Z < 0.8)
		  {
            Red = 1 + (Z - 0.8) * 5;
            Green = 1;
            Blue = 0;
		  } else {
		    if (Z <= 1)
			{
              Red = 1;
              Green = -(Z - 1) * 5;
              Blue = 0;
			} else { 
		      Red = 1;
              Green = 0;
              Blue = 0;
			}
		  }
		}
	  }
    } 
  }
  }
  // RGB im bereich [0:255]
  R = (unsigned int)(Red * 255);
  G = (unsigned int)(Green * 255);
  B = (unsigned int)(Blue * 255);
}

inline unsigned short
TDwell_lb::bins_open()
{
  return bins_per_log() * (max_open_log() - min_open_log());
}

inline unsigned short
TDwell_lb::bins_close()
{
  return bins_per_log() * (max_close_log() - min_close_log());
}

inline char*
TDwell_lb::name()
{
  return dwell_name;
}

inline short
TDwell_lb::min_close_log()
{
  return dwell_min_close_log;
}

inline short
TDwell_lb::max_close_log()
{
  return dwell_max_close_log;
}

inline short
TDwell_lb::min_open_log()
{
  return dwell_min_open_log;
}

inline short
TDwell_lb::max_open_log()
{
  return dwell_max_open_log;
}

inline double
TDwell_lb::min_e(bool sym )	//Tobias 2017 weg. = true
{ 
  // !autorange?
  if (!(e_auto_range()))
  {
    return dwell_min_e;
  } else {
	if (sym)
	{
	  double min_gd = min();
	  double max_gd = max();
	  return -std::max(abs(min_gd), abs(max_gd));
	} else {
      return min();
	}
  }
}

inline double
TDwell_lb::max_e(bool sym)		//Tobias 2017 weg. = true
{
  // !autorange?
  if (!(e_auto_range()))
  {
    return dwell_max_e;
  } else {
	if (sym)
	{
	  double min_gd = min();
	  double max_gd = max();
	  return std::max(abs(min_gd), abs(max_gd));
	} else {
      return max();
	}
  }
}

inline unsigned short
TDwell_lb::bins_per_log()
{
  return dwell_bins_per_log;
}

inline void 
TDwell_lb::name(char* aname)
{
  strcpy(dwell_name, aname);
}

inline void
TDwell_lb::close_range(short amin, short amax)
{
  dwell_min_close_log = amin;
  dwell_max_close_log = amax;
}

inline void
TDwell_lb::open_range(short amin, short amax)
{
  dwell_min_open_log = amin;
  dwell_max_open_log = amax;
}

inline void
TDwell_lb::e_range(double amin, double amax)
{
  dwell_min_e = amin;
  dwell_max_e = amax;
}

inline void
TDwell_lb::bins_per_log(unsigned short abins)
{
  dwell_bins_per_log = abins;
}

inline bool
TDwell_lb::e_auto_range()
{	
  return (dwell_min_e >= dwell_max_e);
}

inline bool
TDwell_lb::sqrt_e()
{
  return dwell_sqrt_e;
}

inline void
TDwell_lb::sqrt_e(bool asqrt)
{
  dwell_sqrt_e = asqrt;
}

void 
TDwell_lb::calcrange(short &log_min_close, short &log_max_close,
					 short &log_min_open, short &log_max_open)
{
// muesste eigentlich zu Spruengen!

  double min_c = 0;
  double max_c = 0;
  double min_o = 0;
  double max_o = 0;
  unsigned short C1, C2;

  if (!Daten.HOHD_Spruenge.empty())
  {
    // min,max laenge ermitteln
    std::vector<TSprung>::iterator iter1;
    std::vector<TSprung>::iterator iter2;
    double laenge;
	  
	iter2 = Daten.HOHD_Spruenge.begin();
        iter1 = iter2++; 
	  
	if (iter2 != Daten.HOHD_Spruenge.end())    // mehr als nur ein sprung
	{
	  // zuerst laenge von SP1 - SP 2
      C1 = (*iter2).C_1();
	  C2 = (*iter2).C_2();
	  if (C1 < C2)
	  {
	    min_c = max_c = ((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval;
	  } else {
	    min_o = max_o = ((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval;
	  }
		  
	  iter1 = iter2++;
		  
	  if (iter2 != Daten.HOHD_Spruenge.end()) // mehr als nur zwei spruenge
	  {
	    // und SP2 - SP3 ermitteln
        C1 = (*iter2).C_1();
        C2 = (*iter2).C_2();
        if (C1 < C2)
		{
		  min_c = max_c = ((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval;
		} else {
		  min_o = max_o = ((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval;
		}
			  
		// dann rest durchlaufen
		while (iter2 != Daten.HOHD_Spruenge.end())
		{
          C1 = (*iter2).C_1();
	      C2 = (*iter2).C_2();
		  if (C1 < C2)
		  {
		    laenge = ((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval; // in us
		    if (laenge > max_c) max_c = laenge;
		    if (laenge < min_c) min_c = laenge;
		  } else {
		    laenge = ((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval; // in us
		    if (laenge > max_o) max_o = laenge;
		    if (laenge < min_o) min_o = laenge;
		  }
		  iter1 = iter2++; // beide einen sprung weiter
		}
	  }
	}
    if (max_c != 0) log_max_close = (short)ceil(log10(max_c));
	if (min_c != 0) log_min_close = (short)floor(log10(min_c));
	if (max_o != 0) log_max_open = (short)ceil(log10(max_o));
	if (min_o != 0) log_min_open = (short)floor(log10(min_o));
  } else {
    log_max_close = log_min_close = 0;
	log_max_open = log_min_open = 0;
  }
}

void
TDwell_lb::init(char* aname, unsigned short abins)	//Tobias 2017 weg: =0
{

short log_min_close = Daten.Fitparameter.log_min_close;	//Tobias 2021
short log_max_close = Daten.Fitparameter.log_max_close;	//Tobias 2021
short log_min_open = Daten.Fitparameter.log_min_open;		//Tobias 2021
short log_max_open = Daten.Fitparameter.log_max_open;		//Tobias 2021
abins = Daten.Fitparameter.bins_per_log;								//Tobias 2021
  // wenn abins == 0 behalte alte anzahl bins!
  //if (abins != 0) bins_per_log(abins);        

  //calcrange(log_min_close, log_max_close, log_min_open, log_max_open);	//Tobias 2021
  init(aname, log_min_close, log_max_close, log_min_open, log_max_open, abins);
}

void
TDwell_lb::init(char* aname,
				short log_min_close, short log_max_close,
				short log_min_open, short log_max_open,
				unsigned short abins)													//Tobias 2017 weg: =0
{
  // wenn abins == 0 behalte alte anzahl bins!
  if (abins != 0) bins_per_log(abins);

  close_range(log_min_close, log_max_close);
  open_range(log_min_open, log_max_open);
  e_range(0, 0);
  name(aname);
}

void
TDwell_lb::GNUplot(bool dosqrt)	//Tobias 2017 weg: =true
{
/*2019
  char cmd[] = "wgnuplot";
  char arg[] = "\"gnuplot.tmp\"";
  errno = 0;
  short returnvalue = _spawnlp (_P_NOWAIT, cmd, cmd, arg, NULL); 

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
TDwell_lb::out(std::ostream& out)
{
// wird von operatoren << fuer 1d und 2d benutzt
	
  // abspeichern ascii
  out << "#dwell_lb"     << std::endl
	  << name()          << std::endl
      << min_close_log() << "   # min_close"    << std::endl
      << max_close_log() << "   # max_close"    << std::endl
      << min_open_log()  << "   # min_open"     << std::endl
      << max_open_log()  << "   # max_open"     << std::endl
	  << min_e(false)    << "   # min_e"        << std::endl
	  << max_e(false)    << "   # max_e"        << std::endl
      << bins_per_log()  << "   # bins per log" << std::endl
      << sqrt_e()        << "   # sqrt_e"       << std::endl;
}

void
TDwell_lb::in(std::istream& in)
{
char     buffer[MaxTextLen];
double  buf_gd, buf2_gd;
short   buf_i16, buf2_i16;
unsigned short  buf_ui16;
bool boolvar;
char     c; 

  in.get(buffer, MaxTextLen, '\n');
  if (!strncmp(buffer, "#dwell_lb",9))			//2019 removed strlwr
  {
	while (in.get(c)) if (c == '\n') break;  // kann laenger als Max sein!
	in.get(buffer, MaxTextLen, '\n');
	name(buffer);
	while (in.get(c)) if (c == '\n') break;
    in >> buf_i16;
	while (in.get(c)) if (c == '\n') break;
    in >> buf2_i16;
	close_range(buf_i16, buf2_i16);
	while (in.get(c)) if (c == '\n') break;
    in >> buf_i16;
	while (in.get(c)) if (c == '\n') break;
    in >> buf2_i16;
	open_range(buf_i16, buf2_i16);
	while (in.get(c)) if (c == '\n') break;
    in >> buf_gd;
	while (in.get(c)) if (c == '\n') break;
    in >> buf2_gd;
	e_range(buf_gd, buf2_gd);
	while (in.get(c)) if (c == '\n') break;
    in >> buf_ui16;
	bins_per_log(buf_ui16);;
	while (in.get(c)) if (c == '\n') break;
    in >> boolvar;
	sqrt_e(boolvar);
	while (in.get(c)) if (c == '\n') break;
  } else {
    in.clear(std::ios::badbit); // status setzen
  }
}

void 
TDwell_lb::load(std::istream& in)
{
// zum laden aus setfile
int     l;
char     buffer[MaxTextLen];
double  buf_gd, buf2_gd;
short   buf_i16, buf2_i16;
unsigned short  buf_ui16;
bool boolvar;

  in.read((char *)&l, sizeof(l));
  in.read((char *)&buffer, l * sizeof(char));
  name(buffer);

  in.read((char *)&buf_i16, sizeof(buf_i16));
  in.read((char *)&buf2_i16, sizeof(buf2_i16));
  close_range(buf_i16, buf2_i16);
  in.read((char *)&buf_i16, sizeof(buf_i16));
  in.read((char *)&buf2_i16, sizeof(buf2_i16));
  open_range(buf_i16, buf2_i16);
  in.read((char *)&buf_gd, sizeof(buf_gd));
  in.read((char *)&buf2_gd, sizeof(buf2_gd));
  e_range(buf_gd, buf2_gd);
  in.read((char *)&buf_ui16, sizeof(buf_ui16));
  bins_per_log(buf_ui16);
  in.read((char *)&boolvar, sizeof(boolvar));
  sqrt_e(boolvar);
}

void 
TDwell_lb::store(std::ostream& out)
{
// abspeichern ins jobfile
int     l;
char     buffer[MaxTextLen];
double  buf_gd;
short   buf_i16;
unsigned short  buf_ui16;
bool boolvar;

  strcpy(buffer, name());
  l = strlen(buffer)+1; // mit term 0
  out.write((char *)&l, sizeof(l));
  out.write((char *)&buffer, l * sizeof(char));

  buf_i16 = min_close_log();
  out.write((char *)&buf_i16, sizeof(buf_i16));
  buf_i16 = max_close_log();
  out.write((char *)&buf_i16, sizeof(buf_i16));
  buf_i16 = min_open_log();
  out.write((char *)&buf_i16, sizeof(buf_i16));
  buf_i16 = max_open_log();
  out.write((char *)&buf_i16, sizeof(buf_i16));
  buf_gd = min_e(false);
  out.write((char *)&buf_gd, sizeof(buf_gd));
  buf_gd = max_e(false);
  out.write((char *)&buf_gd, sizeof(buf_gd));
  buf_ui16 = bins_per_log();
  out.write((char *)&buf_ui16, sizeof(buf_ui16));
  boolvar = sqrt_e();
  out.write((char *)&boolvar, sizeof(boolvar));
}

TDwell_1d::TDwell_1d(unsigned short abins )		//Tobias weg: =50
{
  sqrt_e(true);
  init((char *)" ", abins);
}

void
TDwell_1d::init(char* aname, unsigned short abins)	//Tobias weg: =0
{
  TDwell_lb::init(aname, abins);
  if (!Daten.HOHD_Spruenge.empty())
  {
    // vektoren auf n\F6tige gr\F6sse bringen und loeschen!
	Open.clear();
	Close.clear();
	Open.resize(bins_open(), 0);
    Close.resize(bins_close(), 0);
  }
}

void
TDwell_1d::init(char* aname,
				short log_min_close, short log_max_close,
				short log_min_open, short log_max_open,
				unsigned short abins)					//Tobias weg: =0
{
  TDwell_lb::init(aname,
	              log_min_close, log_max_close, 
	              log_min_open, log_max_open,
				  abins);

  // vektoren auf n\F6tige gr\F6sse bringen und loeschen!
  Open.clear();
  Close.clear();
  Open.resize(bins_open(), 0);
  Close.resize(bins_close(), 0);
}

void
TDwell_1d::calc(std::vector<TSprung>* vec_ptr, double maxevents)		//Tobias weg: =0
{
std::vector<TSprung>::iterator iter1;
std::vector<TSprung>::iterator iter2;
double log_laenge; //log10(sprunglaenge in us)
unsigned short i;
unsigned short  c1, c2;

  if (vec_ptr->empty()) return;
 
  // schon erfolgreich initialisiert?
  if ((min_close_log() == max_close_log()) || (min_open_log() == max_open_log())) init(name());

  // gehe zu zweitem Sprung
  iter2 = vec_ptr->begin();
  iter1 = iter2++;  

  
  // bis SP2 am ende immer alle einen sprung weiter
  while (iter2 != vec_ptr->end())
  {
	// SP1 - laenge1 - SP2
    log_laenge = log10(((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval);  // log(l in us)

    c1 = (*iter2).C_1();
    c2 = (*iter2).C_2();
/*TOBIAS
	if (warn && ((c1 > 1) || (c2 > 1))) 
	{
	  warning("multichannel TDwell_1d::calc");
      warn = false;                             // nur eine warning!
	}
*/
    if (c1 < c2)
	{
	  // jump up, nach close intervall
      // bin ermitteln!
	  i = (unsigned short)((log_laenge - min_close_log()) * bins_close() / (double)(max_close_log() - min_close_log()));
      if (i < bins_close()) Close[i]++;
	} else {
	  // jump down, nach open intervall
      // bin ermitteln!
      i = (unsigned short)((log_laenge - min_open_log()) * bins_open() / (double)(max_open_log() - min_open_log()));
	  if (i < bins_open()) Open[i]++;
	}
    iter1++;
	iter2++; // alle einen sprung weiter
  }
}

TDwell_1d& 
TDwell_1d::operator+=(TDwell_1d& B)
{
  // noch nicht erfolgreich initialisiert?
  if ((min_close_log() == max_close_log()) || (min_open_log() == max_open_log()))
  { 
    sqrt_e(B.sqrt_e());
	init(B.name(), 
		 B.min_close_log(), B.max_close_log(),
		 B.min_open_log(), B.max_open_log(),
		 B.bins_per_log());
  }
  // sonst bleiben die grenzen von this gueltig!

  // sind diese auch kompatibel?
  if (bins_per_log() != B.bins_per_log())
  {
    warning((char *)"bins_per_log are different!");
    return *this;
  }

  // beschraenken auf schnittmenge

  // grenzen der schnittmenge
  short log_min_close = std::max(min_close_log(), B.min_close_log());
  short log_max_close = std::min(max_close_log(), B.max_close_log());
  short log_min_open  = std::max(min_open_log(), B.min_open_log());
  short log_max_open  = std::min(max_open_log(), B.max_open_log());

  // offsets
  short close_off = (log_min_close - min_close_log()) * bins_per_log();
  short open_off  = (log_min_open - min_open_log()) * bins_per_log();
  short B_close_off = (log_min_close - B.min_close_log()) * bins_per_log();
  short B_open_off  = (log_min_open - B.min_open_log()) * bins_per_log();

  // summen berechnen 
  for (unsigned short j=0;j < (log_max_close - log_min_close) * bins_per_log();j++)
  {
    Close[j+close_off] += B.Close[j+B_close_off]; 
  }  

  for (unsigned short j=0;j < (log_max_open - log_min_open) * bins_per_log();j++)
  {
	Open[j+open_off] += B.Open[j+B_open_off]; 
  }  

  return *this;
}

void 
TDwell_1d::GNUplot(bool dosqrt)			//Tobias 2017 weg: = true
{
double log_laenge;

  std::ofstream outO("open.tmp");

  // mit 0 anfangen
  outO << min_open_log() << " " << 0 << std::endl;

  for (unsigned short i=0;i < bins_open();i++)
  {
	log_laenge = ((i * (max_open_log() - min_open_log()) / (double)bins_open() + min_open_log()) +
	              ((i+1) * (max_open_log() - min_open_log()) / (double)bins_open() + min_open_log())) / 2;
 
    outO << log_laenge << " " << Open[i] << std::endl;
  }
  // mit 0 aufhoeren
  outO << max_open_log() << " " << 0 << std::endl;


  std::ofstream outC("close.tmp");

  // mit 0 anfangen
  outC << min_close_log() << " " << 0 << std::endl;
  for (unsigned short i=0;i < bins_close();i++)
  {
	log_laenge = ((i * (max_close_log() - min_close_log()) / (double)bins_close() + min_close_log()) +
	              ((i+1) * (max_close_log() - min_close_log()) / (double)bins_close() + min_close_log())) / 2;

    outC << log_laenge << " " << Close[i] << std::endl;
  }
  // mit 0 aufhoeren
  outC << max_close_log() << " " << 0 << std::endl;

  // run gnuplot
  std::ofstream out("gnuplot.tmp");
  out << "set nologscale xy" << std::endl
      << "set data style points" << std::endl;
  
  // closed
  out << "set xlabel \"log(shut time\\)\" 0,0" << std::endl;

  if (dosqrt)
  {
    out << "set ylabel \"sqrt(Number of events)\" 0,0" << std::endl;
  } else {
    out << "set ylabel \"Number of events\" 0,0" << std::endl;
  }

  out << "plot 'close.tmp' ";

  if (dosqrt)
  {
    out << "using ($1):(sqrt($2)) ";
  } else {
    out << "using ($1):($2) ";
  }

  out << "title \"" << name() << "\", \\" << std::endl
  	  << "     'close.tmp' ";

  if (dosqrt)
  {
    out << "using ($1):(sqrt($2)) ";
  } else {
    out << "using ($1):($2) ";
  }

  out << "smooth sbezier title \"smoothed\"" << std::endl
      << "pause -1" << std::endl;

  // open
  out << "set xlabel \"log(open time\\)\" 0,0" << std::endl
	  << "plot 'open.tmp' ";
	  
  if (dosqrt)
  {
    out << "using ($1):(sqrt($2)) ";
  } else {
    out << "using ($1):($2) ";
  }
	  
  out << "title \"" << name() << "\", \\" << std::endl
	  << "     'open.tmp' ";
	  
  if (dosqrt)
  {
    out << "using ($1):(sqrt($2)) ";
  } else {
    out << "using ($1):($2) ";
  }
	  	  
  out<< "smooth sbezier title \"smoothed\"" << std::endl
      << "pause -1" << std::endl;

  TDwell_lb::GNUplot();
}

double 
TDwell_1d::min()
{
  return 0.0;  // erst sinnvoll, wenn Open und Close getrennt
}

double 
TDwell_1d::max()
{
  return 0.0;  // erst sinnvoll, wenn Open und Close getrennt
}


std::ostream& operator<<(std::ostream& out, TDwell_1d& Dwell)
{
unsigned short i;

  Dwell.out(out);

  out << "# Close:" << std::endl;
  for (i=0;i < Dwell.bins_close();i++)
  {
    out << Dwell.Close[i] << std::endl;
  }
  if (i == 0) out << std::endl;

  out << "# Open:" << std::endl;
  for (i=0;i < Dwell.bins_open();i++)
  {
    out << Dwell.Open[i] << std::endl;
  }
  if (i == 0) out << std::endl;

  return out;
}

std::istream& operator>>(std::istream& in, TDwell_1d& Dwell)
{
char buffer[MaxTextLen];
char c;

  Dwell.in(in);

  in.get(buffer, MaxTextLen, '\n');

  if (!strncmp(buffer, "# close:",8))			//2019 removed strlwr
  {
    Dwell.init(Dwell.name(),
		       Dwell.min_close_log(), Dwell.max_close_log(),
		       Dwell.min_open_log(), Dwell.max_open_log(), 
			   Dwell.bins_per_log());
	  
	while (in.get(c)) if (c == '\n') break;
    
	//close einlesen
	for (unsigned short i=0;i < Dwell.bins_close();i++)
	{
      in >> Dwell.Close[i];
	}


	while (in.get(c)) if (c == '\n') break;

    in.get(buffer, MaxTextLen, '\n');
    
	if (!strncmp(buffer, "# open:",7))		//2019 removed strlwr
    { 
	  while (in.get(c)) if (c == '\n') break;

	  // open einlesen
      for (unsigned short i=0;i < Dwell.bins_open();i++)
	  {
        in >> Dwell.Open[i];
	  }
      while (in.get(c)) if (c == '\n') break;
    } else {
     in.clear(std::ios::badbit); // status setzen
    }
  } else {
    in.clear(std::ios::badbit); // status setzen
  }
  return in;
}

void 
TDwell_1d::load(std::istream& in)
{
// zum laden aus setfile
  TDwell_lb::load(in);

  init(name(),
       min_close_log(), max_close_log(),
	   min_open_log(), max_open_log(), 
	   bins_per_log());

  for (unsigned short i=0;i < bins_close();i++)
  {
    in.read((char *)&Close[i], sizeof(Close[i]));
  }

  for (unsigned short i=0;i < bins_open();i++)
  {
    in.read((char *)&Open[i], sizeof(Open[i]));
  }
}

void 
TDwell_1d::store(std::ostream& out)
{
// abspeichern ins jobfile
  TDwell_lb::store(out);

  for (unsigned short i=0;i < bins_close();i++)
  {
    out.write((char *)&Close[i], sizeof(Close[i]));
  }

  for (unsigned short i=0;i < bins_open();i++)
  {
    out.write((char *)&Open[i], sizeof(Open[i]));
  }
}

TDwell_2d::TDwell_2d(unsigned short abins)		//Tobias 2017 weg: = 10
{
  sqrt_e(true);
  init((char *)" ", abins);
}

void
TDwell_2d::init(char* aname, unsigned short abins )	//Tobias 2017 weg: = 0
{

  TDwell_lb::init(aname, abins);

  if (!Daten.HOHD_Spruenge.empty())
  {
	// matrix auf n\F6tige gr\F6sse bringen und loeschen!
    CO.clear();
    CO.resize(bins_close(), bins_open(), 0);
    COD.clear();
    COD.resize(bins_close(), bins_open(), 0);
 
	
  }
}

void
TDwell_2d::init(char* aname,
				short log_min_close, short log_max_close,
				short log_min_open, short log_max_open,
				unsigned short abins)		//Tobias 2017 weg: = 0
{
  TDwell_lb::init(aname,
	              log_min_close, log_max_close, 
	              log_min_open, log_max_open,
				  abins);
  // matrix auf n\F6tige gr\F6sse bringen und loeschen!
  CO.clear();
  CO.resize(bins_close(), bins_open(), 0);
  
  COD.clear();
  COD.resize(bins_close(), bins_open(), 0);
}

void
TDwell_2d::calc(std::vector<TSprung>* vec_ptr, double maxevents)		//Tobias 2017 weg: = 0
{
std::vector<TSprung>::iterator iter1;
std::vector<TSprung>::iterator iter2;
std::vector<TSprung>::iterator iter3;
double log_laenge1, log_laenge2;
unsigned short  i, j;
unsigned short   c1, c2;


  if (vec_ptr->empty()) return;
  // schon erfolgreich initialisiert?
  
  if ((min_close_log() == max_close_log()) || (min_open_log() == max_open_log())) 
      init(name());
   
  // Tobias 2021 reordered -> jump list 1./2./3. position 
  iter1 = vec_ptr->begin();
  iter2 = vec_ptr->begin();
  iter3 = vec_ptr->begin();
  if (iter1 == vec_ptr->end())
  	return;
  iter2++;
  if (iter2 == vec_ptr->end())
  	return;
  iter3++;
  iter3++;  
  if (iter3 == vec_ptr->end())
  	return;

  if (iter2 == vec_ptr->end())  return;
  
  while ((*iter2).position() == 0)		//Tobias 2021 to get to first event
   {
   	iter1++;
   	iter2++;
   	iter3++;
   }	

 
  double nr_events = events(); // Schattenkonto fuehren!
  // bis SP3 am ende immer alle einen sprung weiter
  while (iter3 != vec_ptr->end())
  {
	// SP1 - laenge1 - SP2
    log_laenge1 = log10(((*iter2).position() - (*iter1).position()) * Daten.Settings.Samplininterval);  // log (l in us)
	// SP2 - laenge2 - SP3
    log_laenge2 = log10(((*iter3).position() - (*iter2).position()) * Daten.Settings.Samplininterval);  // log (l in us)
    //bin ermitteln!
    i = (unsigned short)((log_laenge1 - min_close_log()) * bins_close() / (double)(max_close_log() - min_close_log()));
    j = (unsigned short)((log_laenge2 - min_open_log()) * bins_open() / (double)(max_open_log() - min_open_log()));

	// nur wenn i,j im richtigen intervall liegen,
     
    if ((i < bins_close()) && (j < bins_open()))
	{
	  c1 = (*iter2).C_1();
	  c2 = (*iter2).C_2();
/*Tobias	
      if (warn && ((c1 > 1) || (c2 > 1)))
	  {
         warning((char *)"multichannel TDwell_2d::calc");
         warn = false;  // nur eine warnung
	  }
*/
	  if (c1 < c2)
	     {
	      // zum Ausgeben der jumps, als bin und als Punktl\E4enge
	      //std::cout<<"i: "<<i<<" l\E4nge1: "<<((*iter2).position() - (*iter1).position())<<std::endl;
        //std::cout<<"j: "<<j<<" l\E4nge2: "<<((*iter3).position() - (*iter2).position())<<std::endl;
              
	      // Close - Open
	      CO[i][j]++;
	      nr_events++;
	     } 
	  else 
	     {
              // Open - Close
        CO[j][i]++;							//2021 Tobias
	      nr_events++;      
       }
	}
	iter1++;
  iter2++;
	iter3++; // alle einen sprung weiter
	if ((maxevents > 0) && (nr_events >= maxevents))
	{
	  //restliche spruenge loeschen fuer 1d calc!
	  vec_ptr->erase(iter3, vec_ptr->end());
      break;
	}
  }

}

//Tobias: Procedur zum berechen der Dependency-Plots nach Megleby & Song*************************************************************************************************************
void 
TDwell_2d::calc_dependency()
{
 int  i, j, ii, jj;        	//Laufvariablen
 double Nobs, Nind;   		// aus Moss & Magleby
 double sumclose, sumopen;	// Summe aller geschlossenen bzw. offenen Ereignisse (aus der 2D-Matrix)
 double columnclose, rowopen;	// Zahl der Erreignisse in einer Spalte, bzw. Zeile
 
//Initialisierung  
 sumclose = 0;  
 sumopen  = 0; 

Daten.Dwell_2d_diff.init((char *)"dependency",
                           Daten.Fitparameter.log_min_close,
	                   Daten.Fitparameter.log_max_close,
	                   Daten.Fitparameter.log_min_open,
	                   Daten.Fitparameter.log_max_open,
			   Daten.Fitparameter.bins_per_log);	

// Berechnung der Open und Close  Ereignisse gesamt
 for (i=0; i < ((Daten.Fitparameter.log_max_close - Daten.Fitparameter.log_min_close) * Daten.Fitparameter.bins_per_log); i++)
  for (j=0; j <((Daten.Fitparameter.log_max_open - Daten.Fitparameter.log_min_open) * Daten.Fitparameter.bins_per_log); j++) 
  {
   sumclose = sumclose + CO[i][j];
  }	 
 sumopen = sumclose; //identisch 
 std::cout<<"sumopen:"<<sumopen<<std::endl;
 
 for (i=0; i < ((Daten.Fitparameter.log_max_close - Daten.Fitparameter.log_min_close) * Daten.Fitparameter.bins_per_log); i++)
  {
  std::cout<<std::endl;
   for (j=0; j < ((Daten.Fitparameter.log_max_open - Daten.Fitparameter.log_min_open) * Daten.Fitparameter.bins_per_log); j++) 
    {
     columnclose = 0;	
     rowopen     = 0;	
     for (ii=0; ii < ((Daten.Fitparameter.log_max_close - Daten.Fitparameter.log_min_close) * Daten.Fitparameter.bins_per_log); ii++) 	
       columnclose = columnclose + CO[ii][j];     
   
     for (jj=0; jj < ((Daten.Fitparameter.log_max_open - Daten.Fitparameter.log_min_open) * Daten.Fitparameter.bins_per_log); jj++)  
       rowopen = rowopen + CO[i][jj]; 

     //std::cout<<i<<":"<<j<<" columnclose:"<<columnclose<<" rowopen:"<<rowopen<<std::endl;

     Nind = (columnclose / sumclose) * (rowopen / sumopen) * sumopen;	
     Nobs = CO[i][j]; 	
     //std::cout<<"Nobs:"<<Nobs<<" Nind:"<<Nind<<std::endl;;  
     if (Nind == 0)
        Daten.Dwell_2d_diff.CO[i][j] = 0;
     else  
        Daten.Dwell_2d_diff.CO[i][j] =   ((Nobs-Nind) / Nind);  //

if (Daten.Dwell_2d_diff.CO[i][j] > 2 )
    Daten.Dwell_2d_diff.CO[i][j] = 2;
if (Daten.Dwell_2d_diff.CO[i][j] < -2 )
    Daten.Dwell_2d_diff.CO[i][j] = -2;

//     if (Daten.Dwell_2d_diff.CO[i][j] > 0)
//       Daten.Dwell_2d_diff.CO[i][j] = log(Daten.Dwell_2d_diff.CO[i][j]);
//     else
//       if (Daten.Dwell_2d_diff.CO[i][j] < 0)
//          Daten.Dwell_2d_diff.CO[i][j] = -log(-Daten.Dwell_2d_diff.CO[i][j]);  

     //std::cout<<i<<":"<<j<<":"<<Daten.Dwell_2d_diff.CO[i][j]<<std::endl;
     //std::cout<<"Nobs: "<<Nobs<<" Nind: "<<Nind<<" COD:"<<COD[i][j]<<":i:"<<i<<":j:"<<j<<std::endl; 	
  }
 } 
    
 
//Zeichnen des 2D-Histogramms-Dependency-plots;
Daten.graphtype = dwell_diff;
Daten.paint();

// Setzen des Flags, damit das diff-Bild gezeichnet wird, obwohl die Matrix B leer ist! 
Daten.dependency_plot_draw_diff = true;
Application.set_menu_state();
}	

//*******************************************************************************************************************




void
TDwell_2d::addjumps(TDwell_1d* dwell_1d , double maxevents)		//Tobias 2017 weg: =NULL  = 0
{
std::vector<TSprung>::iterator  sprung_iter;   
unsigned short channels = 0;
unsigned short subchannels = 0;

  for(sprung_iter=Daten.HOHD_Spruenge.begin();sprung_iter != Daten.HOHD_Spruenge.end();sprung_iter++) 
  {
    if ((*sprung_iter).C_2() > channels) channels = (*sprung_iter).C_2();
    if ((*sprung_iter).U_2() > subchannels) subchannels = (*sprung_iter).U_2();
  }
  calc(&Daten.HOHD_Spruenge, maxevents);
  if (dwell_1d != NULL) dwell_1d->calc(&Daten.HOHD_Spruenge, maxevents);




/* Tobias:funktioniert f\FCr mehrere Kan\E4le nicht!
  if (channels == 1)
  {
    calc(&Daten.HOHD_Spruenge, maxevents);
    if (dwell_1d != NULL) dwell_1d->calc(&Daten.HOHD_Spruenge, maxevents);
  } else {
    Daten.separatechannels();
    for(unsigned short i=0;i < channels;i++)
	{
      calc(&Daten.SingleChannels[i], maxevents);
      if (dwell_1d != NULL) dwell_1d->calc(&Daten.SingleChannels[i], maxevents);
	}
  }
*/
  return;
}
//***********************************************************************************************************************************
//Tobias:Procedur zum Schreiben eines Log-Files, das die Sprungl\E4ngen sortiert nach Abfolge bei potentialgesteuerten Zeitreihen enth\E4lt

void
TDwell_2d::Hinkley_Log(TDwell_1d* dwell_1d, double maxevents)		//Tobias 2017 weg: = NULL  = 0
{
 std::vector<TSprung>::iterator  sprung_iter;   
 int i,j;
 int pos1,pos2,sprunglaenge;	
 
//initialisierung
  pos1=0;
  pos2=0;
  i=0;     
//in die Bins sortieren       
  for(sprung_iter=Daten.HOHD_Spruenge.begin();sprung_iter != Daten.HOHD_Spruenge.end();sprung_iter++) 
   {
    i++;	
    pos2 = (*sprung_iter).sprung_position;
    
    if( (i>2) && (i<=12) )
      {	
       sprunglaenge = pos2-pos1;
       //std::cout<<"i:"<<i-2<<" Sprung:"<<sprunglaenge<<std::endl;	
       j = 0;
       while ( log(sprunglaenge) > j * Daten.logx)
         {
          j++;	
         }	
       Daten.dwellmatrix[i-2][j]++;  
      } 
    pos1=pos2;
   }
   
  return;
}

//***********************************************************************************************************************************




TDwell_2d& 
TDwell_2d::operator+=(TDwell_2d& B)
{
  // noch nicht erfolgreich initialisiert?
  if ((min_close_log() == max_close_log()) || (min_open_log() == max_open_log()))
  { 
    sqrt_e(B.sqrt_e()); 
    init(B.name(), 
		 B.min_close_log(), B.max_close_log(),
		 B.min_open_log(), B.max_open_log(),
		 B.bins_per_log());
  }
  // sonst bleiben die grenzen von this gueltig!

  // sind diese auch kompatibel?
  if (bins_per_log() != B.bins_per_log())
  {
    warning((char *)"bins_per_log are different!");
    return *this;
  }

  // beschraenken auf schnittmenge

  // grenzen der schnittmenge
  short log_min_close = std::max(min_close_log(), B.min_close_log());
  short log_max_close = std::min(max_close_log(), B.max_close_log());
  short log_min_open  = std::max(min_open_log(), B.min_open_log());
  short log_max_open  = std::min(max_open_log(), B.max_open_log());

  // offsets
  short close_off = (log_min_close - min_close_log()) * bins_per_log();
  short open_off  = (log_min_open - min_open_log()) * bins_per_log();
  short B_close_off = (log_min_close - B.min_close_log()) * bins_per_log();
  short B_open_off  = (log_min_open - B.min_open_log()) * bins_per_log();

  // summen berechnen 
  for (unsigned short j=0;j < (log_max_close - log_min_close) * bins_per_log();j++)
  {
    for (unsigned short k=0;k < (log_max_open - log_min_open) * bins_per_log();k++)
	{
      CO[j+close_off][k+open_off] += B.CO[j+B_close_off][k+B_open_off]; 
	}  
  }  

  return *this;
}

void 
TDwell_2d::GNUplot(bool dosqrt)		//Tobias 2017 weg: = true
{
double log_laenge1, log_laenge2;

  std::ofstream outCO("closeopen.tmp");

  // falls keine Daten, soll leere fl\E4che gezeigt werden!
  if ((bins_close() == 0) || (bins_open() == 0))
  {
    outCO << "0 0 0" << std::endl << "0 1 0" << std::endl << "0 2 0" << std::endl << std::endl
		  << "1 0 0" << std::endl << "1 1 0" << std::endl << "1 2 0" << std::endl;
  }

  for (unsigned short i=0;i < bins_close();i++)
  {
	log_laenge1 = ((i * (max_close_log() - min_close_log()) / (double)bins_close() + min_close_log()) +
	               ((i+1) * (max_close_log() - min_close_log()) / (double)bins_close() + min_close_log())) / 2;

    for (unsigned short j=0;j < bins_open();j++)
	{
      log_laenge2 = ((j * (max_open_log() - min_open_log()) / (double)bins_open() + min_open_log()) +
	                 ((j+1) * (max_open_log() - min_open_log()) / (double)bins_open() + min_open_log())) / 2;

	  outCO << log_laenge1 << " " 
		    << log_laenge2 << " " 				
			<< CO[i][j] << std::endl;
	}
	outCO << std::endl;					// wenn neue laenge1 leerzeile
  }

  // run gnuplot
  std::ofstream out("gnuplot.tmp");
  out << "set nologscale xyz" << std::endl
      << "set data style lines" << std::endl
      << "set contour base" << std::endl
	  << "set hidden3d" << std::endl
      << "set xlabel \"log(shut time\\)\" 0,0" << std::endl
      << "set ylabel \"log(open time\\)\" 0,0" << std::endl;
  
  if (dosqrt)
  {
	  out << "set zlabel \"sqrt(Number of events)\" 0,0" << std::endl
		<< "set zrange [" << sqrt(min_e(false)) << ":" << sqrt(max_e(false)) << "]" << std::endl;
  } else {
    out << "set zlabel \"Number of events\" 0,0" << std::endl
		<< "set zrange [" << min_e(true) << ":" << max_e(true) << "]" << std::endl;
  }

  out << "set view 60, 30, 1, 1" << std::endl
	  << "set title \"FRONT\"" << std::endl
	  << "splot 'closeopen.tmp' ";
  
  if (dosqrt)
  {
    out << "using ($1):($2):(sqrt($3)) ";
  } else {
    out << "using ($1):($2):($3) ";
  }

  out << "title \"" << name() << "\" " << std::endl
      << "pause -1" << std::endl
      << "set view 60, 210, 1, 1" << std::endl
	  << "set title \"BACK\"" << std::endl
	  << "replot" << std::endl
	  << "pause -1" << std::endl;

  TDwell_lb::GNUplot();
}

double
TDwell_2d::min()
{
  //min ermitteln 
  if (CO.empty()) return 0.0; 
 
  double back = CO[0][0];
  for (unsigned short j=0;j < CO.sizex();j++)
  {
    for (unsigned short k=0;k < CO.sizey();k++)
    {
      if (CO[j][k] < back) back = CO[j][k]; 
    }
  }
  return back;
}
 
double 
TDwell_2d::max()
{
  //min ermitteln
  if (CO.empty()) return 0.0; 

  double back = CO[0][0];
  for (unsigned short j=0;j < CO.sizex();j++)
  {
    for (unsigned short k=0;k < CO.sizey();k++)
    {
      if (CO[j][k] > back) back = CO[j][k]; 
    }
  }
  return back;
}
 
double 
TDwell_2d::events()
{
  //gesamtanzahl ermitteln
  double back = 0;
  for (unsigned short j=0;j < CO.sizex();j++)
  {
    for (unsigned short k=0;k < CO.sizey();k++)
    {
      back += CO[j][k]; 
    }
  }
  return back;
}

TDwell_2d&
TDwell_2d::operator-(TDwell_2d &B)
{				
  // fuer allgemeinere Anwendung muesste Rueckgabetyp TDwell_2d sein
  // dann wuerde Kopie des hier erzeugten Objektes zurueckgegeben
  // so spare ich die kopiererei!
  char   buffer[MaxTextLen];

//  Daten.Dwell_2d_diff.init(buffer);

  // hab ich ueberhaupt zwei datensaetze?
B=Daten.Dwell_2d_B;// Tobias 2018 workaround->not commutative
  if ((CO.empty()) || (B.CO.empty()) )
  {
   if (CO.empty())	
	warning((char *)"matrix CO empty!");
   if (B.CO.empty())	
	warning((char *)"matrix B.CO empty!");	
	
    return Daten.Dwell_2d_diff;
  }
	
  // sind diese auch kompatibel?
  if (bins_per_log() != B.bins_per_log())
  {
    warning((char *)"bins_per_log are different!");
    return Daten.Dwell_2d_diff;
  }
	
  // erstmal skalieren!
  // anschaulicher, wenn differenz im bereich der messdaten!
  // B * factor = A 
  double factor =  events() / (double)B.events();

  // namen eintragen
  sprintf(buffer, "%s - %s", name(), B.name());


  // zrange sichern!
  double min = 0;
  double max = 0;
  if (!(Daten.Dwell_2d_diff.e_auto_range()))
  { 
    min = Daten.Dwell_2d_diff.min_e(false);
    max = Daten.Dwell_2d_diff.max_e(false);
  }
  // vereinigung beider bereiche!
  // diff auf richtige groesse bringen und alle elemente 0!
    		 	

  Daten.Dwell_2d_diff.init(buffer,
                           std::min(min_close_log(), B.min_close_log()),
	                       std::max(max_close_log(), B.max_close_log()),
	                       std::min(min_open_log(), B.min_open_log()),
	                       std::max(max_open_log(), B.max_open_log()),
						   bins_per_log());
  Daten.Dwell_2d_diff.e_range(min, max);
	
  // offsets
  short A_close_off = (Daten.Dwell_2d_diff.min_close_log() - Daten.Dwell_2d_A.min_close_log()) * bins_per_log();
  short A_open_off  = (Daten.Dwell_2d_diff.min_open_log() - Daten.Dwell_2d_A.min_open_log()) * bins_per_log();
  short B_close_off = (Daten.Dwell_2d_diff.min_close_log() - Daten.Dwell_2d_B.min_close_log()) * bins_per_log();
  short B_open_off  = (Daten.Dwell_2d_diff.min_open_log() - Daten.Dwell_2d_B.min_open_log()) * bins_per_log();

  // differenzmatrix berechnen 
  // durch trennung bleibe ich immer im erlaubten bereich!
  // 0 + A

  for (unsigned short j=0;j < CO.sizex();j++)
  {
    for (unsigned short k=0;k < CO.sizey();k++)
    {	
     Daten.Dwell_2d_diff.CO[j-A_close_off][k-A_open_off] += CO[j][k]; 
    }
  }  
  // ... - B * factor
  		
  for (unsigned short j=0;j < B.CO.sizex(); j++)	
  {
   for (unsigned short k=0;k < B.CO.sizey();k++)				
    {
		 Daten.Dwell_2d_diff.CO[j-B_close_off][k-B_open_off] -= B.CO[j][k] * factor; 
    }
  }  
  return Daten.Dwell_2d_diff;
}

double
TDwell_2d::lnlikelihood(TDwell_2d &simul)
{
// bildet loglikelihood aus messwerten und simulierten werten!
// loglikelihood = mess_2d.lnlikelihood(simul_2d)
  
  double back = 0;
  double help; //Tobias
  short xi,yi; //Tobias
  short xmin,ymin,xmax,ymax; //Tobias
  short j,k;
  double Ssimul=0;
  double Smes=0;
  
  
  // sind diese auch kompatibel?
  if (bins_per_log() != simul.bins_per_log())
  {
    warning((char *)"bins_per_log are different!");
    return 0;
  }

  // beschraenken auf schnittmenge
  // wenn mess_2d = 0 -> intervall likelihood = 0
  // wenn simul_2d = 0 -> intervall likelihood = -unendlich -> weglassen

  // grenzen der schnittmenge
  short log_min_close = std::max(min_close_log(), simul.min_close_log());
  short log_max_close = std::min(max_close_log(), simul.max_close_log());
  short log_min_open  = std::max(min_open_log(), simul.min_open_log());
  short log_max_open  = std::min(max_open_log(), simul.max_open_log());
  
  //std::cout<<"log_min_close"<<log_min_close<<std::endl;
  //std::cout<<"log_max_close"<<log_max_close<<std::endl;
  //std::cout<<"log_min_open"<<log_min_open<<std::endl;
  //std::cout<<"log_max_open"<<log_max_open<<std::endl;

  // offsets
  short close_off = (log_min_close - min_close_log()) * bins_per_log();
  short open_off  = (log_min_open - min_open_log()) * bins_per_log();
  short simul_close_off = (log_min_close - simul.min_close_log()) * bins_per_log();
  short simul_open_off  = (log_min_open - simul.min_open_log()) * bins_per_log();

  double faktors = 1.0 /Daten.Fitparameter.multiplikator; //Tobias,  Normierung nicht auf die Anzahl der events, sondern auf die L\E4nge der Zeitreihe


if (Daten.calculate_error == 0)
  {  	
   //Originalcode von OLLI nach Magleby&Weiss
    for ( j=0; j < (log_max_close - log_min_close) * bins_per_log(); j++)
     {
       for ( k=0; k < (log_max_open - log_min_open) * bins_per_log(); k++)
         {
	  if ((simul.CO[j+simul_close_off][k+simul_open_off] != 0) && (CO[j+close_off][k+open_off] != 0))
	    {
             back += log(simul.CO[j+simul_close_off][k+simul_open_off]/  simul.events()  * events() ) * CO[j+close_off][k+open_off];   //Tobias 2018:  without normalization on number of events the fit prefers fast rates, with normalization the term becomes <1, together with 0 as possible results it does not converge / simul.events()
            // std::cout<<j<<":"<<k<<":"<<back<<std::endl;																						  //Magleby & Weiss have simulated a defined number of events equalling the measured data. Otherwise the LLh does not work, because there is no penalty for to much events
	    }																																								  //M & W was modifed by a factor events() incorporating the measured events as normalization. 
	 } 	   
     }
  if (back < 0)																																					//Tobias 2019 required for genetic fit->can not deal with mixed scores
  	back =0;
 	     
  } 

else //einfache quadratische Fehlersumme, 




 if (Daten.calculate_error == 1)  
  {
  for ( j=0; j < (log_max_close - log_min_close) * bins_per_log(); j++)
    {
      for ( k=0; k < (log_max_open - log_min_open) * bins_per_log(); k++)
        {
	
	 //Tobias: Versuch \FCber die quadratische Abweichung
	 help = (simul.CO[j+simul_close_off][k+simul_open_off] * faktors - CO[j+close_off][k+open_off]);
	 
	 Ssimul += simul.CO[j+simul_close_off][k+simul_open_off];
	 Smes   += CO[j+close_off][k+open_off];
	 //std::cout<<"j:"<<j<<" k:"<<k<<" back: "<<back<<" help: "<<help<<" simul: "<<simul.CO[j+simul_close_off][k+simul_open_off]<<" mes: "<<CO[j+close_off][k+open_off]<<" fac: "<<faktors<<std::endl;
	 back += help*help;
        }
    } 

  } 

else //9-fache Fehlersumme, eingef\FChrt von Tobias, um den optischen Eindruck von T\E4lern und H\FCgeln zu ber\FCcksichtigen!

 if (Daten.calculate_error == 2)   
   {
    for ( j=0; j < (log_max_close - log_min_close) * bins_per_log(); j++)
       {
        for ( k=0; k < (log_max_open - log_min_open) * bins_per_log(); k++)
           {
            if (j > 0) xmin = j-1; else xmin = j;
            if (k > 0) ymin = k-1; else ymin = k;
            if (j < ( (log_max_close - log_min_close) * bins_per_log()-1) ) xmax = j+1; else xmax = j;
            if (k < ( (log_max_open - log_min_open) * bins_per_log()-1)   ) ymax = k+1; else ymax = k;
            
            help = 0;
            for (xi = xmin; xi <= xmax; xi++)
              for (yi = ymin; yi <= ymax; yi++)
                {
                 help = help +	(simul.CO[xi+simul_close_off][yi+simul_open_off] * faktors - CO[xi+close_off][yi+open_off]);
                 //std::cout<<"j:"<<j<<"k:"<<k<<":"<<simul.CO[xi+simul_close_off][yi+simul_open_off] * faktors<<":"<<CO[xi+close_off][yi+open_off]<<":"<<help<<std::endl;
                }
            back += help*help;    	
          	
            //std::cout<<"="<<help*help<<" : "<<back<<std::endl;    	
           }      
       }     		
    }

else //25-fache Fehlersumme, eingef\FChrt von Tobias, um den optischen Eindruck von T\E4lern und H\FCgeln zu ber\FCcksichtigen!

 if (Daten.calculate_error == 3)   
   {
    for ( j=0; j < (log_max_close - log_min_close) * bins_per_log(); j++)
       {
        for ( k=0; k < (log_max_open - log_min_open) * bins_per_log(); k++)
           {
            xmin = j-2; if(xmin < 0) xmin=0;
            ymin = k-2; if(ymin < 0) ymin=0;
            xmax = j+2; if(xmax > ( (log_max_close - log_min_close) * bins_per_log()-1 ) ) xmax = (log_max_close - log_min_close) * bins_per_log()-1;
            ymax = k+2; if(ymax > ( (log_max_open - log_min_open) * bins_per_log()-1   ) ) ymax = (log_max_open - log_min_open) * bins_per_log()-1;
                        
            help = 0;
            for (xi = xmin; xi <= xmax; xi++)
              for (yi = ymin; yi <= ymax; yi++)
                {
                 help = help +	(simul.CO[xi+simul_close_off][yi+simul_open_off] * faktors - CO[xi+close_off][yi+open_off]);
                 //std::cout<<"j:"<<j<<"k:"<<k<<":"<<simul.CO[xi+simul_close_off][yi+simul_open_off] * faktors<<":"<<CO[xi+close_off][yi+open_off]<<":"<<help<<std::endl;
                }
            back += help*help;    
            //std::cout<<"="<<help*help<<" : "<<back<<std::endl;    	
           }      
       }     		
    }    

else
std::cout<<"no error sum calculated"<<std::endl;    	     

  
 //std::cout<<"Ergebnis:"<<back<<std::endl;

  return back;
}

double 
TDwell_2d::checkrange(TDwell_2d &simul)
{
// wieviel daten liegen im bereich der simulierten daten?

  double back = 0;

  // sind diese auch kompatibel?
  if (bins_per_log() != simul.bins_per_log())
  {
    warning((char *)"bins_per_log are different!");
    return -1;
  }

  // grenzen der schnittmenge
  short log_min_close = std::max(min_close_log(), simul.min_close_log());
  short log_max_close = std::min(max_close_log(), simul.max_close_log());
  short log_min_open  = std::max(min_open_log(), simul.min_open_log());
  short log_max_open  = std::min(max_open_log(), simul.max_open_log());

  // offsets
  short close_off = (log_min_close - min_close_log()) * bins_per_log();
  short open_off  = (log_min_open - min_open_log()) * bins_per_log();
  //short simul_close_off = (log_min_close - simul.min_close_log()) * bins_per_log();
  //short simul_open_off  = (log_min_open - simul.min_open_log()) * bins_per_log();

  // daten im bereich der simulation aufaddieren 
  for (unsigned short j=0;j < (log_max_close - log_min_close) * bins_per_log();j++)
  {
    for (unsigned short k=0;k < (log_max_open - log_min_open) * bins_per_log();k++)
    {
	  //if ((simul.CO[j+simul_close_off][k+simul_open_off] != 0) &&  //Tobias:macht wegen der nun quadratischen Abweichung keinen Sinn mehr
	  //	  (CO[j+close_off][k+open_off] != 0))
	  {
        back += CO[j+close_off][k+open_off]; 
	  }
    }
  } 

  return back;
}

std::ostream& operator<<(std::ostream& out, TDwell_2d& Dwell)
{
  Dwell.out(out);
  
  out << "# 2d" << std::endl;

  for (unsigned short i=0;i < Dwell.bins_close();i++)
  {
    for (unsigned short j=0;j < Dwell.bins_open();j++)
	{
	  out << Dwell.CO[i][j] << " ";
	}
	out << std::endl;
  }
  return out;
}

std::istream& operator>>(std::istream& in, TDwell_2d& Dwell)
{
char buffer[MaxTextLen];
char c;

  Dwell.in(in);

  in.get(buffer, MaxTextLen, '\n');

  if (!strncmp(buffer, "# 2d",4))			//2019 removed strlwr
  {
    Dwell.init(Dwell.name(),
		       Dwell.min_close_log(), Dwell.max_close_log(),
		       Dwell.min_open_log(), Dwell.max_open_log(), 
			   Dwell.bins_per_log());
	  
	while (in.get(c)) if (c == '\n') break;

    for (unsigned short i=0;i < Dwell.bins_close();i++)
	{
      for (unsigned short j=0;j < Dwell.bins_open();j++)
	  {
		in >>  Dwell.CO[i][j];
	  }
      while (in.get(c)) if (c == '\n') break;
	}
  } else {
    in.clear(std::ios::badbit); // status setzen
  }

  return in;
}

void 
TDwell_2d::load(std::istream& in)
{
// zum laden aus setfile
  TDwell_lb::load(in);

  init(name(),
       min_close_log(), max_close_log(),
	   min_open_log(), max_open_log(), 
	   bins_per_log());

  for (unsigned short i=0;i < bins_close();i++)
  {
    for (unsigned short j=0;j < bins_open();j++)
	{
      in.read((char *)&CO[i][j], sizeof(CO[i][j]));
	}
  }
}

void 
TDwell_2d::store(std::ostream& out)
{
// abspeichern ins jobfile
  TDwell_lb::store(out);


  for (unsigned short i=0;i < bins_close();i++)
  {
    for (unsigned short j=0;j < bins_open();j++)
	{
	  out.write((char *)&CO[i][j], sizeof(CO[i][j]));
	}
  }
}
