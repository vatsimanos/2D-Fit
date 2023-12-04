/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <iomanip>
#include <stdio.h>
// #include <process.h>  //2019

//#include "glib.h"

#include "steuerung.h"
#include "declare.h"
#include "daten.h"
#include "error.h"
#include "round.h"

TChannel_type::TChannel_type(short an_channels, 
                             TRate_matrix aSimulation_Matrix, 
				             TRate_matrix aStartValues,
                             short ai_null, 
				             short ai_channel, 
                             short an_states)
{
  n_Channels = an_channels;
  Simulation_Matrix = aSimulation_Matrix;
  StartValues = aStartValues;
  i_null = ai_null;
  i_channel = ai_channel;
  n_states = an_states;
}

// Information des Setfiles werden in dieser Struktur gespeichert
struct TSetFile
{
  short  version;
  double delta_t;
  char    datafile[fsPathName];
  short  no_channel_types;
  int  no_samples;
  char    noise_file[fsPathName];
  double sigma;
  char    stepresponse[fsPathName];
  std::vector<TChannel_type>  channel_info;
};

//ein Setfile kann Informationen �ber unterschiedliche Kan�le enthalten,
//daher wird eine Liste von PChannel_type-Objeckten angelegt

void Read_Setfile(char* setfile_name, TSetFile &setfile);
void Sort_Matrix(TRate_matrix &Matrix, short n_states);
void Scale_Matrix(TRate_matrix &Matrix, short n_states);

// funktion zum steuern der erstellung der dwelltime-dateien mit vogegebenen werten (aus simulations-set-datei)
void  
Steuerung(char* setfilename, short histopts, double koop, bool test)
{
// TDaten Daten; // ist schon global
char*    p;
TSetFile Setfile;
char     speichername[MaxTextLen];
short   i,j;

// taus
short   levelno;
double  dtus;     // delta_t in Mircosec

  char deftxt[] = "default.txt";
  std::ofstream Textaus(deftxt);
  if (!Textaus)
  {
    char infotext[MaxTextLen];
    sprintf(infotext,"Can't open file: %s",deftxt);

    warning(infotext);           // fehlermeldung wohin auch immer
    return;                      // keine daten lesen!
  }

  Textaus << "Setfile: " << setfilename << std::endl;
  Read_Setfile(setfilename, Setfile);
  // Ausgabe zu Testzwecken:
  Textaus << "version: " << Setfile.version << std::endl;
  Textaus << "datafile: " <<  Setfile.datafile << std::endl;
  Textaus << "delta_t: " <<  Setfile.delta_t << std::endl;
  if (Setfile.channel_info.size() > 0)
  {
    Textaus << "i_null: " << Setfile.channel_info[0].i_null << std::endl;
    Textaus << "i_channel: " << Setfile.channel_info[0].i_channel << std::endl;
    Textaus << "n_channels: " << Setfile.channel_info[0].n_Channels << std::endl;
    Textaus << "n_states: " << Setfile.channel_info[0].n_states << std::endl;
    Textaus << "sigma: " << Setfile.sigma << std::endl;
    Textaus << "rates: " << std::endl;
    // Matrix aus Setfile so umsortieren, da� open Zust�nde immer zuerst kommen. Dies ist f�r die Fit-Routine erforderlich
    Sort_Matrix(Setfile.channel_info[0].Simulation_Matrix,
                Setfile.channel_info[0].n_states);
    // auch fuer startvalues
    Sort_Matrix(Setfile.channel_info[0].StartValues,
                Setfile.channel_info[0].n_states);
    Textaus << "changed order of matrix elements" << std::endl;
    // Matrix Scalieren
    Scale_Matrix(Setfile.channel_info[0].Simulation_Matrix,
                 Setfile.channel_info[0].n_states);
    // auch fuer startvalues
    Scale_Matrix(Setfile.channel_info[0].StartValues,
                 Setfile.channel_info[0].n_states);
    Textaus << "matrix scaled" << std::endl;
    for (i=0;i <= n_states_max-1;i++)
	{
      for (j=0;j <= n_states_max-1;j++)
	  {
        Textaus << std::setw(10) << Setfile.channel_info[0].Simulation_Matrix.r[i][j] << "   "; //:6:3
	  }
      if (Setfile.channel_info[0].Simulation_Matrix.open[i])
	  {
		Textaus << " open" << std::endl;
	  } else {
		Textaus << " close" << std::endl;
	  }
	}

    Textaus << "start values" << std::endl;
    // auch startwerte ausgeben
    for (i=0;i <= n_states_max-1;i++)
	{
      for (j=0;j <= n_states_max-1;j++)
	  {
        Textaus << std::setw(10) << Setfile.channel_info[0].StartValues.r[i][j] << "   "; //:6:3
	  }
      if (Setfile.channel_info[0].StartValues.open[i])
	  {
        Textaus << " open" << std::endl;
	  } else {
		Textaus << " close" << std::endl;
	  }
    }
    // Nur Test des Setfiles auf richtige Syntax?
    if (!test)
	{
      // einlesen der zeitreihe
      Daten.get_data_from_file(Setfile.datafile);
      if (koop != 0)
	  {
        Textaus << "model with cooperativity" << std::endl;
        Daten.Kooperativity = koop;
        Textaus << "cooperativity: " << koop << std::endl;
        Daten.modell = M_Koop;
	  } else {
        Daten.modell = M_5State;
	  }
      Textaus << "timeseries loaded" << std::endl;

	  // !!!Abtastinterval in us!!!
      Daten.Settings.Samplininterval = Setfile.delta_t * 1000000;                      
      // nullniveau
	  Daten.a_fit.i_null = Setfile.channel_info[0].i_null;         
	  // amplitude des stromes = differenz
      Daten.a_fit.i_channel = Setfile.channel_info[0].i_channel;   
      // zahl der Kan�le
      Daten.a_fit.n_channels = Setfile.channel_info[0].n_Channels; 
      // bei benutzung eines filters muss hier ein faktor multipliziert werden

	  // standardabweichung
      Daten.a_fit.sigma = Setfile.sigma;
	  // aufruf des filters
      Daten.hinkley(_HOHD);                     

      Textaus << "jumps detected, ready for fiting" << std::endl;

      // schon mal .jbf-file speichern, da program w�herend des Fit abst�rzen k�nnte
      Daten.store((char *)"default.jbf");

      // abspeichern der dwell-time-files, nicht n�tig, da sofort gefittet wird
      // Daten^.SaveDwellTimeHisto(FileName2); //PASCAL kommentar
      if (histopts != 0)
	  {
	    Daten.SDwellTime.Samples = histopts;
	  } else {
        Daten.SDwellTime.Samples = 100; // anzahl der datenpunkte im histogramm auf max. 250 setzen
	  }
      Textaus << "points in histogram: " << Daten.SDwellTime.Samples << std::endl;
      Daten.SDwellTime.openclosed = true;  // Spr�nge up und down addiert betrachten
      Daten.SDwellTime.ExpAvarage = true;  // Exp. Mittelung, nicht lineare, ausw�hlen
      Daten.ExpEinteil();                  // berechnung des growfaktors
      Daten.mittel();                      // Mittelung durchf�hren
      Textaus << "exp. avarage calculated" << std::endl;

      // uebergabe der ratenkonstanten, in [ms] d.h. fuer [s]: *(1/1000)
      // nur f�r Rikards fit Routine
  
      // PASCAL kommentar
      // Daten^.Rates.CO:=?/1000;
      // Daten^.Rates.OC:=?/1000;
      // Daten^.Rates.OG:=?/1000;
      // Daten^.Rates.GO:=?/1000;
      // Daten^.modell := M_SCO;
  
      // Der Fit wird hier mit den 'wahren' Raten initialisiert, es k�nnten auch die Startwerte sein
      // Daten^.Rate_Matrix:=PChannel_type(SetFile.channel_info^.At(0))^.Simulation_Matrix;
      // Fuer Afshin doch mit startwerten!
      // wenn setfile version<2 werden startwerte aus simulations_matrix generiert s.u.
      Daten.Rate_Matrix = Setfile.channel_info[0].StartValues;

      //Daten.tfit();
      Textaus << "fit done" << std::endl;
      // .jbf-file Speichern
      Daten.store((char *)"default.jbf");
      // Ergebnis abspeichern
  
      // Kooperativit�tsfaktor abspeichern
      if (koop != 0)
	  {
        strcpy(speichername, setfilename); 
        p = strrchr(speichername,'.');
        if (p != NULL) *p = 0; // endung abschneiden 
	    strcat (speichername, ".k");
        std::ofstream Raten(speichername,std::ios::out|std::ios::app);
        if (!Raten)
		{
          Textaus << "file " << speichername << " not found, creating new file" << std::endl;
          Raten.close();
          Raten.open(speichername,std::ios::out);
		  if (!Raten)
		  {
            Textaus << "Can't create file: " << speichername << std::endl;
            return;                      // keine daten lesen!
		  }
		}
		Raten << Daten.Kooperativity << std::endl; //:5:6
        Raten.close();
	  }
      // Fehler abspeichern
      strcpy(speichername, setfilename); 
      p = strrchr(speichername,'.');
      if (p != NULL) *p = 0; // endung abschneiden 
      strcat (speichername, ".err");
      std::ofstream Raten(speichername,std::ios::out|std::ios::app);
      if (!Raten)
	  {
        Textaus << "file " << speichername << " not found, creating new file" << std::endl;
        Raten.close();
        Raten.open(speichername,std::ios::out);
	    if (!Raten)
		{
          Textaus << "Can't create file: " << speichername << std::endl;
          return;                      // keine daten lesen!
		}
	  }
      Raten << Daten.DwellLevel[0][0][openclosing].ErrorFiveState << "   " 
	        << Daten.DwellLevel[0][0][openclosing].iter << std::endl;
      Raten.close();

      // raten speichern
      for (i=0;i <= n_states_max-1;i++)
	  {
        for (j=0;j <= n_states_max-1;j++)
		{
          if (Daten.Rate_Matrix.r[i][j] != 0)
		  {
            strcpy(speichername, setfilename); 
            p = strrchr(speichername,'.');
            if (p != NULL) *p = 0; // endung abschneiden 
            if (Daten.Rate_Matrix.open[i])
			{
			  strcat(speichername, ".o");
			} else {
			  strcat(speichername, ".c");
			}
            char buffer[5]; // eigentlich nur zwei ziffern!
            sprintf(buffer,"%i%i", i, j);
			strcat(speichername, buffer);

            std::ofstream Raten(speichername,std::ios::out|std::ios::app);
            if (!Raten)
			{
              Textaus << "file " << speichername << " not found, creating new file" << std::endl;
              Raten.close();
              Raten.open(speichername,std::ios::out);
	          if (!Raten)
			  {
                Textaus << "Can't create file: " << speichername << std::endl;
                return;                      // keine daten lesen!
			  }
			}
            Raten << Daten.Rate_Matrix.r[i][j] << std::endl; //:5:6)
            Raten.close();
		  }
		}
	  }
      // Zeitkonstanten speichern! fuer Afshin
      // zu jedem Level Zeitkonstanten in jeweils ein File
      // eigentlich reichen ihm Level O und das oberste

      dtus = Setfile.delta_t * 1000000; // delta_t in Mircosec s.a. ExpFitDlg

      levelno = 0;
      while (levelno <= Setfile.channel_info[0].n_Channels)
	  {
        strcpy(speichername, setfilename); 
        p = strrchr(speichername,'.');
        if (p != NULL) *p = 0; // endung abschneiden 
        strcat(speichername, ".l");
        char buffer[7]; // eigentlich nur zwei ziffern!
        sprintf(buffer,"%i", levelno);
	    strcat(speichername, buffer);

        std::ofstream Raten(speichername,std::ios::out|std::ios::app);
        if (!Raten)
		{
          Textaus << "file " << speichername << " not found, creating new file" << std::endl;
          Raten.close();
          Raten.open(speichername,std::ios::out);
          if (!Raten)
		  {
            Textaus << "Can't create file: " << speichername << std::endl;
            return;                      // keine daten lesen!
		  }
          Raten << "taus in mircosec, level: " << levelno 
		        << ". delta_t: " << dtus << " microsec." << std::endl; //:6:3
		}
        // erstes tau kennzeichnen, falls Zeile zu kurz
        Raten << "-> ";
        // Zeile aller konstanten ausgeben
        for (i=1;i <= NMaxTimeConstants;i++)
		{
          // zeitkonstanten ermitteln und ausgeben
          Raten << Daten.DwellLevel[levelno][0][openclosing].tau[i] * dtus //:6:3 // ExpFitDlg
                << "\t"; // tabulator
		}
        Raten << std::endl;

        Raten.close();
        levelno++;
	  } // bis hier fuer afshin
	}
  }
}

void  
Steuerung(char* setfilename, char* dwell, short histopts, double minevents)
{
}

void  
Steuerung(char* setfilename, short histopts)
{
}

void 
Copy_Matrix(TRate_matrix &DestMatrix, TRate_matrix SourceMatrix, short n_states)
{
short  i,j;

  for (i=0;i <= n_states-1;i++)
  {
	for (j=0;j <= n_states-1;j++)
	{
      DestMatrix.r[i][j] = SourceMatrix.r[i][j];
	}
  }
  for (i=0;i <= n_states-1;i++)
  {
    DestMatrix.open[i] = SourceMatrix.open[i];
  }
}

void 
Sort_Matrix(TRate_matrix &Matrix, short n_states)
{
TRate_matrix r;
short       i,j,k;

  k = 0;
  // Zeilen sortieren
  for (i=0;i <= n_states-1;i++)
  {
    if (Matrix.open[i])
	{
      for (j=0;j <= n_states-1;j++)
	  {
        r.r[k][j] = Matrix.r[i][j];
        r.open[k] = true;
	  }
      k++;
	}
  }
  for (i=0;i <= n_states-1;i++)
  {
    if (!Matrix.open[i])
	{
      for (j=0;j <= n_states-1;j++)
	  {
        r.r[k][j] = Matrix.r[i][j];
        r.open[k] = false;
	  }
      k++;
	}
  }
  k = 0;
  // Spalten sortieren
  for (i=0;i <= n_states-1;i++)
  {
    if (Matrix.open[i])
	{
      for (j=0;j <= n_states-1;j++)
	  {
        Matrix.r[j][k] = r.r[j][i];
	  }
        k++;
	}
  }
  for (i=0;i <= n_states-1;i++)
  {
    if (!Matrix.open[i])
	{
      for (j=0;j <= n_states-1;j++)
	  {
        Matrix.r[j][k] = r.r[j][i];
	  }  
      k++;
	}
  }
  // open-close Flags setzen
  for (i=0;i <= n_states-1;i++)
  {
    Matrix.open[i] = r.open[i];
  }
}

void 
Scale_Matrix(TRate_matrix &Matrix, short n_states)
{
short  i,j;

  for (i=0;i <= n_states-1;i++)
  {
    for (j=0;j <= n_states-1;j++)
	{
      Matrix.r[i][j] /= 1000;
	}
  }
}

char* 
Next_line(std::ifstream &fhsetfile, char* line, short &line_counter, bool killcomments = true)
{
char     dum[fsPathName];
bool gefunden;
short   i;

  gefunden = false;
  do
  {
    fhsetfile.getline(dum, MaxTextLen);
    line_counter++;

	if (killcomments)
	{
      char* x = strchr(dum, '#');
      if (x != NULL) *x = 0; // kommentare abschneiden 
	}

    for (i=short('0');i <= short('9');i++)
	{
      if (strchr(dum, char(i)) != NULL) gefunden = true;
	}
    for (i=short('a');i <= short('z');i++)
	{
      if (strchr(dum, char(i)) != NULL) gefunden = true;
	}
    for (i=short('A');i <= short('Z');i++)
	{
      if (strchr(dum, char(i)) != NULL) gefunden = true;
	}
  } while (!gefunden);

  strcpy(line, dum);

  return line;
}

char* 
Get_text(char* text)
{
short  i;
char    b;
char*   x;

  x = text + strlen(text); 

  i = 0;
  do
  {
    b = text[i];
    i++;
    if (((short(b) > short('a')) && (short(b) < short('z'))) || 
		((short(b) > short('A')) && (short(b) < short('Z'))))
	{
      x = text + i-1;
      break;
	}
  } while (short(b) != 0);
  return x;
}

char* 
Get_val(char* line, double &x)
{
short  code,i;
char    zahl[256];

  code = 1;
  i = 0;
  x = 0;
  do
  {
    if (((char(line[i]) >= char('0')) && (char(line[i]) <= char('9'))) || (char(line[i]) == '.'))
	{
      line += i;
      strcpy(zahl, line);
      break;
	}
    i++;
  } while (char(line[i-1]) != 0);
  i=0;
  do
  {
    if (((char(line[i]) < char('0')) || (char(line[i]) > char('9'))) && (char(line[i]) != char('.')))
	{
      line += i;
      i++;
	  break;
	}
    i++;
  } while (char(line[i-1]) != 0);
  i--;
  zahl[i] = 0;
  x = atof(zahl);
  return line;
}

bool 
Read_Matrix(std::ifstream &fhsetfile, TRate_matrix &m, short &line_counter, short &dim)
{
char*    pline;
char     line[256];
short   i,j;
double  dim1,dim2,r;

  pline = Next_line(fhsetfile, line, line_counter);
  pline = Get_val(pline, dim1);
  pline = Get_val(pline, dim2);
  if ((dim1 != 0) && (dim1 != 0))
  {
    dim = round(dim1);
    for (i=1;i <= round(dim1);i++)
	{
      pline = Next_line(fhsetfile, line, line_counter);
      for (j=1;j <= round(dim2);j++)
	  {
        pline = Get_val(pline, r);
        if ((i > n_states_max) || (j > n_states_max)) break; // Vorsichtsma�nahme gegen Zugriffe auserhalb des Arrays
        m.r[i-1][j-1] = r;
	  }
	}
	return true;
  } else {
    char infotext[MaxTextLen];
    sprintf(infotext, "error reading matrix dimensions in line %i", line_counter);
    return false;
  }
}

bool 
Read_Row(std::ifstream &fhsetfile, double row[], short &line_counter)
{
char*   pline;
char    line[256];
short  i;
double j;
double dim1,dim2;

  pline = Next_line(fhsetfile, line, line_counter);
  pline = Get_val(pline, dim1);
  pline = Get_val(pline, dim2);
  if ((dim1 == 1) && (dim2 >=2))
  {
    pline = Next_line(fhsetfile, line, line_counter);
    for (i=1;i <= round(dim2);i++)
	{
      pline = Get_val(pline, j);
      row[i-1] = j;
	}
	return true;
  } else {
    char infotext[MaxTextLen];
    sprintf(infotext, "error reading matrix dimensions in line %i", line_counter);
    return false;
  }
}

void 
Read_Setfile(char* setfile_name, TSetFile &setfile)
{
char    filename[MaxTextLen];
char    line[MaxTextLen];
char    infotext[MaxTextLen+30];
char*   p;
char*   pline;
char*   token;
char    seps[] = " \t=";  // space, tabs und '=' als sep in setfile
short  i,j,line_counter;
short  error_count;
short  no_channel_type_1 = 0;
short  no_channels;
short  channel, null;
TRate_matrix  k;
TRate_matrix  k_start;
TRate_matrix  bla;
double row[n_states_max];
double assignment[n_states_max]; 
double i_current[n_states_max]; // max 5 Zust�nde
short  n_states, dum;
  strcpy(filename, setfile_name); 
  p = strrchr(filename,'.');
  if (p != NULL) *p = 0; // endung abschneiden 
  strcat (filename, ".set");
  std::ifstream fhsetfile(filename,std::ios::in);
  if (!fhsetfile)
  {
	sprintf(infotext, "Can't open setfile: %s", filename);
    warning (infotext);
    return;                      // keine daten lesen!
  }  
  // Versions einlesen
  fhsetfile.getline(line, MaxTextLen);
  token = strstr(line, "# sv");			//2019 removed strlwr
  if (token != NULL)
  {
	setfile.version = atoi(token+5);
  } else {
    setfile.version = 1;
  }

  // header einlesen 
  line_counter = 1;
  for (i=1;i <= 7;i++)
  {
    error_count = 0;

    pline = Next_line(fhsetfile, line, line_counter, false);

    // Establish string and get the first token:
    token = strtok(line, seps);

    if (!strcmp(token, "delta_t"))			//2019 removed strlwr
	{
	  token = strtok(NULL, seps);
      setfile.delta_t = atof(token);
	  error_count++;
	}
    if (!strcmp(token, "datafile"))		//2019 removed strlwr
	{
	  token = strtok(NULL, seps);
      strcpy(setfile.datafile, token);
      error_count++;
	}
    if (!strcmp(token, "no_channel_types"))	//2019 removed strlwr
	{
	  token = strtok(NULL, seps);
      setfile.no_channel_types = atoi(token);
      error_count++;
	}
    if (!strcmp(token, "channel_no"))		//2019 removed strlwr
	{
	  token = strtok(NULL, seps);
      no_channel_type_1 = atoi(token);
      error_count++;
	}
    if (!strcmp(token, "no_samples"))		//2019 removed strlwr
	{
	  token = strtok(NULL, seps);
      setfile.no_samples = atoi(token);
      error_count++;
	}
    if (!strcmp(token, "noise_file"))		//2019 removed strlwr
	{
	  token = strtok(NULL, seps);
      strcpy(setfile.noise_file, token);
      error_count++;
	}
    if ((!strcmp(token, "sigma_noise")) ||	//2019 removed strlwr
        (!strcmp(token, "sigma")))
	{
	  token = strtok(NULL, seps);
      setfile.sigma = atof(token);
      error_count++;
	}
    if (!strcmp(token, "stepresponse"))		//2019 removed strlwr
	{
	  token = strtok(NULL, seps);
      strcpy(setfile.stepresponse, token);
      error_count++;
	}
	
    if (error_count != 1)
	{
	  sprintf(infotext, "error while reading header of setfile in line %i", line_counter);
      warning(infotext);
      return;
	}
  }

  if (setfile.version == 1) setfile.no_channel_types = 1;
  for (i=1;i <= setfile.no_channel_types;i++)
  {
    if (setfile.version >= 2)
	{
      pline = Next_line(fhsetfile, line, line_counter, false);
      pline = strchr(line, '=');
      if (pline == NULL)
	  {
	    sprintf(infotext, "no '=' character in line %i", line_counter);
        warning(infotext);
        return;
	  }
      no_channels = atoi(pline+1);
	} else {
      no_channels = no_channel_type_1;
	}

    // setfile version 1 hat keine kommentare vor ueberschrift!)
    Next_line(fhsetfile, line, line_counter, false);
    if (!Read_Matrix(fhsetfile, k, line_counter, n_states))
	{
	  sprintf(infotext, "error reading simualtion matrix in line %i", line_counter);
      warning(infotext);
      return;
	}
    if (setfile.version >= 2)
	{
      if (!Read_Matrix(fhsetfile, k_start, line_counter, dum))
	  {
        sprintf(infotext, "error reading start values for fitting in line %i", line_counter);
        warning(infotext);
        return;
	  }
    } else {
      // setfile version < 2 -> start values automatisch aus simulation matrix
      Copy_Matrix(k_start, k, n_states);
	}

    Next_line(fhsetfile, line, line_counter, false);
    if (!Read_Matrix(fhsetfile, bla, line_counter, dum))
	{ 
	  // indices of the parameters in k to be fitted
      sprintf(infotext, "error reading indices of the parameters in K to be fitted in line %i", line_counter);
      warning(infotext);
      return;
  	}

    Next_line(fhsetfile, line, line_counter, false);
    if (!Read_Row(fhsetfile, row, line_counter))
	{
      sprintf(infotext, "error reading initial state distibution in line %i", line_counter);
      warning(infotext);
      return;
	}

    Next_line(fhsetfile, line, line_counter, false);
    if (!Read_Row(fhsetfile, assignment, line_counter))
	{
      sprintf(infotext, "error reading assignment states <-> output symbols in line %i", line_counter);
      warning(infotext);
      return;
	}

    Next_line(fhsetfile, line, line_counter, false);
    if (!Read_Row(fhsetfile, i_current, line_counter))
	{
      sprintf(infotext, "error reading output levels in line %i", line_counter);
      warning(infotext);
      return;
	}

    // Auswertung des Setfiles, nicht unbedingt allgemeing�ltig
    // Nullstrom und Kanalstrom setzen
    if (i_current[0] > i_current[1])
	{
      null = round(i_current[1]);
      channel = round(i_current[0] - i_current[1]);
	} else {
      null = round(i_current[0]);
      channel = round(i_current[1] - i_current[0]);
    }
    // Open und Closed Zust�nde finden
    for (j=0;j <= n_states-1;j++)
	{
      if (round(assignment[j]) == 1)
	  {
        if (i_current[0] < i_current[1]) 
		{
		  k.open[j] = false;
		} else {
		  k.open[j] = true;
		}
	  } else {
        if (i_current[0] < i_current[1])
		{
		  k.open[j] = true;
		} else {
		  k.open[j] = false;
		}
	  }
	}
    // Open und Closed Zust�nde auch fuer start values finden
    for (j=0;j <= n_states-1;j++)
	{
      if (round(assignment[j]) == 1)
	  {
        if (i_current[0] < i_current[1])
		{
		  k_start.open[j] = false;
		} else {
		  k_start.open[j] = true;
		}
	  } else {
        if (i_current[0] < i_current[1])
		{
		  k_start.open[j] = true;
		} else {
		  k_start.open[j] = false;
		}
	  }
	}

    TChannel_type ct(no_channels, k, k_start, null, channel, n_states);
    setfile.channel_info.push_back(ct);
  }
}
