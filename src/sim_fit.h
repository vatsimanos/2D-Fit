/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  2000 Tobias Huth 
 *  based on uwe kirst's daytest  
 *  CAU KIEL 
 * Genetic Fit: Copyright 1995-1996 Massachusetts Institute of Technology
 *******************************************************************/


#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include "ga/ga.h"
#include "gtkfunc.h"
#include "setfile.h"



//f�r den simplex
#define TINY 1.0e-10  //A small number.
#define GET_PSUM \
		  for (j=1;j<=ndim;j++) {\
		  for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
		  psum[j]=sum;}
#define SWAPtoby(a,b) {swap=(a);(a)=(b);(b)=swap;}

//von olli kopiert
inline void 
store_path(char* FileName)
{
  // nur dir, ohne filename sichern
  strcpy(Daten.cdir, FileName);
  char* p = strrchr(Daten.cdir,'\\');
  if (p != NULL)
  {
    p[1] = 0; 
  } else {
	Daten.cdir[0]=0; 
  }
}
 
// enth�lt die Parameter f�r Simulation & Fit
//void on_CM_SimSettings_activate(GtkMenuItem* menuitem, gpointer user_data);
void load_setfile();
//void load_setfile_OK(GtkWidget* widget);
//Funktion zum Darstellen von Dependency-Plots (Magleby & Song)
//void
//on_CM_2D_Dependency_activate(GtkMenuItem* menuitem, gpointer user_data);

//void
//on_CM_MaximumLikelihood_activate(GtkMenuItem* menuitem, gpointer user_data);
//BEGIN KARSTENFIT
//void
//on_end_MaximumLikelihood(GtkMenuItem* menuitem, gpointer user_data);
//void
//on_start_MaximumLikelihood(GtkMenuItem* menuitem, gpointer user_data);
//void
//on_stop_MaximumLikelihood(GtkMenuItem* menuitem, gpointer user_data);    
void
Save_Timeseries_and_Amplitudenhistogramm();       
//END KARSTENFIT
void destroy_fit_window();
//void destroy_simulat_window(GtkMenuItem* menuitem, gpointer user_data);      
//void destroy_settings_window(GtkMenuItem* menuitem, gpointer user_data);   
// �ffnet Fenster f�r Simulationsparameter
//void
//on_CM_Simulation_activate(GtkMenuItem* menuitem, gpointer user_data);
// akiviert den 2DFit
//void
//on_CM_2DFit_activate(GtkMenuItem* menuitem, gpointer user_data);
// Prozedur f�r den Fit
//void
//on_CM_2DFit(GtkMenuItem* menuitem, gpointer user_data);
//Aufruf des Genetic-Fits
//void 
//on_CM_GeneticFit(GtkMenuItem* menuitem, gpointer user_data);
//void on_load_setfile_activate(GtkMenuItem* menuitem, gpointer user_data);
//Tobias 2013
//void on_load_ATFfile_activate(GtkMenuItem* menuitem, gpointer user_data);
//void load_ATFfile_OK(GtkWidget* widget);
//void destroy_ATF_import_window(GtkMenuItem* menuitem, gpointer user_data);

//void on_load_noisefile_activate(GtkMenuItem* menuitem, gpointer user_data);
//void load_noisefile_OK(GtkWidget* widget);

//void on_load_stepfile_activate(GtkMenuItem* menuitem, gpointer user_data);
//void load_stepfile_OK(GtkWidget* widget);

// Funktion f�r Lotgbuch der Spr�nge nach Abfolge sortiert
//void
//on_CM_HinkleyLog_activate(GtkMenuItem* menuitem, gpointer user_data);


// Simuliert Zeitreihen
//void simulate_timeseries(GtkMenuItem* menuitem, gpointer user_data);

//leistet die Vorarbeit zum jeweiligen (simplex, genetic)Fit
void 
ZweiDFit();


// simuliert mit nr_const Ratenkonstanten, bis >=events simuliert sind (die Ratenkonstanten sind in *rconst enthalten)
double simulation(int multiplikator, int nr_const, double *rconst, double rtol);
		      
                                  
void create_fit_window();       

bool Simulat_window();

void out_time();

unsigned int GeneticFit_SimpleGA();
float objective(GAGenome &);
void PopulationInitializer(GAPopulation &);
void PopulationEvaluator(GAPopulation &);

void
multiprocessor_fit();		//Tobias 2013

//void
//on_CM_log_dwell_add_marked_activate(GtkMenuItem* menuitem, gpointer user_data);

short
load_setfile(char *Filename);

void
init_2Dfit();

void
out_2D_settings();

void
generation_out_2Dfit(GASimpleGA &,GABin2DecGenome &);

void
summary_out_2Dfit(GASimpleGA &, GABin2DecGenome &);

void
multiprocessor_master();

void
multiprocessor_slave();

void
import_ATF();

void
multiprocessor_master_save_timeseries(); // Tobias 2020

void
multiprocessor_slave_save_timeseries();		// Tobias 2020


