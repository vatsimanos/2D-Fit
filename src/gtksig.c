/*******************************************************************
 *  Patch Clamp 
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <iostream>
#include <stdio.h>
#include <errno.h>
#include <time.h>

#include "gtkapp.h"
#include "declare.h"
#include "daten.h"
#include "gtkfunc.h"
#include "round.h"
#include "steuerung.h"
#include "error.h"
#include "sim_fit.h" 


int
main(int argc, char *argv[])
{
char     setfilename[MaxTextLen];
short   histopts = 0;
double  minevents = 100000;
double  koop = 0;
char*    dwell = NULL;
bool test = false;
short   i=0;

	

//BEGIN KARSTENASCI
Daten.AD_fromwert=0;
Daten.AD_towert=4095;
//END KARSTENASCI

   
  // jetzt settings.ini laden
   std::ifstream in;				
   in.open("settings.ini");
  in >> Daten.Settings;
  // ist richtige datei vorhanden?
  if (!in) std::cout<<"no settings.ini"<<std::endl;
  // egal was drin steht zu beginn immer 0!
  in.close();
  Daten.Settings.SamplesToSkip = 0;
  Daten.Settings.SamplesToProcess = 0;


 // Tobias 2013 init MPI=================================================================================							
int myid;
myid = 0;	//in case MPI not supported
int size;
char processor_name[255];
int namelen;
time_t nu;


if (Daten.Geneticfitparameter.mpi_support)
   {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); 
    MPI_Comm_size(MPI_COMM_WORLD,&size);  
    MPI_Get_processor_name( processor_name, &namelen ); 
   }
  	 	       							
if (Daten.Geneticfitparameter.init_random == 0)		//randomized timer 
	{
	 int seed =	time(&nu)*myid; // Efthymios 2023 (previously not properly initialized when using parallel computing)// old: int seed = time(&nu)+ myid + myid << 24;
   Daten.random_engine.seed(seed);											//random initializing of generator
   GARandomSeed(Daten.uni_01(Daten.random_engine));												//same for GA
  }
else
	{
	 Daten.random_engine.seed(Daten.Geneticfitparameter.init_random);	//determined random sequence
   GARandomSeed(Daten.Geneticfitparameter.init_random);	
  }	   

		  								
// Tobias 2013 automatic start=================================================================================


if (Daten.Geneticfitparameter.automated_loading)
   {
    Daten.get_data_from_file (Daten.Geneticfitparameter.file_time_series);
    strcpy(Daten.SimSettings.loadname,Daten.Geneticfitparameter.file_setfile);
    Daten.reset_level();  
   }
   
// Tobias 2013 init MPI=================================================================================						
if (Daten.Geneticfitparameter.mpi_support)
   {	
    multiprocessor_fit();
    return 0;
   }						
// Tobias 2013 init MPI=================================================================================    
    
    
      Steuerung(setfilename, dwell, histopts, minevents);
  return 0;
}



void
exec_capscreen(char* cmdstring, char* optstring)
{
}

