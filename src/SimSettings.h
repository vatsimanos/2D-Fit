/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  2000 Erweiterung von Tobias Huth
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/
 struct TImportSettings			//Tobias 2013
   {
    unsigned int ts_start;
    unsigned int EPB_start;
    unsigned int EPB_length;
    unsigned int ts_end;
   };
 
 struct TSimSettings
  {
   char loadname[fsPathName];
   bool actionpotentials;
   int   repeats;
   int   markercounts;
   unsigned int  sampling_frequency; 
   unsigned int  filter_frequency;
   unsigned int  noise_length;
   // char    noisefile[255];  TObias 2022 noisefile is now loaded via settings.ini
   bool resample;
   bool dosbox;		//Toby 2013
   bool show_2D_Diff; 	//Toby 2013
   char check[fsPathName]; 				//Tobias 2013 Levelfit checksum
   bool write_result_file;
   char result_file_name[fsPathName];
  };
  
 struct TFitparameter
  {
    double  spannweite;     	    // Spannweite des Simplex (f�r die Initialisierung der Startwerte)
    double  ftol;                  //Genauigkeit->Abbruchbedingung f�r den Fit
    int   multiplikator;         // Anzahl der wiederholungen der Zeitreihe / Fit 
    int   maxmultiplikator;      // die Genauikeit (Anzahl der Simulationen) werden von min nach max gesteigert! 
    int   minmultiplikator;
    int   maxiterations;         //Maximum allowed number of function evaluations.  
    short   log_min_close;
    short   log_max_close;
    short   log_min_open;
    short   log_max_open;
    short   bins_per_log;
  };
  
 struct TGeneticfitparameter //Struktur f�r die Settings des genetic fits
  {
   bool levelfit;			//Tobias 2012 Levelfit
   bool proportionalfit;		//Tobias 2013 Levelfit
   double stromvektor[11];		//Tobias 2012 Levelfit
   short  state_conductance[11];	//Tobias 2013 Levelfit
   double levelfit_min[3];		//Tobias 2013 Levelfit
   double levelfit_max[3];		//Tobias 2013 Levelfit
   double level1;			//Tobias 2013 Levelfit
   double level2;			//Tobias 2013 Levelfit
   bool automated_loading;		//Tobias 2013 Levelfit
   char files_directory[fsPathName];		//Tobias 2013 Levelfit
   char file_time_series[fsPathName];		//Tobias 2013 Levelfit
   char file_noise_series[fsPathName];		//Tobias 2013 Levelfit
   short noise_generation;								//2023
   bool file_noise_series_loaded;					//Tobias 2022
   char file_setfile[fsPathName];		//Tobias 2013 Levelfit
   bool mpi_support;			//Tobias 2013 Levelfit
   int nr_calc_mean_result;			//Tobias 2013 
   bool elitism;				//Tobias 2014
   	
   bool logscale;	
   short fit_methode;
   int popsize;                           
   int ngen;
   double pmut;
   double pcross;
   short bits;
   int nConvergence;
   double startwertmin[101];
   double startwertmax[101];
   short methode;
   double likelihood_2d;			//2D_error global f�r �bergabe
   double histogram_diff;			//amplitude histogram error f�r �bergabe
   int init_random;
  }; 	 
  

  
 struct  TMultiFit       //struct f�r multifit-fit f�r mehrere Setfiles hintereiander
  {
    char  setfile[255];   //setfiles
    char  logfile[255];  //logfiles
  };  	 
  
