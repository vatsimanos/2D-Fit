/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  2000 Erweiterung von Tobias Huth
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

// initilaisiert die Parameter

 markerstart = false;
 lauf = Daten.Marker.begin();
 markernumber = 1;
 iwvisible = true;
 nrselected = 0;
 paintinaction = false;
 TwoD_fitinaction = false;
 loeschanzahl = 0;
 loeschvalue  = 2048;
 geneticfit = true;
 calculate_error = 2;  //0=-lnlikelihood, 1=1bin, 2=filter 9bin QA, 3=filter 25bin QA
 firstsimul = true;
 simulate_memory = false;
 errormessage[0] = 0; //sprintf(errormessage, "");
 Daten.logx = log(5000)/50; 
 cut_marked_save =  0;
 cut_start = 0;
    
    //Tobias 2013
    ImportSettings.ts_start = 1;
    ImportSettings.EPB_start = 20656;
    ImportSettings.EPB_length = 20000;
    ImportSettings.ts_end = 42000;
 
    SimSettings.loadname[0] = 0;
    SimSettings.sampling_frequency = 100000; //Hz Dalanco_max
    SimSettings.filter_frequency = 10000; // Hz Mess_maxzu Dalanco_max
    SimSettings.noise_length =1000000;     //default 1M
    SimSettings.actionpotentials = false;
    SimSettings.repeats = 1;							
    SimSettings.markercounts = 0;
    SimSettings.resample = true;
    SimSettings.dosbox = true;  	//Toby 2013
    SimSettings.show_2D_Diff = false; 	//Toby 2013

    // SimSettings.noisefile[0] = 0; //Tobias 2022
    sprintf(SimSettings.check, "error");	//Toby 2013
    SimSettings.write_result_file = true;
    sprintf(SimSettings.result_file_name, "fit_results.log");

       
    
    Fitparameter.spannweite =    5000;
    Fitparameter.ftol =          1;
    Fitparameter.minmultiplikator = 1;
    Fitparameter.maxmultiplikator = 1;
    Fitparameter.multiplikator = 1;
    Fitparameter.maxiterations = 1000;
    Fitparameter.log_min_close = 0;
    Fitparameter.log_max_close = 5;
    Fitparameter.log_min_open  = 0;
    Fitparameter.log_max_open  = 5;
    Fitparameter.bins_per_log  = 10;
    
    Geneticfitparameter.logscale = true;
    Geneticfitparameter.fit_methode = 1; //default = 1, 2D fit without amplitude fit
    Geneticfitparameter.popsize  = 1000;
    Geneticfitparameter.ngen     = 10000;
    Geneticfitparameter.pmut   = 0.01;
    Geneticfitparameter.pcross = 0.6;
    Geneticfitparameter.bits = 16;
    Geneticfitparameter.nConvergence = 1;	//Generations to evaluate before fit without improvement terminates 
    Geneticfitparameter.methode = 1;
    Geneticfitparameter.likelihood_2d = 0;			//2D_error global f�r �bergabe
    Geneticfitparameter.histogram_diff = 0;			//amplitude histogram error f�r �bergabe
    for (unsigned int i=1; i<=100; i++)
     {
      Geneticfitparameter.startwertmin[i] = 10000;
      Geneticfitparameter.startwertmax[i] = 10000;	
     }	
    Geneticfitparameter.proportionalfit = false; 	//Tobias 2013 Levelfit 
    Geneticfitparameter.levelfit = false;	 	//Tobias 2012 Levelfit
    Geneticfitparameter.stromvektor[1] = 0;	 	//Tobias 2013 Levelfit
    Geneticfitparameter.stromvektor[2] = 0;
    Geneticfitparameter.stromvektor[3] = 0;
    Geneticfitparameter.stromvektor[4] = 0;
    Geneticfitparameter.stromvektor[5] = 0;
    Geneticfitparameter.stromvektor[6] = 0;
    Geneticfitparameter.stromvektor[7] = 0;
    Geneticfitparameter.stromvektor[8] = 0;
    Geneticfitparameter.stromvektor[9] = 0;
    Geneticfitparameter.stromvektor[10] = 0;
    Geneticfitparameter.state_conductance[1] = 0;
    Geneticfitparameter.state_conductance[2] = 0;
    Geneticfitparameter.state_conductance[3] = 0;
    Geneticfitparameter.state_conductance[4] = 0;
    Geneticfitparameter.state_conductance[5] = 0;
    Geneticfitparameter.state_conductance[6] = 0;
    Geneticfitparameter.state_conductance[7] = 0;
    Geneticfitparameter.state_conductance[8] = 0;
    Geneticfitparameter.state_conductance[9] = 0;
    Geneticfitparameter.state_conductance[10] = 0;
    Geneticfitparameter.levelfit_min[0] = 0;
    Geneticfitparameter.levelfit_min[0] = 0;
    Geneticfitparameter.levelfit_min[1] = 150;
    Geneticfitparameter.levelfit_max[1] = 250;
    Geneticfitparameter.levelfit_min[2] = 150;
    Geneticfitparameter.levelfit_max[2] = 250;
    
    Geneticfitparameter.automated_loading = false;			//Tobias 2013 Levelfit
    Geneticfitparameter.files_directory[0] = 0;		//Tobias 2013 Levelfit
    Geneticfitparameter.file_time_series[0] = 0;		//Tobias 2013 Levelfit
    Geneticfitparameter.file_noise_series[0] = 0;	//Tobias 2013 Levelfit
    Geneticfitparameter.file_noise_series_loaded = false; //Tobias 2022  	
    Geneticfitparameter.file_setfile[0] = 0;				//Tobias 2013 Levelfit
    Geneticfitparameter.nr_calc_mean_result = 10;				//Tobias 2013
    Geneticfitparameter.elitism = false;								//Tobias 2014
    Geneticfitparameter.init_random = 10;								//Tobias 2018 0=randomize timer /10=default
   
    
    Geneticfitparameter.mpi_support = false;
            
    sprintf(logfilename,"default.log");
    sprintf(logdwellfilename,"DwellHinkley.log");
    
    int i;
    
    dependency_plot_draw_diff = false;
		matlab_file_loaded = false;						//Tobias 2022 to prevent loading of matlab_file after first load
		double* stepresponse_loaded;
		Stuetzpunkte_loaded = 0;	 

    
