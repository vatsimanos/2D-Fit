/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  2000 Tobias Huth 
 *  based on uwe kirst's daytest  
 *  CAU KIEL
    FAU Erlangen-Nürnberg
 *  Genetic Fit: Copyright 1995-1996 Massachusetts Institute of Technology
 *******************************************************************/
  
#include "sim_fit.h"
#include "sim_fit_dialogs.c"
#include <thread>		//2019
#include <chrono>
#include <cstring>
#include "MultiLevel2D.h"



//**************************************************************************************************************************
void ZweiDFit()
{
}

//**************************************************************************************************************************
double simulation(int multiplikator, int nr_const, double *rconst, double rtol)
 {		
  TMultiLevel2D multilevel2D(Daten.starting_level, Daten.level_increment, Daten.number_of_levels);		
  char textbuf[MaxTextLen];
  int j,i;
  int ii;
  int counts;
  double xtol;

  double ampl_simul[65536];
  double histogram_diff;
  double likelihood_2d;
  double error = 0;	 

 Setfile.iterations++;

if (Daten.Geneticfitparameter.init_random != 0)
  {
   Daten.random_engine.seed(Daten.Geneticfitparameter.init_random);
  } 
  
  //for(i=1; i<= nr_const; i++)
    //{
    //if ( (rconst[i] < Daten.Geneticfitparameter.startwertmin[i]) || (rconst[i] > Daten.Geneticfitparameter.startwertmax[i]) )
    // 	rconst[i] = Daten.Geneticfitparameter.startwertmin[i]+KISS%int(Daten.Geneticfitparameter.startwertmax[i]-Daten.Geneticfitparameter.startwertmin[i]);;
    //Setfile.set_rconst(i, rconst[i]);
    //}  	
  
  xtol = rtol; //rtol bei startwerten = 0
  if (xtol == 0)
    xtol = -1;
 
  //Berechnung der neuen Anzahl der Simulationen aus rtol(Konvergenz), minimalem und maximalen Multiplikator
  //if (multiplikator == 0)
  //  {
  //    Daten.Fitparameter.multiplikator= (int) (Daten.Fitparameter.maxmultiplikator * Daten.Fitparameter.ftol*10 /rtol);
  //    if (Daten.Fitparameter.multiplikator > Daten.Fitparameter.maxmultiplikator)
  //   		Daten.Fitparameter.multiplikator = Daten.Fitparameter.maxmultiplikator;
  //    if (Daten.Fitparameter.multiplikator < Daten.Fitparameter.minmultiplikator)
  //      Daten.Fitparameter.multiplikator = Daten.Fitparameter.minmultiplikator;  
  //   }   
  // else
  //   Daten.Fitparameter.multiplikator = multiplikator;
  

 
  //initialisierung der 2D-Matrix	
  Daten.Dwell_2d_ptr->init((char *)"simulated", Daten.Fitparameter.log_min_close,
  				        Daten.Fitparameter.log_max_close, 
  				        Daten.Fitparameter.log_min_open, 
  				        Daten.Fitparameter.log_max_open, 
  				        Daten.Fitparameter.bins_per_log);
				        
  //Initialisierung	
  Daten.reset();	
  Daten.Wert = new unsigned short[Setfile.get_no_samples()];
  Daten.TwoD_fitinaction = true;

  // Anzahl der Wiederholungen(Marker in der Zeitreihe * Multiplikator f�r die Genauigkeit 
  counts = multiplikator * Daten.SimSettings.markercounts; 	

  if (counts == 0)
    counts = multiplikator;

  for (ii=1; ii <= counts; ii++)	 
   {
	

     //�bergeben der PARAMETER aus dem Setfile 
     Daten.Anz_Samples = Setfile.get_no_samples();
     Daten.a_fit.sigma = Setfile.get_sigma();
     Daten.a_fit.i_null = Setfile.get_zero_current();
     Daten.a_fit.i_channel = Setfile.get_single_channel_current();
     Daten.a_fit.n_channels = Setfile.get_no_channels_total();
     Daten.Settings.SamplesToSkip = 0;
     Daten.Settings.SamplesToProcess = Setfile.get_no_samples();
     
    //std::cout<<"Daten.Anz_Samples :"<<Daten.Anz_Samples<<std::endl;
    // std::cout<<"Daten.a_fit.sigma :"<<Daten.a_fit.sigma<<std::endl;
    // std::cout<<"Daten.a_fit.i_null :"<<Daten.a_fit.i_null<<std::endl;
    // std::cout<<"Daten.a_fit.i_channel :"<<Daten.a_fit.i_channel<<std::endl;
    // std::cout<<"Daten.a_fit.n_channels :"<<Daten.a_fit.n_channels<<std::endl;
    // std::cout<<"Daten.Settings.SamplesToSkip :"<<Daten.Settings.SamplesToSkip<<std::endl;
    // std::cout<<"Daten.Settings.SamplesToProcess :"<<Daten.Settings.SamplesToProcess<<std::endl;     
     
     //simuliert die Zeitreihe
     Setfile.simulate_timeseries(&TobyTimeseries,&Rausch, (char *)"", 0);
     //for (j=1; j <= Setfile.get_no_samples(); j++)
     //   Daten.Wert[j-1] = Daten.SimulatWert[j];  //Wert startet bei 0, SimulatWert bei 1
    	     
     memcpy( &Daten.Wert[0], &Daten.SimulatWert[1], Setfile.get_no_samples()*sizeof(Daten.Wert[1]));
     
     //for (int j=0; j < Setfile.get_no_samples(); j++)
     //    std::cout<<Daten.Wert[j]<<std::endl;
    
  
     //setzten der number Anfangswerte auf offset/value zum fitten gemessener Daten
     if (Daten.loeschanzahl > 0)
       for (j=0; j < Daten.loeschanzahl; j++)
         Daten.Wert[j] = Daten.loeschvalue;
     // ruft den HOHD auf
     Daten.hinkley(_HOHD);
     //addiert die bins zum 2D
     Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);   
   }
  

  Setfile.total_nr_events = Setfile.total_nr_events + Daten.Dwell_2d_ptr->events();
  //std::cout<<"Daten.Dwell_2d_ptr->events() :"<<Daten.Dwell_2d_ptr->events()<<std::endl; 
        
  //Tobias 2011 Berechnung der Amplitutenhistogram-Abweichung
  
  //initial
  histogram_diff = 0;
  for (i=0;i < 65536; i++)
        ampl_simul[i] = 0;
  
  for (i=0;i < Daten.Anz_Samples;i++)
        ampl_simul[Daten.Wert[i]]++;

   //double temp = 0;
  for (i=0;i < 65536; i++){								//2014 Test level default 0 -> <65536
        histogram_diff = histogram_diff +(ampl_simul[i] - Daten.original_ampl_histogram[i])*(ampl_simul[i] - Daten.original_ampl_histogram[i]);
        //temp = temp +Daten.original_ampl_histogram[i];
  }
         //std::cout<<"sum:"<<temp<<std::endl; 
  

  if (Daten.Geneticfitparameter.fit_methode == 4) 	//create multi level 2d-histograms of experimental time series
 		{
 	 	 likelihood_2d = multilevel2D.return_2d_matrix_LLH();	
 	 	 error = likelihood_2d;
       //std::cout<<"error :"<<error<<std::endl; 
    }	 
  else
  	{
     likelihood_2d = Daten.Dwell_2d_A.lnlikelihood(Daten.Dwell_2d_B);  
     //std::cout<<"likelihood_2d :"<<likelihood_2d<<std::endl; 
    }
    
  if (Daten.Geneticfitparameter.fit_methode == 1)
    error = likelihood_2d;
   //std::cout<<"error :"<<error<<std::endl; 

  if (Daten.Geneticfitparameter.fit_methode == 2)
    error = histogram_diff;
    
  if (Daten.Geneticfitparameter.fit_methode == 3)
  	 if (Daten.calculate_error == 0){
      //std::cout<<"likelihood_2d :"<<likelihood_2d<<std::endl;
      //td::cout<<"histogram_diff :"<<histogram_diff<<std::endl; 
      //double factor = std::pow(10,8); //Efthymios 2024 to make the two loss scores comparable
  	 	 //error = likelihood_2d + factor/ histogram_diff;
       error = 1 + (likelihood_2d / (2*Daten.Anz_Samples*log(2*Daten.Anz_Samples))) - (histogram_diff / (2*Daten.Anz_Samples*Daten.Anz_Samples)); //Efthymios 2024 normalized individual losses to make them comparable
       //std::cout<<"upper :"<<(likelihood_2d / (2*Daten.Anz_Samples*log(2*Daten.Anz_Samples)))<<std::endl;  
       //std::cout<<"lower :"<<(histogram_diff / (2*Daten.Anz_Samples*Daten.Anz_Samples))<<std::endl;  
       //std::cout<<"error :"<<error<<std::endl;  

      }								// Tobias 2018 maximizing likelihood_2d now 
  	 else	  
       error = likelihood_2d * histogram_diff;						   //Tobias 2018 other option: log(likelihood_2d) + log(histogram_diff);

       
       
    
  Daten.Geneticfitparameter.likelihood_2d = likelihood_2d;			//2D_error global f�r �bergabe
  Daten.Geneticfitparameter.histogram_diff = histogram_diff;  
    
  
 
 //Ausgabe der Ratenkonstanten/Ergebniss in die DOS-Box 
 unsigned int parameter;
  if (Daten.SimSettings.dosbox)
   {
    std::cout<<"iterations: "<<std::setw(5)<<Setfile.iterations<<" multiplikator "<< Daten.Fitparameter.multiplikator<<" ";
    for (i=1; i<=Setfile.get_no_rconst(); i++)
         {						
         	parameter = round(Setfile.get_rconst(i));																			//sprintf(textbuf,"%5.0f ",floor(Setfile.get_rconst(i)+0.5));
 	        std::cout<<std::setw(6)<<parameter;																						//std::cout<<textbuf;
         }  
    if (Daten.Geneticfitparameter.levelfit == true)			//Tobias 2013 Levelfit
          {
           std::cout<<"l1:"<<floor(Daten.Geneticfitparameter.level1)<<" ";
           if (Daten.Geneticfitparameter.levelfit_max[2] != 0)
           		std::cout<<"l2:"<<floor(Daten.Geneticfitparameter.level2)<<" ";
          }  
	
    std::cout<<std::setprecision(10)<<"  2D: "<<std::setw(10)<<likelihood_2d<<" Ampl: "<<std::setw(10)<<histogram_diff<<" fit: "<<std::setw(10)<<error<<" "<<std::setprecision(4);  
    std::cout<< std::endl;
   }
              
	
  return error;
 }
 

//**************************************************************************************************************************
void out_time()
 {
  //double seconds = 0;	
  std::ofstream logout(Daten.logfilename,std::ios::app);	 	
 	
  //seconds = 0;//Tobias 2021 g_timer_elapsed (timer, micros);
  
  auto timer_end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = timer_end - timer_start;
  
  std::cout <<"time/s "<<elapsed.count()<<std::endl;	
  if (Daten.SimSettings.dosbox == true)	
  	std::cout <<std::endl;
 }
 

//***********************************************************************************************************************************************************************************************************
unsigned int GeneticFit_SimpleGA()
 {
  int i;
  
  //Setzen von BIN B (simulierte Daten)
   Daten.Dwell_2d_ptr = &Daten.Dwell_2d_B; 
   Daten.Dwell_1d_ptr = &Daten.Dwell_1d_B;  
   Daten.Dwell_2d_ptr->init((char *)"simulated", Daten.Fitparameter.log_min_close,
  				         Daten.Fitparameter.log_max_close, 
  				         Daten.Fitparameter.log_min_open, 
  				         Daten.Fitparameter.log_max_open, 
  				         Daten.Fitparameter.bins_per_log);

  GABin2DecPhenotype map;
 
  for(i=1; i<=Setfile.get_no_rconst();i++)
    { 
     if (Daten.Geneticfitparameter.logscale == true)
  	    map.add(Daten.Geneticfitparameter.bits, std::log(Daten.Geneticfitparameter.startwertmin[i]), std::log(Daten.Geneticfitparameter.startwertmax[i]));
     else	
        map.add(Daten.Geneticfitparameter.bits, Daten.Geneticfitparameter.startwertmin[i], Daten.Geneticfitparameter.startwertmax[i]);
    }
// 2012 Tobias leveldetection here - currently only state 2

  if (Daten.Geneticfitparameter.levelfit == true) 
  	{
  	 map.add(Daten.Geneticfitparameter.bits,Daten.Geneticfitparameter.levelfit_min[1],Daten.Geneticfitparameter.levelfit_max[1]); 	//1. open level
  	 if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
  	    map.add(Daten.Geneticfitparameter.bits,Daten.Geneticfitparameter.levelfit_min[2],Daten.Geneticfitparameter.levelfit_max[2]);  //2. open level
  	}	
  
  
  
//Create the template genome using the phenotype map we just made.
  GABin2DecGenome genome(map, objective);

// Now create the GA using the genome and run it.  We'll use sigma truncation
// scaling so that we can handle negative objective scores.      
  GASimpleGA ga(genome); 
  
  GASigmaTruncationScaling scaling;
  ga.populationSize(Daten.Geneticfitparameter.popsize);
  ga.nGenerations(Daten.Geneticfitparameter.ngen);
  ga.pMutation(Daten.Geneticfitparameter.pmut);
  ga.pCrossover(Daten.Geneticfitparameter.pcross);
  ga.scaling(scaling);
  ga.nConvergence(Daten.Geneticfitparameter.nConvergence);
  ga.initialize();

 while( ((ga.statistics().convergence() == 0) || (ga.statistics().convergence()>1)) && (ga.statistics().generation() < Daten.Geneticfitparameter.ngen) )
   {
    generation_out_2Dfit(ga, genome);
    ++ga;
   }

 summary_out_2Dfit(ga,genome);   
 return 0;
} 

//*************************************************************************************************************************************** 

float objective(GAGenome & c)
{
 int i;
 double *Ratenkonstanten;	
 double ergebnis;
 
 GABin2DecGenome & genome = (GABin2DecGenome &)c;

 Ratenkonstanten = new double[Setfile.get_no_rconst()+1];

 for (i=1; i<=Setfile.get_no_rconst(); i++)
   {
    if (Daten.Geneticfitparameter.logscale == true)
       Ratenkonstanten[i]=std::exp(genome.phenotype(i-1));	 
    else	
       Ratenkonstanten[i]=genome.phenotype(i-1);
   }
 for (i=1; i <= Setfile.get_no_rconst(); i++)
    {   
     Setfile.set_rconst(i,Ratenkonstanten[i]);   
    }
 if (Daten.Geneticfitparameter.levelfit == true)  //Tobias 2013 Levelfit
     for (i=1; i <= 10; i++)
       {
       	if (Daten.Geneticfitparameter.state_conductance[i] == 0)
       	  Daten.Geneticfitparameter.stromvektor[i] = 0;
       	if (Daten.Geneticfitparameter.state_conductance[i] == 1)
       	  Daten.Geneticfitparameter.stromvektor[i] = genome.phenotype(Setfile.get_no_rconst());
       	if (Daten.Geneticfitparameter.state_conductance[i] == 2)
       	  Daten.Geneticfitparameter.stromvektor[i] = genome.phenotype(Setfile.get_no_rconst()+1);  
       	  
       	Daten.Geneticfitparameter.level1 = genome.phenotype(Setfile.get_no_rconst());
       	if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
       		Daten.Geneticfitparameter.level2 = genome.phenotype(Setfile.get_no_rconst()+1);  
       }	
 
 ergebnis = simulation(Daten.Fitparameter.multiplikator, 0, Ratenkonstanten, 0);
 delete[] Ratenkonstanten;
 return -ergebnis;
}

//***************************************************************************************************************************************

void init_2Dfit()
{
  TMultiLevel2D multilevel2D(Daten.starting_level, Daten.level_increment, Daten.number_of_levels);	
  int i;
  int nextposition;
  char textbuf[MaxTextLen];
  double version;
  TDoublefeld matlab_test(0,(char *)"stepresponse");
  int test;
  char* test_filename;
  Daten.SimSettings.resample = true;
  Daten.SimSettings.noise_length =  Daten.Anz_Samples;
  //Zeitreihe vorhanden ??
  if (Daten.Anz_Samples == 0)
    {
     std::cout<<"No timeseries found"<<std::endl;	
     sprintf(textbuf,"No timeseries found");
     Application.statustext(textbuf);
     return;
    } 

//starten des Timers
   //Tobias 2021 timer = g_timer_new();

  strcpy(textbuf,Daten.SimSettings.loadname);
  strcat(textbuf,".log");
  strcpy(Daten.logfilename,textbuf);

 
  //laden des Setfiles!          
  version = Setfile.read_setfile(Daten.SimSettings.loadname);
      if (version == 0)
       {
        std::cout<<Daten.errormessage<<std::endl;
        strcpy(textbuf,Daten.errormessage);
        Application.statustext(textbuf);
        return;
       }	 
   

  test_filename = Setfile.get_stepresponse_filename(NULL); 
  test  = matlab_test.load_matlab_file(test_filename);     
  if (test == 1)
    {
     std::cout<<Daten.errormessage<<std::endl;
     //strcpy(textbuf,Daten.errormessage);
     //Application.statustext(textbuf);
     return;
    }   	

  if (Daten.Anz_Samples != Setfile.get_no_samples())
  	{
  	 std::cout<<" length of loaded time series and length of time series in setfile do not match, aborting!----------------------->"<<std::endl;
  	 exit(0);		
    }  
//initialisierung 
   Daten.SimulatWert = new unsigned short[Setfile.get_no_samples()+1];  //reserviert den Speicher f�r die simulierte Zeitreihe   
   Setfile.total_nr_events = 0;      
  
   
//�bergeben der PARAMETER aus dem Setfile 

   Daten.a_fit.sigma = Setfile.get_sigma();
   Daten.a_fit.i_null = Setfile.get_zero_current();
   Daten.a_fit.i_channel = Setfile.get_single_channel_current();
   Daten.a_fit.n_channels = Setfile.get_no_channels_total();
   
   
//Zeitreihe im Speicher auswerten wenn noch nicht geschehen 
    	
       Daten.TwoD_fitinaction = true;
       Daten.Dwell_2d_ptr = &Daten.Dwell_2d_A;  //Setzen von BIN A (gemessenen Daten) 
       Daten.Dwell_1d_ptr = &Daten.Dwell_1d_A;
  
       Daten.Dwell_2d_ptr->init((char *)"measured", Daten.Fitparameter.log_min_close,
        				    Daten.Fitparameter.log_max_close, 
   	        			    Daten.Fitparameter.log_min_open, 
  				            Daten.Fitparameter.log_max_open, 
  				            Daten.Fitparameter.bins_per_log);
//Auswahl kontinuierliche Zeireihe <-> markierte Einzelreihen

       if (Daten.SimSettings.actionpotentials == false)
         {
          Daten.Settings.SamplesToSkip = 0;
          Daten.Settings.SamplesToProcess = Daten.Anz_Samples;        
          Daten.hinkley(_HOHD);

   
//addiert die bins zum 2D
          Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);

         }
       else //Begin Schleife f�r markierte Zeitreihen
         {
          for (Daten.lauf = Daten.Marker.begin(); Daten.lauf !=Daten.Marker.end(); Daten.lauf++)
             {	
              if(Daten.lauf+1 !=Daten.Marker.end()) //n�tig f�r das korrekte Ende der Zeitreihe
                   nextposition = (Daten.lauf+1)->startpos;
                else
                   nextposition = Daten.Anz_Samples;	
      
              Daten.Settings.SamplesToSkip = Daten.lauf->startpos;
              Daten.Settings.SamplesToProcess = nextposition - Daten.lauf->startpos;        
              if (Daten.Settings.SamplesToProcess > 2)  //nur Hinkley, wenn Me�datenabschnitt, nicht bei Zwischenabschenitt
                {  
                 Daten.hinkley(_HOHD);
//addiert die bins zum 2D
                 Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);
                 Daten.SimSettings.markercounts +=1; 
                } 
                
             }  //end Schleife
      
        }  
      for (i=0;i < 65536;i++)
        Daten.original_ampl_histogram[i]=0;          
      for (i=0;i < Daten.Anz_Samples;i++)
       {
        Daten.original_ampl_histogram[Daten.Wert[i]]++;
       }        
       
 if (Daten.Geneticfitparameter.fit_methode == 4) 	//create multi level 2d-histograms of experimental time series
 	{
 	 multilevel2D.generate_matrix_2D_exp();	
  }	       
       
       
 //Daten.Fitparameter.multiplikator = Daten.Fitparameter.minmultiplikator; 

}
//*************************************************************************************************************************************** 
void
out_2D_settings()
{
	
  int i;	
  char textbuf[MaxTextLen];
  std::ofstream logout(Daten.logfilename,std::ios::app);
	
//Ausgabe der Systemzeit
   time_t rawtime;
   time ( &rawtime );
   sprintf ( textbuf,"%s", ctime (&rawtime) );
  
   std::cout<<"***************************************************************"<<std::endl;
   logout<<"***************************************************************"<<std::endl; 
   
   std::cout<<"start 2D-fit at: "<<textbuf<<std::endl;
   logout<<"start 2D-fit at: "<<textbuf<<std::endl;	
   
   std::cout<<"timeseries: "<<Daten.FName<<std::endl;
   std::cout<<"simulating setfile: "<<Daten.SimSettings.loadname<<std::endl;
   std::cout<<"simulating noisefile: "<<Daten.Geneticfitparameter.file_noise_series<<std::endl;		//Tobias 2022
   std::cout<<"fit settings:***************************************************************"<<std::endl; 
   std::cout<<"simulating timeseries "<<Daten.SimSettings.repeats<<" times"<<std::endl;
   std::cout<<"sampling Frequency set to: "<<Daten.SimSettings.sampling_frequency<<" kHz"<<std::endl;  
   std::cout<<"filter Frequency set to:   "<<Daten.SimSettings.filter_frequency<<" kHz"<<std::endl; 
   std::cout<<"multiplikator min /max: "<<Daten.Fitparameter.minmultiplikator<<"/"<<Daten.Fitparameter.maxmultiplikator<<std::endl;
   std::cout<<"actionpotentials: "<<Daten.SimSettings.actionpotentials<<std::endl;
   std::cout<<"Set first "<<Daten.loeschanzahl<<" points to "<<Daten.loeschvalue<<std::endl;
   std::cout<<"2D-histogram settings:***************************************************************"<<std::endl; 
   std::cout<<"error calculation( 0=LogLikelihood., 1=1bin error, 2=filter 9bin error, 3=filter 25binerror, 4= multi level / multiple 2D-histograms): "<<Daten.calculate_error<<std::endl;
   std::cout<<"Bpl: "<<Daten.Fitparameter.bins_per_log<<std::endl;
   std::cout<<"Close: "<<Daten.Fitparameter.log_min_close<<" to "<<Daten.Fitparameter.log_max_close<<std::endl;
   std::cout<<"Open: "<<Daten.Fitparameter.log_min_open<<" to "<<Daten.Fitparameter.log_max_open<<std::endl;  
   std::cout<<"Show 2D-histogram (0=false; 1=true): "<<Daten.SimSettings.show_2D_Diff<<std::endl; 
   std::cout<<"algorithm settings:***************************************************************"<<std::endl; 

   logout<<"timeseries: "<<Daten.FName<<std::endl;
   logout<<"simulating setfile: "<<Daten.SimSettings.loadname<<std::endl;
   logout<<"simulating noisefile: "<<Daten.Geneticfitparameter.file_noise_series<<std::endl;			//Tobias 2022
   logout<<"fit settings:***************************************************************"<<std::endl; 
   logout<<"simulating timeseries "<<Daten.SimSettings.repeats<<" times"<<std::endl;
   logout<<"sampling Frequency set to: "<<Daten.SimSettings.sampling_frequency<<" kHz"<<std::endl;  
   logout<<"filter Frequency set to:   "<<Daten.SimSettings.filter_frequency<<" kHz"<<std::endl; 
   logout<<"multiplikator min /max: "<<Daten.Fitparameter.minmultiplikator<<"/"<<Daten.Fitparameter.maxmultiplikator<<std::endl;
   logout<<"actionpotentials: "<<Daten.SimSettings.actionpotentials<<std::endl;
   logout<<"Set first "<<Daten.loeschanzahl<<" points to "<<Daten.loeschvalue<<std::endl;
   logout<<"2D-histogram settings:***************************************************************"<<std::endl; 
   logout<<"error calculation( 0=LogLikelihood., 1=1bin error, 2=filter 9bin error, 3=filter 25binerror, 4= multi level / multiple 2D-histograms): "<<Daten.calculate_error<<std::endl;
   logout<<"Bpl: "<<Daten.Fitparameter.bins_per_log<<std::endl;
   logout<<"Close: "<<Daten.Fitparameter.log_min_close<<" to "<<Daten.Fitparameter.log_max_close<<std::endl;
   logout<<"Open: "<<Daten.Fitparameter.log_min_open<<" to "<<Daten.Fitparameter.log_max_open<<std::endl; 
   logout<<"Show 2D-histogram (0=false; 1=true): "<<Daten.SimSettings.show_2D_Diff<<std::endl;   
   logout<<"algorithm settings:***************************************************************"<<std::endl; 
   

    std::cout<<"error calculation (1=2D, 2=amplitude histogram, 3=2D*amplitude histogram): "<<Daten.Geneticfitparameter.fit_methode<<std::endl; 
    std::cout<<"proportional-Fit of rate konstants (0=false, 1=true): "<<Daten.Geneticfitparameter.proportionalfit<<std::endl;
    std::cout<<"population size:             "<<Daten.Geneticfitparameter.popsize<<std::endl;
    std::cout<<"generations for Convergence: "<<Daten.Geneticfitparameter.nConvergence<<std::endl;
    std::cout<<"elitism setting(1 = yes, 0 = No): "<<Daten.Geneticfitparameter.elitism<<std::endl;
   	std::cout<<"maximum generations:         "<<Daten.Geneticfitparameter.ngen<<std::endl;
   	std::cout<<"mutation rate:  "<<Daten.Geneticfitparameter.pmut<<std::endl;
   	std::cout<<"crossover rate: "<<Daten.Geneticfitparameter.pcross<<std::endl; 
  	std::cout<<"bits / rate constant: "<<Daten.Geneticfitparameter.bits<<std::endl;
  	std::cout<<"log scaling of rate  constants(1 = yes, 0 = No): "<<Daten.Geneticfitparameter.logscale<<std::endl;
  	std::cout<<"times simulation for error calculation n = "<<Daten.Geneticfitparameter.nr_calc_mean_result<<std::endl;
    	
   	logout<<"error calculation (1=2D, 2=amplitude histogram, 3=2D*amplitude histogram): "<<Daten.Geneticfitparameter.fit_methode<<std::endl;
       	logout<<"proportional-Fit of rate konstants (0=false, 1=true): "<<Daten.Geneticfitparameter.proportionalfit<<std::endl;
    	logout<<"population size:             "<<Daten.Geneticfitparameter.popsize<<std::endl;
    	logout<<"generations for Convergence: "<<Daten.Geneticfitparameter.nConvergence<<std::endl;
    	logout<<"elitism setting(1 = yes, 0 = No): "<<Daten.Geneticfitparameter.elitism<<std::endl;
   	logout<<"maximum generations:         "<<Daten.Geneticfitparameter.ngen<<std::endl;
   	logout<<"mutation rate:  "<<Daten.Geneticfitparameter.pmut<<std::endl;
   	logout<<"crossover rate: "<<Daten.Geneticfitparameter.pcross<<std::endl;
  	logout<<"bits / rate constant: "<<Daten.Geneticfitparameter.bits<<std::endl;
  	logout<<"log scaling of rate  constants(1  = yes, 0 = No): "<<Daten.Geneticfitparameter.logscale<<std::endl;
  	logout<<"times simulation for error calculation n = "<<Daten.Geneticfitparameter.nr_calc_mean_result<<std::endl;

   	for (i=1; i<=Setfile.get_no_rconst(); i++)
   	    {
   	    std::cout<<"rate constant k"<<i<<" set to min:"<<Daten.Geneticfitparameter.startwertmin[i]<<" max:"<<Daten.Geneticfitparameter.startwertmax[i]<<std::endl;
   	    logout<<"rate constant k"<<i<<" set to min:"<<Daten.Geneticfitparameter.startwertmin[i]<<" max:"<<Daten.Geneticfitparameter.startwertmax[i]<<std::endl;
   	    }
   	if (Daten.Geneticfitparameter.levelfit == true)
   	    {
   	     std::cout<<"levelfit activated***************************************************************"<<std::endl;
   	     std::cout<<"level conductance matrix: ";
   	     for (i=1; i<=Setfile.get_matrix_dimension(); i++)
   	       std::cout<<Daten.Geneticfitparameter.state_conductance[i]<<"  ";
   	     std::cout<<std::endl;  
   	     std::cout<<"level 1 conductance set to min: "<<Daten.Geneticfitparameter.levelfit_min[1]<<" max:"<<Daten.Geneticfitparameter.levelfit_max[1]<<std::endl;
   	     if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
   	     	std::cout<<"level 2 conductance set to min: "<<Daten.Geneticfitparameter.levelfit_min[2]<<" max:"<<Daten.Geneticfitparameter.levelfit_max[2]<<std::endl;
   	     
   	     logout<<"levelfit activated***************************************************************"<<std::endl;
   	     logout<<"level conductance matrix: ";
   	     for (i=1; i<=Setfile.get_matrix_dimension(); i++)
   	       logout<<Daten.Geneticfitparameter.state_conductance[i]<<"  ";
   	     logout<<std::endl;  
   	     logout<<"level 1 conductance set to min: "<<Daten.Geneticfitparameter.levelfit_min[1]<<" max:"<<Daten.Geneticfitparameter.levelfit_max[1]<<std::endl;
   	     if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
   	     	logout<<"level 2 conductance set to min: "<<Daten.Geneticfitparameter.levelfit_min[2]<<" max:"<<Daten.Geneticfitparameter.levelfit_max[2]<<std::endl;
        }
   	
   	   
       
   std::cout<<std::endl;
   logout<<std::endl;
   std::cout<<"Nr of intervalls in timeseries: "<<Daten.SimSettings.markercounts<<std::endl;
   logout<<"Nr of intervalls in timeseries: "<<Daten.SimSettings.markercounts<<std::endl;	
   
}


//*************************************************************************************************************************************** 
void
generation_out_2Dfit(GASimpleGA & ga, GABin2DecGenome & genome)
{
 int i;
 std::ofstream logout(Daten.logfilename,std::ios::app);
 	if (Daten.SimSettings.dosbox == true)
    		std::cout<<"*********************************************************************"<<std::endl;	
 genome = ga.statistics().bestIndividual();

    std::cout<<"gen "<<std::setw(4)<<ga.statistics().generation(); 
    std::cout<<" conv "<<std::setprecision(5)<<std::setw(6)<<ga.statistics().convergence()<<" / "<<std::setw(2)<<ga.statistics().nConvergence();
    std::cout<<" err "<<std::scientific<<std::setprecision(6)<<std::setw(13)<<genome.score()<<std::defaultfloat;
    std::cout<<" kij ";	
    

    for (i=1; i<=Setfile.get_no_rconst(); i++)
      if (Daten.Geneticfitparameter.logscale == true)
         std::cout<<std::setw(6)<< round(exp(genome.phenotype(i-1)))<<" ";
       else  
         std::cout<<std::setw(6)<<round(genome.phenotype(i-1))<<" ";

         	
    if (Daten.Geneticfitparameter.levelfit == true)			// Tobias 2013 Levelfit
	{
         Daten.Geneticfitparameter.level1 = genome.phenotype(Setfile.get_no_rconst());
         if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
       	 	Daten.Geneticfitparameter.level2 = genome.phenotype(Setfile.get_no_rconst()+1); 
    	 std::cout<<" fittet lvl1: "<<Daten.Geneticfitparameter.level1<<" ";
    	 if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
    	 	std::cout<<" fittet lvl2: "<<Daten.Geneticfitparameter.level2;
    	} 
   std::cout<<" ";
     
    out_time();	
}
//*************************************************************************************************************************************** 
void
summary_out_2Dfit(GASimpleGA & ga,GABin2DecGenome & genome)
{
 char textbuf[MaxTextLen];	
 int i;
 int ii;
 double *Ratenkonstanten;
 Ratenkonstanten = new double[Setfile.get_no_rconst()+3];	 
 int ampl[4096]; 
 int init_nr;
 
 double result[4];
 double ergebnis;
 double sigma;
 MPI_Status status;
 
 double *erg;
 erg = new double[Daten.Geneticfitparameter.nr_calc_mean_result+1];
 double *error_likelihood_2d;
 error_likelihood_2d = new double[Daten.Geneticfitparameter.nr_calc_mean_result+1];
 double *error_histogram_diff;
 error_histogram_diff = new double[Daten.Geneticfitparameter.nr_calc_mean_result+1];
 
 double ergebnis_likelihood_2d;
 double sigma_likelihood_2d;
 double ergebnis_histogram_diff;
 double sigma_histogram_diff;

  
//besten (kleinste Evaluation) Parametersatz ausgeben und 2D-Histogramm anzeigen
  std::cout<<"simulating best result ("<<Daten.Geneticfitparameter.nr_calc_mean_result<<"x):"<<std::endl;
 
  genome = ga.statistics().bestIndividual();
 
  for (i=1; i<=Setfile.get_no_rconst(); i++)
    {
    if (Daten.Geneticfitparameter.logscale == true)
      Ratenkonstanten[i]=exp(genome.phenotype(i-1)); 
    else	
      Ratenkonstanten[i]=genome.phenotype(i-1);    	
    }
  for (i=1; i <= Setfile.get_no_rconst(); i++)      
    if (Daten.Geneticfitparameter.logscale == true)
      Setfile.set_rconst(i,exp(genome.phenotype(i-1))); 
    else
      Setfile.set_rconst(i,genome.phenotype(i-1)); 


  if (Daten.Geneticfitparameter.levelfit == true)  //Tobias 2013 Levelfit
    {
     Ratenkonstanten[Setfile.get_no_rconst()+1] = genome.phenotype(Setfile.get_no_rconst());
     if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
     	Ratenkonstanten[Setfile.get_no_rconst()+2] = genome.phenotype(Setfile.get_no_rconst()+1);	
     for (i=1; i <= 10; i++)
         {
       	  if (Daten.Geneticfitparameter.state_conductance[i] == 0)
       	    Daten.Geneticfitparameter.stromvektor[i] = 0;
       	  if (Daten.Geneticfitparameter.state_conductance[i] == 1)
       	    Daten.Geneticfitparameter.stromvektor[i] = genome.phenotype(Setfile.get_no_rconst());
       	  if (Daten.Geneticfitparameter.state_conductance[i] == 2)
       	    Daten.Geneticfitparameter.stromvektor[i] = genome.phenotype(Setfile.get_no_rconst()+1);   
         }
    }
      
//error and amplitude histogram calculation ==============
  for (i=0; i < 4096; i++) 
    {
    	ampl[i]=0;    // 2011 initialisieren  
    }
    
  ergebnis = 0;
  sigma = 0;
  ergebnis_likelihood_2d = 0;
  sigma_likelihood_2d = 0;
  ergebnis_histogram_diff = 0;
  sigma_histogram_diff = 0;
  init_nr = 0;
  
  if (Daten.Geneticfitparameter.mpi_support == true)
     {
      for (i=1; i<= Daten.Geneticfitparameter.nr_calc_mean_result; i++)
         {
       	  Ratenkonstanten[0] = i;
          if (init_nr < numprocs-1)									//send out as many genomes as procs avaiable
              {
               init_nr++;
               MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, init_nr, 1, MPI_COMM_WORLD);
              }
          else												//threreafter receive and send to same proc
              {	     
               MPI_Recv(&result, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);     
               erg[int(result[0])] = result[1];
               error_likelihood_2d[int(result[0])]=result[2];
               error_histogram_diff[int(result[0])]=result[3];
               ergebnis = ergebnis + erg[int(result[0])];
               ergebnis_likelihood_2d = ergebnis_likelihood_2d + error_likelihood_2d[int(result[0])];
               ergebnis_histogram_diff = ergebnis_histogram_diff + error_histogram_diff[int(result[0])];
               MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD); 
              }        
          }
  
      for (i=1; i <= init_nr; i++)										//collect remaining results from procs
       {
        MPI_Recv(&result, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);     
        erg[int(result[0])] = result[1];
        error_likelihood_2d[int(result[0])]=result[2];
        error_histogram_diff[int(result[0])]=result[3];
        ergebnis = ergebnis + erg[int(result[0])];
        ergebnis_likelihood_2d = ergebnis_likelihood_2d + error_likelihood_2d[int(result[0])];
        ergebnis_histogram_diff = ergebnis_histogram_diff + error_histogram_diff[int(result[0])];
       }	
      }
  else  
      for (i=1; i<= Daten.Geneticfitparameter.nr_calc_mean_result; i++)
          {     
           erg[i] = simulation(Daten.Fitparameter.multiplikator, 0, Ratenkonstanten, 0);
           error_likelihood_2d[i]=Daten.Geneticfitparameter.likelihood_2d;
           error_histogram_diff[i]=Daten.Geneticfitparameter.histogram_diff;
      
           ergebnis = ergebnis + erg[i];
           ergebnis_likelihood_2d = ergebnis_likelihood_2d + error_likelihood_2d[i];
           ergebnis_histogram_diff = ergebnis_histogram_diff + error_histogram_diff[i];
      
           for (ii=0;ii < Daten.Anz_Samples;ii++)		//2011
              ampl[Daten.Wert[ii]]++;
          } 
     
  ergebnis = ergebnis / Daten.Geneticfitparameter.nr_calc_mean_result;
  ergebnis_likelihood_2d = ergebnis_likelihood_2d / Daten.Geneticfitparameter.nr_calc_mean_result;
  ergebnis_histogram_diff = ergebnis_histogram_diff / Daten.Geneticfitparameter.nr_calc_mean_result;
  
   for (i=1; i<= Daten.Geneticfitparameter.nr_calc_mean_result; i++)
   	{
          sigma = sigma + (ergebnis- erg[i]) * (ergebnis- erg[i]);
          sigma_likelihood_2d = sigma_likelihood_2d + (ergebnis_likelihood_2d - error_likelihood_2d[i]) * (ergebnis_likelihood_2d - error_likelihood_2d[i]);
          sigma_histogram_diff = sigma_histogram_diff + (ergebnis_histogram_diff - error_histogram_diff[i]) * (ergebnis_histogram_diff - error_histogram_diff[i]);
        }   
   sigma = sigma / (Daten.Geneticfitparameter.nr_calc_mean_result-1);  
   sigma =sqrt(sigma);
   sigma_likelihood_2d = sigma_likelihood_2d / (Daten.Geneticfitparameter.nr_calc_mean_result-1);
   sigma_likelihood_2d = sqrt (sigma_likelihood_2d);
   sigma_histogram_diff = sigma_histogram_diff / (Daten.Geneticfitparameter.nr_calc_mean_result-1);
   sigma_histogram_diff = sqrt (sigma_histogram_diff);
   
//====================================================   
  std::this_thread::sleep_for(std::chrono::milliseconds(10));//wait until all slaves have posted
  std::cout<<std::endl<<std::endl<<"result:"<<std::endl;
  
if (Daten.Geneticfitparameter.fit_methode == 1)
  {
   std::cout<<"2D_Dwelltime deviation: "<<ergebnis_likelihood_2d<<"+-"<< sigma_likelihood_2d<<std::endl;  
  }
if (Daten.Geneticfitparameter.fit_methode == 2)
  {
   std::cout<<"Amplitude histogram deviation: "<<ergebnis_histogram_diff<<"+-"<<sigma_histogram_diff<<std::endl;  
  }  
if (Daten.Geneticfitparameter.fit_methode == 3)
  {
   std::cout<<"2D_Dwelltime x Amplitude histogram deviation: "<<ergebnis<<"+-"<<sigma<<std::endl;  
  }  
std::cout   <<"2D_Dwelltime deviation +-SD:    "<<ergebnis_likelihood_2d<<"   "<<sigma_likelihood_2d<<std::endl;
std::cout   <<"Amplitude histogram deviation +- SD:  "<<ergebnis_histogram_diff<<"  "<<sigma_histogram_diff<<std::endl;
     
std::cout<<std::endl<<"best kij: "<<std::endl;     

  for (i=1; i<=Setfile.get_no_rconst(); i++)
     { 
      if (Daten.Geneticfitparameter.logscale == true)	
        {
         Setfile.set_rconst(i,exp(genome.phenotype(i-1)));	
         std::cout   <<exp(genome.phenotype(i-1))<<"   ";
        }	
      else
        {  
         Setfile.set_rconst(i,exp(genome.phenotype(i-1)));	
         std::cout   << genome.phenotype(i-1)<<"   ";
        } 
     }
   if (Daten.Geneticfitparameter.levelfit == true)			// Tobias 2013 Levelfit
	     {
         Daten.Geneticfitparameter.level1 = genome.phenotype(Setfile.get_no_rconst());
         if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
       	 	Daten.Geneticfitparameter.level2 = genome.phenotype(Setfile.get_no_rconst()+1);
       	  std::cout<<std::endl; 
    	 std::cout<<"best lvl1: "<<Daten.Geneticfitparameter.level1<<" ";
    	 if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
    	 	std::cout<<"best lvl2: "<<Daten.Geneticfitparameter.level2<<" ";
    	}      
   std::cout<<std::endl<<std::endl;

//result out ==============================================================================================================   
   std::ofstream result_out(Daten.SimSettings.result_file_name,std::ios::app);
// system time
   time_t timestamp;
   tm *now;
   timestamp = time(0);
   now = localtime(&timestamp);
   result_out<<now->tm_year+1900<<"/";     
   if (now->tm_mon < 10)
   	 result_out<<"0";
   result_out<<now->tm_mon+1<<"/";
   if (now->tm_mday < 10)
   	 result_out<<"0";
   result_out<<now->tm_mday<<" ";
   if (now->tm_hour < 10)
   	   result_out<<"0";   
   result_out<<now->tm_hour<<":";
   if (now->tm_min < 10)
   	  result_out<<"0";   
   result_out<<now->tm_min<<":";
   if (now->tm_sec < 10)
   	 result_out<<"0";   
   result_out<<now->tm_sec<<" ";

//MPI size
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   result_out<<"MPI_size "<<numprocs<<" ";
//time elapsed
   //double seconds;
   //seconds = 0; //Tobias 2021 g_timer_elapsed (timer, micros);
   
   auto timer_end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed = timer_end - timer_start;
   result_out<<"sim_time "<<elapsed.count()<<" ";
// simulations
   result_out<<"nr_of_simulations "<<Setfile.iterations<<" ";   
//total fit error
   result_out<<std::setprecision(8)<<std::setw(9)<<"total_error "<<ergebnis<<" "<<sigma<<" ";
   result_out<<std::setprecision(6);	


//kij
   result_out<<"kij ";
   for (i=1; i<=Setfile.get_no_rconst(); i++)
     { 
      if (Daten.Geneticfitparameter.logscale == true)	
        {
         Setfile.set_rconst(i,exp(genome.phenotype(i-1)));	
         result_out<<exp(genome.phenotype(i-1))<<" ";	
        }	
      else
        {  
         Setfile.set_rconst(i,exp(genome.phenotype(i-1)));	
         result_out<< genome.phenotype(i-1)<<" ";
        }
     }    
  if (Daten.Geneticfitparameter.levelfit == true)			// Tobias 2013 Levelfit
	   {    
		  result_out<<"lvl ";  
      Daten.Geneticfitparameter.level1 = genome.phenotype(Setfile.get_no_rconst());
    	result_out<<Daten.Geneticfitparameter.level1<<" ";
    	if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
    	result_out<<Daten.Geneticfitparameter.level2<<" ";
     }  
//setfile
   result_out<<"setfile "<<Daten.SimSettings.loadname<<" ";
//time series
   result_out<<"time_series "<<Daten.FName<<" ";
//noisefile
   result_out<<"noisefile "<<Daten.Geneticfitparameter.file_noise_series<<" ";    //Tobias 2022
//nr error calculation
   result_out<<"nr_mean_calc.: "<<Daten.Geneticfitparameter.nr_calc_mean_result<<" ";
//2D error
result_out<<std::setprecision(8);
   result_out<<" 2D_error "<<ergebnis_likelihood_2d<<" "<<sigma_likelihood_2d<<" ";
//amplitude error
   result_out<<" ampl_error "<<ergebnis_histogram_diff<<" "<<sigma_histogram_diff<<" ";	
result_out<<std::setprecision(5);   
//times simulating
	 result_out<<"simulating_times "<<Daten.SimSettings.repeats<<" ";   //uncertain for fit
//sampling frequency
   result_out<<"sampling "<<Daten.SimSettings.sampling_frequency<<" ";  
//filter frequency
   result_out<<"filter "<<Daten.SimSettings.filter_frequency<<" "; 
//simulating multiplikator
   result_out<<"min_max_multi "<<Daten.Fitparameter.minmultiplikator<<" "<<Daten.Fitparameter.maxmultiplikator<<" ";
//action potentials?
   result_out<<"APs "<<Daten.SimSettings.actionpotentials<<" ";
//sweeps
   result_out<<"sweeps "<<Daten.SimSettings.markercounts<<" ";	
//set first n point to 0
   result_out<<"set_first "<<Daten.loeschanzahl<<" to "<<Daten.loeschvalue<<" ";
//error calculation
   result_out<<"error_calculation "<<Daten.calculate_error<<" ";
//bins per log
   result_out<<"Bpl "<<Daten.Fitparameter.bins_per_log<<" ";
//close
   result_out<<"close "<<Daten.Fitparameter.log_min_close<<" "<<Daten.Fitparameter.log_max_close<<" ";
//open
   result_out<<"open "<<Daten.Fitparameter.log_min_open<<" "<<Daten.Fitparameter.log_max_open<<" "; 
//setting 2D & ampl fit
   result_out<<"2D&Ampl "<<Daten.Geneticfitparameter.fit_methode<<" ";
//proportional fit
   result_out<<"proportional-Fit "<<Daten.Geneticfitparameter.proportionalfit<<" ";
//population size
   result_out<<"pop_size "<<Daten.Geneticfitparameter.popsize<<" ";
//generations for convergence
   result_out<<"gen_convergence "<<Daten.Geneticfitparameter.nConvergence<<" ";
//elitism
   result_out<<"elitism "<<Daten.Geneticfitparameter.elitism<<" ";
//maximum generation
   result_out<<"max_gen "<<Daten.Geneticfitparameter.ngen<<" ";
//mutation rate
   result_out<<"mutation_r "<<Daten.Geneticfitparameter.pmut<<" ";
//crossover rate
   result_out<<"crossover_r "<<Daten.Geneticfitparameter.pcross<<" ";
//bits/rate constant
  result_out<<"bits/kij "<<Daten.Geneticfitparameter.bits<<" ";
//log scaling kij
  result_out<<"log_scale_kij "<<Daten.Geneticfitparameter.logscale<<" ";
//min max settings kij
  result_out<<"kij_min_max ";
  for (i=1; i<=Setfile.get_no_rconst(); i++)
   	  result_out<<i<<" "<<Daten.Geneticfitparameter.startwertmin[i]<<" "<<Daten.Geneticfitparameter.startwertmax[i]<<" ";

   result_out<<std::endl;
   result_out.close();
   	 
   out_time();	
   // Setfile.write_setfile(Daten.SimSettings.loadname);		//new setfiles are bogus
   
   delete[] erg; 
   delete[] error_likelihood_2d;
   delete[] error_histogram_diff;
   delete[] Ratenkonstanten;
   
}
//*************************************************************************************************************************************** 
void multiprocessor_fit()
{
 MPI_Comm_rank(MPI_COMM_WORLD,&myID);
 MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
 MPI_Get_processor_name( processor_name, &namelen );

 //std::cout<<"numprocs"<<numprocs<<" myID:"<<myID<<std::endl;
 if (myID == 0)
        { 	          
         if (Daten.SaveTimeseries.automated >= 1  and Daten.SaveTimeseries.automated <= 4)
         	  multiprocessor_master_save_timeseries(); // Tobias 2020
         else	  
         		multiprocessor_master();
         
         MPI_Finalize();
        } 
 else       
        {
         if (Daten.SaveTimeseries.automated >= 1  and Daten.SaveTimeseries.automated <= 4)
         	  multiprocessor_slave_save_timeseries(); // Tobias 2020
         else	  	  	
       	    multiprocessor_slave();
       	    
       	 MPI_Finalize();
        }    
  delete[] Daten.stepresponse_loaded;			//Tobias 2022 delete storage for stepresponses                
}
//*************************************************************************************************************************************** 
void
multiprocessor_master()
{
	
 int i;
 double *Ratenkonstanten;
 int size;
//std::this_thread::sleep_for(std::chrono::milliseconds(2000));			//Tobias test 2022
MPI_Comm_size(MPI_COMM_WORLD,&size);  
strcpy(Daten.SimSettings.loadname, Daten.Geneticfitparameter.file_setfile); 
load_setfile(Daten.SimSettings.loadname); //Tobias 2013
 std::this_thread::sleep_for(std::chrono::milliseconds(size % 40));		// Tobias test 2022 std::this_thread::sleep_for(std::chrono::milliseconds(size*10+50));
 std::cout<<std::endl;
 std::cout<<"myID: "<<myID<<" on "<<processor_name<<"of  size "<<size<<" with setfile "<<Daten.SimSettings.loadname<<" master ready"<<std::endl;
 std::cout<<std::endl;
 std::this_thread::sleep_for(std::chrono::milliseconds(50));

init_2Dfit();								//viel ist hier wohl un�tig f�r master
//std::this_thread::sleep_for(std::chrono::milliseconds(500));	//Tobias test 2022

out_2D_settings();
std::this_thread::sleep_for(std::chrono::milliseconds(50));		//Tobias test 2022 std::this_thread::sleep_for(std::chrono::milliseconds(500));

 //Setzen von BIN B (simulierte Daten)														//kann weg? //viel ist hier wohl un�tig f�r master
   Daten.Dwell_2d_ptr = &Daten.Dwell_2d_B; 
   Daten.Dwell_1d_ptr = &Daten.Dwell_1d_B;  
   Daten.Dwell_2d_ptr->init((char *)"", Daten.Fitparameter.log_min_close,
  				         Daten.Fitparameter.log_max_close, 
  				         Daten.Fitparameter.log_min_open, 
  				         Daten.Fitparameter.log_max_open, 
  				         Daten.Fitparameter.bins_per_log);


 GABin2DecPhenotype map;		//genome coding
 for(i=1; i<=Setfile.get_no_rconst(); i++)
    { 
     if (Daten.Geneticfitparameter.logscale == true)
  	   map.add(Daten.Geneticfitparameter.bits, log(Daten.Geneticfitparameter.startwertmin[i]), log(Daten.Geneticfitparameter.startwertmax[i]));
     else	
        map.add(Daten.Geneticfitparameter.bits, Daten.Geneticfitparameter.startwertmin[i], Daten.Geneticfitparameter.startwertmax[i]);
    }
    
// 2012 Tobias leveldetection here
 if (Daten.Geneticfitparameter.levelfit == true) 
  	{
  	 map.add(Daten.Geneticfitparameter.bits,Daten.Geneticfitparameter.levelfit_min[1],Daten.Geneticfitparameter.levelfit_max[1]);  //1. open level
  	 if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
  	 	map.add(Daten.Geneticfitparameter.bits,Daten.Geneticfitparameter.levelfit_min[2],Daten.Geneticfitparameter.levelfit_max[2]);  //2. open level
  	}	

//Create the template genome using the phenotype map 
 GABin2DecGenome genome(map);
 GAPopulation pop(genome, Daten.Geneticfitparameter.popsize);
 pop.initializer(PopulationInitializer);
 pop.evaluator(PopulationEvaluator);

// Now create the GA  
  GASimpleGA ga(pop);

  if (Daten.Geneticfitparameter.fit_methode == 1)		//2D fit
  	{
     if (Daten.calculate_error == 0)
  	   ga.maximize();			//LogLikelihood
     else	
       ga.minimize();			//Chi square
    }        
  if (Daten.Geneticfitparameter.fit_methode == 2)		//amplitude histogram fit
  	ga.minimize();				
  if (Daten.Geneticfitparameter.fit_methode == 3){		//2022 LLH is devided by amplitude error while CHI2 is multiplied   
   if (Daten.calculate_error == 0)
  	   ga.maximize();			//LogLikelihood         //Efthymios 2024 fixed to maximize objecitve function when using combi of  LLH and amphist 
     else	
       ga.minimize();			//Chi square   
    }    
  
  
  GANoScaling scaling;
  ga.scaling(scaling);
  
  GATournamentSelector selector;
  ga.selector(selector);
  
  genome.mutator (GA1DBinaryStringGenome::FlipMutator);
  genome.crossover(GA1DBinaryStringGenome::EvenOddCrossover);

  ga.nGenerations(Daten.Geneticfitparameter.ngen);
  ga.pMutation(Daten.Geneticfitparameter.pmut);
  ga.pCrossover(Daten.Geneticfitparameter.pcross);

  ga.nConvergence(Daten.Geneticfitparameter.nConvergence);
  ga.initialize();
  if (Daten.Geneticfitparameter.elitism)			//Tobias 2014, new elitism setting
  	ga.elitist(gaTrue); 		
  else	
  	ga.elitist(gaFalse);
  	
  generation_out_2Dfit(ga, genome);
  
 double convergence = 0; 
 while( ((convergence == 0) || (convergence>1)) && (ga.statistics().generation() < Daten.Geneticfitparameter.ngen) )				//( ((ga.statistics().convergence() == 0) || (ga.statistics().convergence()>1)) && (ga.statistics().generation() < Daten.Geneticfitparameter.ngen) )
 { 		 
  ga.step();
  convergence = ga.statistics().convergence();
  if (Daten.calculate_error == 0)
  	if (convergence != 0)
  		convergence = 1/convergence;
   generation_out_2Dfit(ga, genome);	
 } 
 
  
 Ratenkonstanten = new double[Setfile.get_no_rconst()+3];		// 1st termination signal for nodes
 Ratenkonstanten[0] = -1;
 for (i=1; i < numprocs; i++)
     MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
 
 summary_out_2Dfit(ga, genome);						
 
 for (i=1; i < numprocs; i++)																  // 2nd termination signal for nodes
     MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
 
 delete[] Ratenkonstanten; 
}

//*************************************************************************************************************************************** 
void
multiprocessor_slave()
{	
 int ii;
 MPI_Status status;
 double *Ratenkonstanten;		
 double ergebnis[4];

 strcpy(Daten.SimSettings.loadname, Daten.Geneticfitparameter.file_setfile); 
 load_setfile(Daten.SimSettings.loadname); //Tobias 2013
 std::this_thread::sleep_for(std::chrono::milliseconds(myID % 40));		//Tobias test 2022 std::this_thread::sleep_for(std::chrono::milliseconds(myID*10));
 std::cout<<"myID: "<<myID<<" on "<<processor_name<<" setfile "<<Daten.SimSettings.loadname<<" timeseries: "<<Daten.FName<<std::endl;
 	
 init_2Dfit();
 Ratenkonstanten = new double[Setfile.get_no_rconst()+3]; 
 //Setzen von BIN B (simulierte Daten)
   Daten.Dwell_2d_ptr = &Daten.Dwell_2d_B; 
   Daten.Dwell_1d_ptr = &Daten.Dwell_1d_B;  
   Daten.Dwell_2d_ptr->init((char *)"", Daten.Fitparameter.log_min_close,
  				         Daten.Fitparameter.log_max_close, 
  				         Daten.Fitparameter.log_min_open, 
  				         Daten.Fitparameter.log_max_open, 
  				         Daten.Fitparameter.bins_per_log); 
 while (1)
   { 
//std::this_thread::sleep_for(std::chrono::milliseconds(1500));  	//Tobias test 2022


	
    MPI_Recv(Ratenkonstanten,        			/* message buffer */
             Setfile.get_no_rconst()+3,			    /* one data item */
             MPI_DOUBLE,        			         /* of type double real */
             0,						                          /* receive from any sender */
             1,       					                     /* any type of message */
             MPI_COMM_WORLD,    			/* default communicator */
             &status);          			                  /* info about the received message */	  
                      
    if (Ratenkonstanten[0] < 0)
       break;                                                
		 
   if (Daten.Geneticfitparameter.proportionalfit == true)  
       {
   	for (ii=1; ii <= Setfile.get_no_rconst(); ii++)    			//Neu ratio
     		{
       		 Setfile.set_rconst(ii,Ratenkonstanten[ii]);   
       		 Setfile.set_rconst(ii+1,Ratenkonstanten[ii]/100*Ratenkonstanten[ii+1]);
       		 ii++;
     		} 
      }	
     		
   else
      {
     	for (ii=1; ii <= Setfile.get_no_rconst(); ii++)	
     	   Setfile.set_rconst(ii,Ratenkonstanten[ii]); 	
      }
         
   if (Daten.Geneticfitparameter.levelfit == true)  //Tobias 2013 Levelfit
      {
       Daten.Geneticfitparameter.level1 = Ratenkonstanten[Setfile.get_no_rconst()+1];
       if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
       		Daten.Geneticfitparameter.level2 = Ratenkonstanten[Setfile.get_no_rconst()+2]; 	
       for (ii=1; ii <= 10; ii++)
         {
       	  if (Daten.Geneticfitparameter.state_conductance[ii] == 0)
       	     Daten.Geneticfitparameter.stromvektor[ii] = 0;
       	  if (Daten.Geneticfitparameter.state_conductance[ii] == 1)
       	     Daten.Geneticfitparameter.stromvektor[ii] = Daten.Geneticfitparameter.level1;
       	  if (Daten.Geneticfitparameter.state_conductance[ii] == 2)
       	     Daten.Geneticfitparameter.stromvektor[ii] = Daten.Geneticfitparameter.level2;   	   
       }      	  
     }
     
     if (Daten.SaveTimeseries.random_noise)			//Tobias 2021 randomize noise
    		{
    	 	 Daten.SimSettings.resample = true;	
    	 	 Setfile.set_sigma (Daten.SaveTimeseries.min_noise + Daten.uni_01(Daten.random_engine)*(Daten.SaveTimeseries.max_noise - Daten.SaveTimeseries.min_noise));
    		}  
     
      
    if (Daten.SimSettings.dosbox) 
    	std::cout<<"myID:"<<std::setw(5)<<myID<<" on "<<std::setw(8)<<processor_name<<" ";
    		
    		
    ergebnis[1] = simulation(Daten.Fitparameter.multiplikator, 0, Ratenkonstanten, 0);    
    ergebnis[0] = Ratenkonstanten[0];							//receive & transfer genome number via Ratenkonstanten[0]
    
    MPI_Send(&ergebnis,             			/* message buffer */
                 2,             			/* nr of data items */
                 MPI_DOUBLE,           			/* data item is an integer */
             	 0,              			/* destination process rank */
                 1,           				/* user chosen message tag */
                 MPI_COMM_WORLD);            
	 }
   
  while (1)										              //calculate results
     {    	
      MPI_Recv(Ratenkonstanten,        			/* message buffer */
             Setfile.get_no_rconst()+3,			         /* one data item */
             MPI_DOUBLE,        			              /* of type double real */
             0,						                               /* receive from any sender */
             1,       					                          /* any type of message */
             MPI_COMM_WORLD,    			    /* default communicator */
             &status);          			                       /* info about the received message */	 
     if (Ratenkonstanten[0] < 0)
        break;   
  
   
   if (Daten.Geneticfitparameter.proportionalfit == true)
      {  
   	for (ii=1; ii <= Setfile.get_no_rconst(); ii++)    			//Neu ratio
     		{
       		 Setfile.set_rconst(ii,Ratenkonstanten[ii]);   
       		 Setfile.set_rconst(ii+1,Ratenkonstanten[ii]/100*Ratenkonstanten[ii+1]);
       		 ii++;
     		} 
      }		
   else
      {
     	for (ii=1; ii <= Setfile.get_no_rconst(); ii++)	
     	   Setfile.set_rconst(ii,Ratenkonstanten[ii]); 	
      } 
       
         
     if (Daten.Geneticfitparameter.levelfit == true)  //Tobias 2013 Levelfit
        {
         Daten.Geneticfitparameter.level1 = Ratenkonstanten[Setfile.get_no_rconst()+1];
         if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
         	  Daten.Geneticfitparameter.level2 = Ratenkonstanten[Setfile.get_no_rconst()+2]; 	
         for (ii=1; ii <= 10; ii++)
            {
       	     if (Daten.Geneticfitparameter.state_conductance[ii] == 0)
       	        Daten.Geneticfitparameter.stromvektor[ii] = 0;
       	     if (Daten.Geneticfitparameter.state_conductance[ii] == 1)
       	        Daten.Geneticfitparameter.stromvektor[ii] = Daten.Geneticfitparameter.level1;
       	     if (Daten.Geneticfitparameter.state_conductance[ii] == 2)
       	        Daten.Geneticfitparameter.stromvektor[ii] = Daten.Geneticfitparameter.level2;   	   
            }
        }     
      if (Daten.SimSettings.dosbox) 
    		std::cout<<"myID:"<<std::setw(5)<<myID<<" on "<<std::setw(8)<<processor_name<<" ";
    		
      ergebnis[1] = simulation(Daten.Fitparameter.multiplikator, 0, Ratenkonstanten, 0);
      ergebnis[2] = Daten.Geneticfitparameter.likelihood_2d;
      ergebnis[3] = Daten.Geneticfitparameter.histogram_diff;
      ergebnis[0] = Ratenkonstanten[0];							//receive & transfer genome number via Ratenkonstanten[0]
      
      MPI_Send(&ergebnis,             			/* message buffer */
                 4,             								/* nr of data items */
                 MPI_DOUBLE,           			/* data item is an integer */
             	   0,              								/* destination process rank */
                 1,           								/* user chosen message tag */
                 MPI_COMM_WORLD);               
     }
   
   
 std::cout<<Setfile.iterations<<" simulations were performed by "<<myID<<" on "<<processor_name<<std::endl;
 delete[] Ratenkonstanten;
 //std::this_thread::sleep_for(std::chrono::milliseconds(100));																	//Tobias test 2022
}

//*************************************************************************************************************************************** 

void PopulationInitializer(GAPopulation & pop)
{
 int i;
 
 for (i=0; i < pop.size();i++)
   pop.individual(i).initialize();	
 std::cout<<"population initialized"<<std::endl;	
}
//*************************************************************************************************************************************** 


//important-> reevaluation of every gene in contrast to non-multiprocessor-fit
void PopulationEvaluator(GAPopulation & pop)
{
 int i,ii;
 int init_nr;
 double *Ratenkonstanten;
 Ratenkonstanten = new double[Setfile.get_no_rconst()+3];
 Ratenkonstanten[0] = -1;
 double ergebnis[2];
 MPI_Status status;		
  		

 GABin2DecPhenotype map;		//genome coding																									//kan man sich die map sparen da in pop?
 for(i=1; i<=Setfile.get_no_rconst();i++)
     if (Daten.Geneticfitparameter.logscale == true)
  	map.add(Daten.Geneticfitparameter.bits, log(Daten.Geneticfitparameter.startwertmin[i]), log(Daten.Geneticfitparameter.startwertmax[i]));
     else	
        map.add(Daten.Geneticfitparameter.bits, Daten.Geneticfitparameter.startwertmin[i], Daten.Geneticfitparameter.startwertmax[i]);
// 2012 Tobias leveldetection here
 if (Daten.Geneticfitparameter.levelfit == true) 
  	{
  	 map.add(Daten.Geneticfitparameter.bits,Daten.Geneticfitparameter.levelfit_min[1],Daten.Geneticfitparameter.levelfit_max[1]);  //1. open level
  	 if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
  	 	map.add(Daten.Geneticfitparameter.bits,Daten.Geneticfitparameter.levelfit_min[2],Daten.Geneticfitparameter.levelfit_max[2]);  //2. open level
  	}	
 GABin2DecGenome genome(map);

auto timer_start = std::chrono::system_clock::now(); //Tobias 2022 
 init_nr = 0;	
 for (ii=0; ii < pop.size(); ii++)
      {   		
      	
      	Setfile.iterations++;
      	Ratenkonstanten[0] = ii;
        genome = pop.individual(ii);
         
        for (i=1; i<=Setfile.get_no_rconst(); i++)
           {
            if (Daten.Geneticfitparameter.logscale == true)
              Ratenkonstanten[i]=exp(genome.phenotype(i-1)); 
            else	
              Ratenkonstanten[i]=genome.phenotype(i-1);
           }  
        if (Daten.Geneticfitparameter.levelfit == true)  //Tobias 2013 Levelfit
               {
       	        Ratenkonstanten[Setfile.get_no_rconst()+1] = genome.phenotype(Setfile.get_no_rconst());
       	        if(Daten.Geneticfitparameter.levelfit_max[2] != 0)
       	        	Ratenkonstanten[Setfile.get_no_rconst()+2] = genome.phenotype(Setfile.get_no_rconst()+1);  
               }           

        if (init_nr < numprocs-1)										//send out as many genomes as procs available
             {
              init_nr++;
              MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, init_nr, 1, MPI_COMM_WORLD);
             }
        else												//threreafter receive and send to same proc
             {	     
              MPI_Recv(&ergebnis, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);   
              pop.individual(int(ergebnis[0])).score(ergebnis[1]);     	
              MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD); 
             }        
      }
  
  for (i=1; i <= init_nr; i++)										//collect remaining results from procs
       {
        MPI_Recv(&ergebnis, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);          
        pop.individual(int(ergebnis[0])).score(ergebnis[1]);   
       }

     
          
 delete[] Ratenkonstanten;
}

//#####################################################################################################################################################

void multiprocessor_master_save_timeseries() // Tobias 2020
{
 long int i,ii;
 long int init_nr;
 long int header_size, ts_size, dth_size, total_size;
 double *Ratenkonstanten;
 double ergebnis[2];
 int size;
 MPI_Status status;	
 MPI_Comm_size(MPI_COMM_WORLD,&size);  
 std::string save_name,dummy;
 

 
 strcpy(Daten.SimSettings.loadname, Daten.Geneticfitparameter.file_setfile); 
 load_setfile(Daten.SimSettings.loadname); //Tobias 2013


 std::this_thread::sleep_for(std::chrono::milliseconds(50+size*10));
 std::cout<<std::endl;
 std::cout<<"myID: "<<myID<<" on "<<processor_name<<"of  size "<<size<<" with setfile "<<Daten.SimSettings.loadname<<" master ready"<<std::endl;
 std::cout<<std::endl;
 std::this_thread::sleep_for(std::chrono::milliseconds(50));
  	
 init_2Dfit();			
 Daten.Anz_Samples = Setfile.get_no_samples();		//Tobias 2022					
 std::this_thread::sleep_for(std::chrono::milliseconds(500));
 Ratenkonstanten = new double[Setfile.get_no_rconst()+3];		

 init_nr = 0;
 header_size = 1000;				//const if changed in slave must be here as well
 ts_size = sizeof(Daten.Wert[0])*Daten.Anz_Samples;
 dth_size = sizeof(Daten.Dwell_2d_B.CO[0][0])*
   					(Daten.Fitparameter.log_max_close-Daten.Fitparameter.log_min_close)*Daten.Fitparameter.bins_per_log*
   					(Daten.Fitparameter.log_max_open-Daten.Fitparameter.log_min_open)*Daten.Fitparameter.bins_per_log;		
   					
 save_name = "";	 																																//prepare name for save
 //save_name.append(Daten.SaveTimeseries.folder);  																	// add path
 dummy = Daten.Geneticfitparameter.file_setfile;
 save_name.append(dummy.substr(0,dummy.find(".")));	   					
 save_name.append("_samples_");
 save_name.append(std::to_string(Daten.SaveTimeseries.number));  					

 std::cout<<"++++++++++++ saving "<<Daten.SaveTimeseries.number<<" timeseries with nr of "<<Daten.Anz_Samples<<" samples to "<<save_name;
 if (Daten.SaveTimeseries.automated == 1)
 	  {
 	   std::cout<<" timeseries"<<std::endl;
 	   total_size = header_size + ts_size;
 	   save_name.append("_header_");
 	   save_name.append(std::to_string(header_size)); 
 	   save_name.append("_timeseries_");
 	   save_name.append(std::to_string(ts_size));
 	   save_name.append(".dat");	 	   	
 	  } 	 	
 else if (Daten.SaveTimeseries.automated == 2)
 	  {
 	   std::cout<<" 2D-histograms"<<std::endl;
 	   total_size = header_size + dth_size;
 	   save_name.append("_header_");
 	   save_name.append(std::to_string(header_size)); 
 	   save_name.append("_2DtH_");
 	   save_name.append(std::to_string(dth_size));
 	   save_name.append(".dat");	 			
 	  } 	 
 else if (Daten.SaveTimeseries.automated == 3)
 	  {
 	   std::cout<<" timeseries and 2D histograms"<<std::endl;
 	   total_size = header_size + ts_size + dth_size;
 	   save_name.append("_header_");
 	   save_name.append(std::to_string(header_size));
 	   save_name.append("_timeseries_");
 	   save_name.append(std::to_string(ts_size)); 
 	   save_name.append("_2DtH_");
 	   save_name.append(std::to_string(dth_size));
 	   save_name.append(".dat");	 				 	
 	  }
 else if (Daten.SaveTimeseries.automated == 4)
 	  {
 	   std::cout<<" 2D-histograms"<<std::endl;
 	   total_size = header_size + 9*dth_size;
 	   save_name.append("_header_");
 	   save_name.append(std::to_string(header_size)); 
 	   save_name.append("_2DtH_");
 	   save_name.append(std::to_string(9*dth_size));
 	   save_name.append(".dat");	 			
 	  } 	 
 else 
 	  {
   	 std::cout<<"option for dataset generation not supported"<<std::endl; 	
   	 return;	
   	}	  
  
   		   
 std::cout<<"MPI transfer buffer size: "<<total_size<<std::endl;
 unsigned char* buffer = new unsigned char[total_size];  		  	 	
 std::cout<<std::endl;
 std::cout<<std::endl;	 		 	
 
 char* char_save_name;
 char_save_name = &save_name[0];
 FILE* stream = std::fopen(char_save_name, "wb");					//open file for binary output	
 if (stream == NULL)
 	 {
 	 	std::cout<<"Failed to open save file "<<char_save_name<<std::endl;
 	 	std::this_thread::sleep_for(std::chrono::milliseconds(5000)); 		
    Ratenkonstanten[0] = -1;
    for (i=1; i < numprocs; i++)
      MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);			//termination
    MPI_Finalize();	 	
 	 }				
 
 unsigned long int received = 1;
  
 for (ii=0; ii < Daten.SaveTimeseries.number; ii++)		//Tobias 2020 manual edit
      {   		
      	
      	Ratenkonstanten[0] = ii;
         
        for (i=1; i<=Setfile.get_no_rconst(); i++)
           {
            Ratenkonstanten[i] = Setfile.get_rconst(i);
					 }           

        if (init_nr < numprocs-1)										//send out as many genomes as procs available 
             {
              init_nr++;
              MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, init_nr, 1, MPI_COMM_WORLD);
             }
        else												//threreafter receive and send to same proc
             {	     
              MPI_Recv(buffer, total_size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
              if (std::fwrite(buffer, 1, total_size, stream) < total_size)
              	std::cout<<"##################error writing into save file"<<std::endl; 
              
        		  int progress = received % (Daten.SaveTimeseries.number / 100+1);  	
              if (progress == 0)		               
                std::cout<<"data saved: "<<received<<std::endl; 		
              received++;   
              //std::cout<<"**********************"<<received<<std::endl; 	
              MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD); 
              //std::cout<<"**********************"<<received<<std::endl; 	
             }        
      }
   
  for (i=1; i <= init_nr; i++)										//collect remaining results from procs
       {
        MPI_Recv(buffer, total_size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);               
        if (std::fwrite(buffer, 1, total_size, stream) < total_size)
          std::cout<<"##################error writing into save file"<<std::endl;
        int progress = received % (Daten.SaveTimeseries.number / 100+1);  	
        if (progress == 0)		               
          std::cout<<"data saved: "<<received<<std::endl;  	
        received++;	     
       }     
 
 
 Ratenkonstanten[0] = -1;
 for (i=1; i < numprocs; i++)
     MPI_Send(Ratenkonstanten, Setfile.get_no_rconst()+3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);			//termination
 	
 std::fclose(stream);							//close save file
 delete[] Ratenkonstanten;
 delete[] buffer;
 std::cout<<"total of "<<received - 1<<" data sets were created"<<std::endl;
 std::cout<<"time needed for completion: ";
 out_time();	 
 std::this_thread::sleep_for(std::chrono::milliseconds(100)); 	
}

//*************************************************************************************************************************************** 


void multiprocessor_slave_save_timeseries()		// Tobias 2020
{
 long int ii;
 long unsigned int laenge;
 std::string save_name, save_name_ts, save_name_2D, save_name_ts2D, dummy;
 MPI_Status status;
 double *Ratenkonstanten;		
 double ergebnis[4];
 double helper;
 long unsigned int header_size, ts_size, dth_count, dth_size, total_size;
double random_SNR, random_sigma;
double max_level;


 strcpy(Daten.SimSettings.loadname, Daten.Geneticfitparameter.file_setfile); 
 load_setfile(Daten.SimSettings.loadname); //Tobias 2013
 init_2Dfit();
 Daten.Anz_Samples = Setfile.get_no_samples();		//Tobias 2022	
 Ratenkonstanten = new double[Setfile.get_no_rconst()+3];	
    
 Daten.Dwell_2d_ptr = &Daten.Dwell_2d_B; 
 Daten.Dwell_1d_ptr = &Daten.Dwell_1d_B;  
 Daten.Dwell_2d_ptr->init((char *)"",
 									 Daten.Fitparameter.log_min_close,
  				         Daten.Fitparameter.log_max_close, 
  				         Daten.Fitparameter.log_min_open, 
  				         Daten.Fitparameter.log_max_open, 
  				         Daten.Fitparameter.bins_per_log);		//Tobias 2021 "0"
 
     

  while (1)
   {
    MPI_Recv(Ratenkonstanten,        			/* message buffer */
             Setfile.get_no_rconst()+3,			    /* one data item */
             MPI_DOUBLE,        			         /* of type double real */
             0,						                          /* receive from any sender */
             1,       					                     /* any type of message */
             MPI_COMM_WORLD,    			/* default communicator */
             &status);          			                  /* info about the received message */	         
    if (Ratenkonstanten[0] < 0)
       break;                                             

   if (Daten.SaveTimeseries.random_rate)									//Tobias 2020 randomize rconsts according to the boundaries in the setfile
   	 {
   	 	for (ii=1; ii <= Setfile.get_no_rconst(); ii++)
   	 	  {
   	 	   helper = std::log(Daten.Geneticfitparameter.startwertmin[ii])+Daten.uni_01(Daten.random_engine)*(std::log(Daten.Geneticfitparameter.startwertmax[ii])-std::log(Daten.Geneticfitparameter.startwertmin[ii]));
   	 	   //std::cout<<"Daten.Geneticfitparameter.startwertmin[ii]: "<<Daten.Geneticfitparameter.startwertmin[ii]<<" Daten.Geneticfitparameter.startwertmax[ii]: "<<Daten.Geneticfitparameter.startwertmax[ii]<<" helperlog: "<<helper<<" helper: "<<std::exp(helper)<<std::endl;	
   	 	   helper = std::exp(helper);
   	 	   Setfile.set_rconst(ii,helper);  	
   	 	  }	   
     }
   else
   	 {  
   	  for (ii=1; ii <= Setfile.get_no_rconst(); ii++)    			
       	  	 Setfile.set_rconst(ii,Ratenkonstanten[ii]);   
     }	    
   
   if (Daten.SaveTimeseries.random_lvl)  //Tobias 2020 randomize current level
        {
         Daten.Geneticfitparameter.level1 = Daten.Geneticfitparameter.levelfit_min[1] +Daten.uni_01(Daten.random_engine)*(Daten.Geneticfitparameter.levelfit_max[1]-Daten.Geneticfitparameter.levelfit_min[1]);
         Daten.Geneticfitparameter.level2 = Daten.Geneticfitparameter.levelfit_min[2] +Daten.uni_01(Daten.random_engine)*(Daten.Geneticfitparameter.levelfit_max[2]-Daten.Geneticfitparameter.levelfit_min[2]);
         for (ii=1; ii <= 10; ii++)
            {
       	     if (Daten.Geneticfitparameter.state_conductance[ii] == 0)
       	        Daten.Geneticfitparameter.stromvektor[ii] = 0;
       	     if (Daten.Geneticfitparameter.state_conductance[ii] == 1)
       	        Daten.Geneticfitparameter.stromvektor[ii] = Daten.Geneticfitparameter.level1;
       	     if (Daten.Geneticfitparameter.state_conductance[ii] == 2)
       	        Daten.Geneticfitparameter.stromvektor[ii] = Daten.Geneticfitparameter.level2;   	   
            }
        }   
  
    if (Daten.SaveTimeseries.random_noise)			//Tobias 2021 randomize noise
    	{
    	 Daten.SimSettings.resample = true;	
	 		 random_SNR = Daten.SaveTimeseries.min_noise + Daten.uni_01(Daten.random_engine)*(Daten.SaveTimeseries.max_noise - Daten.SaveTimeseries.min_noise);
	     if (Daten.Geneticfitparameter.level1 > Daten.Geneticfitparameter.level2)
	     		max_level = Daten.Geneticfitparameter.level1;
	     else
	     		max_level = Daten.Geneticfitparameter.level2;
	     		

	     random_sigma = max_level / random_SNR;
    	 Setfile.set_sigma (random_sigma);
std::cout<<"min SNR:"<<Daten.SaveTimeseries.min_noise<<" max SNR: "<<Daten.SaveTimeseries.max_noise<<" simulated SNR: "<<random_SNR<<
	" level1, Level2, max level: "<<Daten.Geneticfitparameter.level1<<", "<<Daten.Geneticfitparameter.level2<<", "<<max_level<<" random_sigma: "<<random_sigma<<std::endl;
    	} 		

  
    if (Daten.SimSettings.dosbox) 
    	std::cout<<"myID:"<<std::setw(5)<<myID<<" on "<<std::setw(8)<<processor_name<<" ";


      ergebnis[1] = simulation(Daten.Fitparameter.multiplikator, 0, Ratenkonstanten, 0);

 		


    save_name = "";	 

    dummy = Daten.Geneticfitparameter.file_setfile;
    char* dummy1;
    dummy1 = Setfile.get_stepresponse_filename(NULL); 

    //std::cout<<"dummy: "<<dummy<<std::endl;
    save_name.append(dummy.substr(0,dummy.find(".")));																				// add filename from setfile
    save_name.append("_kij");																																	// add kij to name
    for (ii=1; ii <= Setfile.get_no_rconst(); ii++)    			//N
     		{
     			 save_name.append("_");
       		 save_name.append(std::to_string(int(Setfile.get_rconst(ii)+0.5)));   
     		} 
   save_name.append("_s_");																							//add sigma noise to name
   save_name.append(std::to_string((int)round(Daten.a_fit.sigma)));				//Tobias 2021

	save_name.append("_lvl");	
	for (ii=1; ii <= Setfile.get_matrix_dimension(); ii++)    			//N
     		{
         save_name.append("_");																									//add lvl
         save_name.append(std::to_string(int(Daten.Geneticfitparameter.stromvektor[ii])));
			  }

//Efthymios 2023 START denote step response and noise type in header


   save_name.append("_NoiseType_");

   if (Daten.Geneticfitparameter.noise_generation == 2)																						//load recorded noise
   {			
      save_name.append("fromTS");
   }     
   
else if (Daten.Geneticfitparameter.noise_generation == 1)																				//generate noise with spectrum
   {			
         save_name.append("spectral");
   }      
   
   
else if (Daten.Geneticfitparameter.noise_generation == 0)																				//generate simple simulated noise
	{
	   save_name.append("lpWhite");
   }

save_name.append("_StepResp_");


    std::string dummy2 = "";
    for (int i = 0; i < sizeof(dummy1) - 5; i++) {  //Dateiendung vom string entfernen
        dummy2 = dummy2 + dummy1[i];
    }

save_name.append(dummy2);



   //std::cout<<"myID: "<<myID<<" on "<<processor_name<<" transferring data: "<<save_name<<std::endl; 	
    
   if (Daten.SaveTimeseries.automated == 1) 	// transfer header + time series to master
    	{

      save_name_ts = save_name;   
      save_name_ts.append("_ts");
      save_name_ts.append(".000000000");
      dummy = std::to_string(int(Ratenkonstanten[0]));
      save_name_ts.replace(save_name_ts.length()-dummy.length(), dummy.length(), dummy);
      save_name_ts.resize(1000,'+');
      ts_size = sizeof(Daten.Wert[0])*Daten.Anz_Samples;

      header_size = save_name_ts.length();
      //std::cout<<"header_size: "<<header_size<<std::endl;
      //std::cout<<"ts_size: "<<ts_size<<std::endl;

    	 total_size = header_size + ts_size;
    	 //std::cout<<"total_size: "<<total_size<<std::endl;
    	 unsigned char *buffer = new unsigned char[total_size];
    	 for (ii=0; ii < header_size; ii++)
    	    buffer[ii] = save_name_ts.at(ii);																	//transfer header into transfer buffer
    	 std::memcpy (&buffer[header_size], &Daten.Wert[0], ts_size); 				//copy time series into buffer
    	 
    	 MPI_Send(buffer,             			/* message buffer */
                 total_size,             			/* nr of data items */
                 MPI_UNSIGNED_CHAR,           			/* data item is a char */
             	   0,              			/* destination process rank */
                 1,           				/* user chosen message tag */
                 MPI_COMM_WORLD);  	
    	 delete[] buffer;  	
    	} 	

   else if (Daten.SaveTimeseries.automated == 2)  // transfer header + 2D histogram to master
    	{    		

      save_name_2D = save_name;   
      save_name_2D.append("_2D");
      save_name_2D.append(".000000000");
      dummy = std::to_string(int(Ratenkonstanten[0]));
      save_name_2D.replace(save_name_2D.length()-dummy.length(), dummy.length(), dummy);		// add number
      save_name_2D.resize(1000,'+');
      header_size = save_name_2D.length();
      dth_count = (Daten.Fitparameter.log_max_close-Daten.Fitparameter.log_min_close)*Daten.Fitparameter.bins_per_log*
   					(Daten.Fitparameter.log_max_open-Daten.Fitparameter.log_min_open)*Daten.Fitparameter.bins_per_log;

      dth_size = dth_count*sizeof (Daten.Dwell_2d_B.CO[0][0]);	



    		total_size = header_size + dth_size;
    		unsigned char *buffer = new unsigned char[total_size];
    		//std::cout<<"total_size: "<<total_size<<std::endl;
    		for (ii=0; ii < header_size; ii++)
    	     buffer[ii] = save_name_2D.at(ii);															//transfer header into trasnfer buffer

				double *temp_buffer = new double[dth_count];										//temp buffer for dwell-time data
				ii=0;
				for (int id=0;id < (Daten.Fitparameter.log_max_close-Daten.Fitparameter.log_min_close)*Daten.Fitparameter.bins_per_log;id++)			//transfer dwell-time data into temp buffer
    				 	 for (int jd=0;jd < (Daten.Fitparameter.log_max_open-Daten.Fitparameter.log_min_open)*Daten.Fitparameter.bins_per_log;jd++)
									{
	  						 	 temp_buffer[ii] = Daten.Dwell_2d_B.CO[id][jd];
	  						 	 ii++;
									}	
    		std::memcpy(&buffer[header_size], &temp_buffer[0], dth_size);	

    	  MPI_Send(buffer,             			/* message buffer */
                 total_size,             			/* nr of data items */
                 MPI_UNSIGNED_CHAR,           			/* data item is a char */
             	   0,              			/* destination process rank */
                 1,           				/* user chosen message tag */
                 MPI_COMM_WORLD);            

    		delete[] buffer; 
    		delete[] temp_buffer; 	
    	}
    			
   else if (Daten.SaveTimeseries.automated == 3) // transfer header + timeseries + 2D histogram to master


    	{

         save_name_ts2D = save_name; 
         save_name_ts2D.append("_ts2D");
         save_name_ts2D.append(".000000000");
         header_size = save_name_ts2D.length();
         ts_size = sizeof(Daten.Wert[0])*Daten.Anz_Samples;
         dth_count = (Daten.Fitparameter.log_max_close-Daten.Fitparameter.log_min_close)*Daten.Fitparameter.bins_per_log*
   					(Daten.Fitparameter.log_max_open-Daten.Fitparameter.log_min_open)*Daten.Fitparameter.bins_per_log;
         dth_size = dth_count*sizeof (Daten.Dwell_2d_B.CO[0][0]);	  

         //std::cout<<"header_size: "<<total_size<<std::endl;  
         //std::cout<<"ts_size: "<<ts_size<<std::endl;  
         //std::cout<<"dth_size: "<<dth_size<<std::endl;      
    		total_size = header_size + ts_size + dth_size;

    		unsigned char *buffer = new unsigned char[total_size];
    		for (ii=0; ii < header_size; ii++)
    	     buffer[ii] = save_name_ts2D.at(ii);
    	  std::memcpy (&buffer[header_size], &Daten.Wert[0], ts_size);
    	  
				double *temp_buffer = new double[dth_count];										//temp buffer for dwell-time data
				ii=0;
				for (int id=0;id < (Daten.Fitparameter.log_max_close-Daten.Fitparameter.log_min_close)*Daten.Fitparameter.bins_per_log;id++)			//transfer dwell-time data into temp buffer
    				 	 for (int jd=0;jd < (Daten.Fitparameter.log_max_open-Daten.Fitparameter.log_min_open)*Daten.Fitparameter.bins_per_log;jd++)
									{
	  						 	 temp_buffer[ii] = Daten.Dwell_2d_B.CO[id][jd];
	  						 	 ii++;
									}	
    		std::memcpy(&buffer[header_size + ts_size], &temp_buffer[0], dth_size);	//Efthymios 2023 changed from: std::memcpy(&buffer[header_size], &temp_buffer[0], dth_size); 
    			
    		//std::cout<<"total_size: "<<total_size<<std::endl;	
    	  MPI_Send(buffer,             			/* message buffer */
                 total_size,             			/* nr of data items */
                 MPI_UNSIGNED_CHAR,           			/* data item is a char */
             	   0,              			/* destination process rank */
                 1,           				/* user chosen message tag */
                 MPI_COMM_WORLD); 	
    		delete[] buffer;  	
    		delete[] temp_buffer; 	 	
    	}
  
  
  
   else if (Daten.SaveTimeseries.automated == 4)  // transfer header + 2D multi level histogram to master
    	{    		

      save_name_2D = save_name;   
      save_name_2D.append("_2Dmult");
      save_name_2D.append(".000000000");
      dummy = std::to_string(int(Ratenkonstanten[0]));
      save_name_2D.replace(save_name_2D.length()-dummy.length(), dummy.length(), dummy);		// add number
      save_name_2D.resize(1000,'+');
      header_size = save_name_2D.length();
      dth_count = (Daten.Fitparameter.log_max_close-Daten.Fitparameter.log_min_close)*Daten.Fitparameter.bins_per_log*
         (Daten.Fitparameter.log_max_open-Daten.Fitparameter.log_min_open)*Daten.Fitparameter.bins_per_log;
            
      dth_size = 9*dth_count*sizeof (Daten.Dwell_2d_B.CO[0][0]);
      double levels[9] = {-0.4, -0.3, -0.2, -0.1, -0.0, 0.1, 0.2, 0.3, 0.4};

         
         double original_base = Daten.a_fit.i_null;
         double original_amp = Daten.a_fit.i_channel;
         total_size = header_size + dth_size;
    		unsigned char *buffer = new unsigned char[total_size];

         double *temp_buffer = new double[9*dth_count];		



    		//std::cout<<"total_size: "<<total_size<<std::endl;
    		for (ii=0; ii < header_size; ii++)
    	     buffer[ii] = save_name_2D.at(ii);															//transfer header into trasnfer buffer
							//temp buffer for dwell-time data


               int i;
             for (i=0;i < 9;i++)	//9 levels
         {
           
 
  	      Daten.a_fit.i_null =  original_base + original_amp*levels[i];
		   Daten.a_fit.i_channel = original_amp - original_amp*levels[i];    
     	   Daten.hinkley(_HOHD);		
     	   Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);
         //std::cout<<"level: "<<Daten.a_fit.i_channel<<" nr_of events: "<<Daten.Dwell_2d_ptr->events()<<std::endl;
            //std::cout<<"level: "<<Daten.a_fit.i_channel<<" nr_of events: "<<Daten.Dwell_2d_ptr->events()<<std::endl;
            									//temp buffer for dwell-time data
				ii=0;
				for (int id=0;id < (Daten.Fitparameter.log_max_close-Daten.Fitparameter.log_min_close)*Daten.Fitparameter.bins_per_log;id++)			//transfer dwell-time data into temp buffer
    				 	 for (int jd=0;jd < (Daten.Fitparameter.log_max_open-Daten.Fitparameter.log_min_open)*Daten.Fitparameter.bins_per_log;jd++)
									{
	  						 	 temp_buffer[i*dth_count + ii] = Daten.Dwell_2d_B.CO[id][jd];

                         //if(i == 2){
                          // if(temp_buffer[i*dth_count + ii]>0){
                         //std::cout<<"dummy: "<<ii<<" : "<<temp_buffer[i*dth_count + ii]<<std::endl;
                         //std::cout<<"dummy2: "<<ii<<" : "<<Daten.Dwell_2d_B.CO[id][jd]<<std::endl;
                          // }
                         //}
                         ii++;
									}
         }	

    		std::memcpy(&buffer[header_size], &temp_buffer[0], dth_size);	

    	  MPI_Send(buffer,             			/* message buffer */
                 total_size,             			/* nr of data items */
                 MPI_UNSIGNED_CHAR,           			/* data item is a char */
             	   0,              			/* destination process rank */
                 1,           				/* user chosen message tag */
                 MPI_COMM_WORLD);

                 //std::cout<<"buffer: "<<buffer<<std::endl;	
                 //std::cout<<"total_size: "<<total_size<<std::endl;	
                 //std::cout<<"MPI_UNSIGNED_CHAR: "<<MPI_UNSIGNED_CHAR<<std::endl;	
                 //std::cout<<"MPI_COMM_WORLD: "<<MPI_COMM_WORLD<<std::endl;		       

    		delete[] buffer; 
    		delete[] temp_buffer;
         	
    	}


   else
   	  {
   	  std::cout<<"option for dataset generation not supported"<<std::endl; 	
   	  return;	
   	  }	                          
   }

    
   
}

//Efthymios 2023 END denote step response and noise type in header