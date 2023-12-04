/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  2000 Tobias Huth 
 *  based on uwe kirst's daytest  
 *  CAU KIEL
 *  Genetic Fit: Copyright 1995-1996 Massachusetts Institute of Technology
 *******************************************************************/
  #include "chrono"
  
   TTimeseries TobyTimeseries; 
   TRausch Rausch;
   TSetfile Setfile;

   double gewichtsfaktor;
   bool MaximumLikelihoodfit;
   double abknickfrequenz;
   double markerzustand;	// Hier noch zu kurzem integer umwandeln!!!
   char mlh_logfile[255];

   auto timer_start = std::chrono::system_clock::now();
   int myID;
   int numprocs;
   char processor_name[255];
   int namelen;



short
load_setfile(char *Filename)		//Tobias 2013
 {
  double version;
  char textbuf[255];
  char* textbuf1;
  char* test_filename;
  int test;
 
  TDoublefeld matlab_test(0,(char*)"stepresponse");
  //std::ofstream logout(Daten.logfilename,std::std::ios::app);       
  version = Setfile.read_setfile(Daten.SimSettings.loadname);     
 
    if (version == 0)
     {
      std::cout<<Daten.errormessage<<std::endl;
      strcpy(textbuf,Daten.errormessage);
      Application.statustext(textbuf);
      return -1;
     }
   
    Daten.a_fit.sigma = Setfile.get_sigma();
    Daten.a_fit.i_null = Setfile.get_zero_current();	
    Daten.a_fit.i_channel = Setfile.get_single_channel_current();	
    Daten.a_fit.n_channels = Setfile.get_no_channels_total();	
    Daten.SimSettings.sampling_frequency = (unsigned int)floor(1/Setfile.get_delta_t()+0.5);
  
/*    Tobias 2022
    textbuf1 = Setfile.get_rausch_filename(NULL); 

    if (strspn(textbuf1,"none") > 2)
       {
        Daten.SimSettings.noisefile[0]=0;
        if (!Daten.Geneticfitparameter.mpi_support)
        	  std::cout<<"generating artificial noise"<<std::endl; 
       }
    else
       {
    	  strcpy(Daten.SimSettings.noisefile, textbuf1);
    	 }
*/ 
    
    test_filename = Setfile.get_stepresponse_filename(NULL); 
    test  = matlab_test.load_matlab_file(test_filename); 
    if (test == 1)
  	{
   	 std::cout<<Daten.errormessage<<std::endl;
   	 strcpy(textbuf,Daten.errormessage);
   	 Application.statustext(textbuf);
   	 return-1;
        } 
   
    //else    
    	 //sprintf(Daten.SimSettings.stepresponsefile,test_filename);

    	if (!Daten.Geneticfitparameter.mpi_support)	
    		std::cout<<"Setfile loaded: "<<Daten.SimSettings.loadname<<std::endl;
    //std::cout<<"Noisefile loaded: "<<Daten.SimSettings.noisefile<<std::endl;
    //std::cout<<"Stepresponsefile loaded: "<<test_filename<<std::endl;
    //logout<<"Setfile loaded: "<<Daten.SimSettings.loadname<<std::endl;
    //logout<<"Noisefile loaded: "<<Daten.SimSettings.noisefile<<std::endl;
    //logout<<"Stepresponsefile loaded: "<<test_filename<<std::endl; 
    
  return 0;
 }
 
