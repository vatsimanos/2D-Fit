  

TTimeseries *TSetfile::simulate_timeseries(TTimeseries *ts,						//Tobias 2017 defaults weg
					   TRausch *rausch,
					   char *signal_filename,
					   double dT_filter)
{
  //Tobias unsigned short int *timeseries;  
  int i, j, k=0;
  double D;
  double test_freq=0, test_freq_ampl=0 , test_freq_phase=0;
  double **Stromvektor;
  double **init_sv;
  double strom_alt=0;
  //,strom_neu=0;
  int Zust_array[256];
  //Tobias:long int nu=0;
  int add_signal_flag=0;
  int g_flag=0;
  int I;
  int AD_Bits = 16;		//Tobias 2022
  int Stuetzpunkte = (int)pow(2, AD_Bits);
  short unsigned int Min = (short unsigned int) pow(2, AD_Bits), Max = 0;
  unsigned long int low_boundary = 0, high_boundary =0;
  
  unsigned long int slength;
  enum e_channel_type {one_pole, digfi};
  enum e_channel_type filter_type = digfi;
  TKanal **Kanal;
  std::ifstream signal_file;
  //ofstream tjumps;
  //char *FNtjumps;
  TDoublefeld stepresponse(0,(char *)"stepresponse");
  //TDoublefeld Amplhis(Stuetzpunkte, (char *)"Amplhis");				//Tobias 2022 removed
  unsigned long int nr_error = 0;			//Tobias 2018
  stepresponse.load_matlab_file(stepresponse_file); 	
   
#ifdef debug_simulat
  cout << rausch << endl;
  for (unsigned long int ll=1;ll<=stepresponse.Stuetzpunkte;ll++)
  cout << ll << " : " << stepresponse[ll] << endl;	
#endif
  /*
  if ((ts!=NULL)||(noise!=NULL)){
    err.warning(NOT_CALLED_WITH_NULL,"TSetfile::simulate_timeseries","ts or noise");
  }
  */
  
  ts->set_num_samples(no_samples);
  ts->set_delta_t(delta_t);
  ts->set_filename(datafile);
  ts->set_f_3dB(f_3dB_digfi);
  rausch->makenoise(Daten.SimSettings.noise_length, 
  		    Daten.SimSettings.sampling_frequency, 
  		    Daten.SimSettings.filter_frequency, 
  		    Daten.SimSettings.resample,
  		    Daten.Geneticfitparameter.file_noise_series);  		//Tobias 2022


  if (filter_type == digfi)
    {
      double i_max,f_3dB_alt,i_3dB,tau_alt,delta_i;
      int n=1;
#ifdef debug_mem_simulat
      cerr << "TSetfile::simulate_timeseries new TKanal* " 
	   << sizeof(TKanal_digfi)*no_channels_total << endl;
#endif
      //Kanal=new TKanal*[no_channels_total];
    
      Kanal=new TKanal*[sizeof(TKanal_digfi)*no_channels_total];
      // Beruecksichtigung der Eckfrequenz (ab SV3.0) 
      // Verhaeltnis von 3dB Freq zu Abtastfrequenz
      // Hoehe des Sprunges im Matlabfile
      i_max=stepresponse[stepresponse.Stuetzpunkte];
      // Hoehe des Stromes bei der Eckfrequenz ??
      i_3dB=i_max*(1-1/sqrt(2));
#ifdef debug_simulat
      cout << "i_max : " << i_max << "\ti_3dB : " << i_3dB << endl;
#endif
      // i_3dB wird zwischen dem n und (n+1)ten Punkt errreicht
      while (stepresponse[n]<i_3dB){
#ifdef debug_simulat
	cout << "n : " << n << "\tstepresponse[n] : " << stepresponse[n] 
	     << endl;
#endif
	n++;
      }
      n--;
      delta_i=stepresponse[n+1]-stepresponse[n];
#ifdef debug_simulat
      cout << "n : " << n << "\tr[n] : " << stepresponse[n] 
	   << "\tdelta_i : " << delta_i << endl;
#endif
      // Berechnung der Eckfrequenz des Filters das durch das
      // Matlabfile beschrieben wird ist analytisch nicht ohne weiteres 
      // moeglich !!!!!!!!!!!! (geht nur beim 1-pooligen)
      // Beim 4 poligen Besselfilter dauert es 7 mal so lange 
      // bis bei der stepresponse der 1/sqrt(2)te Teil der Hoehe 
      // errreicht ist wie beim einpoligen Filter
      tau_alt=delta_t*(n+(i_3dB-stepresponse[n])/delta_i)/7;
#ifdef debug_simulat
      cout << "tau_alt : " << tau_alt << "\tf_3dB_alt : " 
	   << 1/(2*M_PI*tau_alt) << endl;
#endif 

      f_3dB_alt=1/(2*M_PI*tau_alt);
      // dT_filter ist der Zeitabstand der zwischen zwei Stuetzpunkten des 
      // Matlabfiles angenommen werden muss, damit man auf die vom
      // benutzer festgelegte Eckfrequenz kommt
      // wenn dT_filter != 0 , dann wurde die -T option nicht benutzt
      // fuer eine Eckfrequenz von 50kHz betraegt der Zeitabstand zwischen
      // zwei Punkten im bessel_4.mat File 1.0754e-6 s (laut Auskunft von Philipp)
      // er hat sich die stepresponseen mit Matlab erzeugt und die Werte dann
      // ausgemessen
      //abknickfrequenz des Filters
      if (dT_filter==0)
	  
      dT_filter = (1 / (double)Daten.SimSettings.sampling_frequency);  // Efthymios 2023 sampling periode der gemessenen step response eingeben.
      ////old: dT_filter=1.0754e-6  * 50000 * 0.25 / ((double)Daten.SimSettings.filter_frequency/(double)Daten.SimSettings.sampling_frequency) / f_3dB_digfi;

  //cout<<"f_3db_digfi: "<<f_3dB_digfi<<endl;
  //cout<<"dT_filter: "<<dT_filter/1.0754e-6<<endl;

    }
  else{
#ifdef debug_mem_simulat
    cerr << "TSetfile::simulate_timeseries new TKanal* " 
	 << sizeof(TKanal_gefiltert)*no_channels_total << endl;
#endif
    Kanal=new TKanal*[sizeof(TKanal_gefiltert)*no_channels_total];}
 

  // Fuer jede Kanalsorte ...
  for (i = 1; i <= no_channel_types; i++)
    {
      // Stromvektor ist Matrix mit einer einer Reihe und
      // channel_type[i]->no_states Spalten
      // die Indizes fangen bei 1 an
      // in ihm wird der Strom in Abhaengigkeit vom Zustand gespeichert 
      Stromvektor = real_matrix(1,1,1,channel_type[i]->no_states);
      for (j=1; j<=channel_type[i]->no_states; j++) 
       {
       	if (Daten.Geneticfitparameter.levelfit == true)				//Tobias 2012
       	   Stromvektor[1][j] = Daten.Geneticfitparameter.stromvektor[j];	//Tobias 2012
       	else   
            Stromvektor[1][j] = channel_type[i]->level[j]; 
        //cout<<endl;    
        //cout<<"Stromvektor"<<" "<<j<<":"<<Stromvektor[1][j]<<endl;		//Tobias 2012
       }     									//Tobias 2012

      init_sv = real_matrix(1,1,1,channel_type[i]->no_states);
      for (j=1; j<=channel_type[i]->no_states; j++) 
        init_sv[1][j] = channel_type[i]->initial_state_value[j];
      
      // Fuer jeden Kanal einer Sorte ...
#ifdef debug_simulat
      cout << "dT_filter : " << dT_filter << endl;
#endif
      for (j=1; j<=no_channels_type[i]; j++)
       {
	k++;

	if (filter_type == digfi)
	{ 		
	  //channel_type[i]->k sind die Ratenkonstanten des Kanals
	  //delta_t * Verhaeltnis vielleicht 
	  //Verhaeltnis von Abtastrate und 
	  //3dB Freq. des digitalen Filters ausrechnen (siehe
	  // Uebrholt, dT_filter ist richtig siehe oben
#ifdef debug_mem_simulat
	  cerr << "TSetfile::simulate_timeseries new Kanal[" << i << "] \n"; 
#endif
	  Kanal[k] = new TKanal_digfi(channel_type[i]->no_states,
				      Stromvektor, channel_type[i]->k,
				      dT_filter, stepresponse,
				      init_sv); 		      			      
       }     
	else{

#ifdef debug_mem_simulat
	  cerr << "TSetfile::simulate_timeseries new Kanal[" << i << "] \n"; 
#endif

	  Kanal[k] = new TKanal_gefiltert(channel_type[i]->no_states,
					  Stromvektor, channel_type[i]->k,
					  init_sv, f_3dB_digfi);}					  
#ifdef debug_simulat
	cout << "Kanal " << k << " faengt im Zustand " << Kanal[k]->Zust() 
	     << " an " << endl; 
#endif
      } // endfor i
      free_real_matrix(Stromvektor, 1, 1, 1, channel_type[i]->no_states);
      free_real_matrix(init_sv,1,1,1,channel_type[i]->no_states);
    }
  // Es gibt jetzt insgesamt k Kanaele
  if (no_channels_total>255)
    {
      std::cout << "it is not possible to simulate more than 255 channels" << std::endl;
      exit(5);
    }
  strom_alt=zero_current;
  
     
  for (i=1;i<=no_channels_total;i++)
    {
      Zust_array[i]=Kanal[i]->Zust();
#ifdef debug_simulat
      cout << "Zust_array[" << i << "] = " << Zust_array[i] << endl;
#endif
      strom_alt+=Kanal[i]->Stromvektor(Kanal[i]->Zust());
    }
  
  if (strcmp(signal_filename, "") != 0)
    {
      signal_file.open(signal_filename, std::ios::in

		       );
      if (!signal_file)
	{
	  std::cout << "error opening signal file";
	  exit(1);
	}
      signal_file.read((char *)&slength, sizeof(unsigned long int));
      add_signal_flag = 1;
    }
  //Amplhis.Reset();			//Tobias 2022 removed


//Tobias: n�chste Zeilen ausgenommen und durch D=UNI ersetzt;
//  if (random_offset == 0)
//    I = KISS%32767; //TObias:- time(NULL);
//  else
//    I = (random_offset>0)? -random_offset : random_offset;
//  D = ran1(&I);
  
  D = Daten.uni_01(Daten.random_engine);
  // end Init()

  // Simulieren()
  unsigned long int ri, r_pos, s_pos = 0;
  short unsigned int sidummy, signal;
  double t;
  float fdummy;
  
  //Tobias: variabler Beginn f�r das Rauschen, damit bei mehrfacher Simulation kein gerichteter Fehler auftritt
std::uniform_int_distribution<unsigned long int> uni_start{0,rausch->length_noiseseries}; 
r_pos = uni_start(Daten.random_engine);

  
  for (ri = 1, t = delta_t; ri <= no_samples; ri++, t += delta_t) {
    sidummy = 0;
    for (j = 1; j <= no_channels_total; j++) {
      if (filter_type == digfi) {    	
	sidummy += ((TKanal_digfi *) Kanal[j])->Strom(t);
#ifdef debug_simulat
	/*if (sidummy <0) cout <<n"Strom : " << sidummy << "\tKanal : " 
	  << j <<endl;*/
#endif
      } else {
	sidummy += ((TKanal_gefiltert *) Kanal[j])->Strom(t);
      } /* endif */      
    } /* end for j <= no_channels_total */    
    if (add_signal_flag) {
      if (s_pos>=slength) {
	s_pos = 0;
	signal_file.seekg(4,std::ios::beg);
      } /* endif */
      signal_file.read((char *)&signal,2);
      s_pos++;
      sidummy = (short unsigned int) ((sidummy * (signal-2048)) / 1000.0);
      std::cout<<"sidummy3:"<<sidummy<<std::endl;;
    } /* endif */
    /*if (sidummy <0) 
      cout << "Signal : " << sidummy << endl;*/
    if (sigma != 0) {
      if (r_pos >= rausch->length_noiseseries) 
	r_pos = 0;
      if (!g_flag){
	fdummy=rausch->tnoiseseries[r_pos];  //Tobias rausch->rauschcopy
	fdummy *= sigma;
	r_pos++;
      }
      else{
fdummy = Daten.snd_01(Daten.random_engine);		//Tobias 2018 gasdev replaced
	// cout << "gasdev : " << fdummy << endl;
	fdummy *= sigma;
      }
    } else {
      fdummy = 0;
    }
    if (test_freq != 0) {
      fdummy += test_freq_ampl*sin(2*M_PI*t*test_freq +
				   test_freq_phase*M_PI/180);
      } /* endif */
      
    sidummy += (short unsigned int)(fdummy + zero_current + 0.5);
    /*if (sidummy <0) 
      cout << "fdummy : " << sidummy << endl;*/
    if (sidummy <= Stuetzpunkte && sidummy > 0) 
      {
	Daten.SimulatWert[ri]=sidummy;
	Max = (sidummy > Max)? sidummy : Max;
	Min = (sidummy < Min)? sidummy : Min;
	//Amplhis[sidummy]++;
	//cout << sidummy << endl;
      } 
    else 
      { 
      	nr_error++;
      //std::cout<<"sidummy:"<<sidummy<<" Min:"<<Min<<" Max:"<<Max<<" fdummy:"<<fdummy<<" zero_current:"<<zero_current<<std::endl;
       // std::cerr << "current value out of range: " << sidummy << "\ti:\t" << ri << std::endl;
	ri--;
      } /* endif */
      
  if (sidummy < 500)
  	low_boundary++;
  if (sidummy >	65036)
  	high_boundary++;    
    /* .tjumps nicht mit abspeichern
    // Zum abspeichern der Spruenge ...
    // in tjumps
    for (j=1;j<=no_channels_total;j++)
      {
	if(Zust_array[j]!=Kanal[j]->Zust())
	  {
	    strom_neu=strom_alt+(Kanal[j]->Stromvektor(Kanal[j]->Zust())-
				 Kanal[j]->Stromvektor(Zust_array[j]));
#ifdef debug_simulat
	    cout << "Kanal " << j << " : " << Zust_array[j] << "\t" 
		 << Kanal[j]->Zust() << "\t" 
		 << strom_alt << "\t"
		 << strom_neu << "\t"
		 << t << endl;
#endif
	    tjumps << t << "\t\t" << j << "\t" <<  Zust_array[j] << "\t" 
		   << Kanal[j]->Zust() << "\t" 
		   << strom_alt << "\t"
		   << strom_neu << "\t" << ri
		   << endl;
	    strom_alt=strom_neu;
	    Zust_array[j]=Kanal[j]->Zust();
	  }
	
      }
      */
  } /* endfor */
  // end Simulieren
  // ts->set_value(Daten.SimulatWert,no_samples);		//Tobias 2018 kann wohl weg
#ifdef debug_mem_simulat
  cerr << "TSetfile::simulate_timeseries delete[](timeseries) \n"; 
#endif
  //Tobias: delete[](Daten.SimulatWert);
#ifdef debug_mem_simulat
  cerr << "TSetfile::simulate_timeseries delete[](noiseseries) \n"; 
#endif

  for (int i=1;i<=no_channels_total;i++)
    if (Kanal[i]!=NULL){
#ifdef debug_mem_simulat
      cerr << "TSetfile::simulate_timeseries delete(Kanal[" << i << "]) \n"; 
#endif
      if (filter_type == digfi)
	((TKanal_digfi *)Kanal[i])->clear();
// die n�chste Zeile beseitigt ein MEMLEAK, da der Destruktor von TKanal_digfi bei delete scheinbar nicht aufgerufen wird,
// und so matrix nicht freigegeben wird!	
	((TKanal_digfi *)Kanal[i])->Differenzenmatrix.~Matrix();
	
      delete(Kanal[i]);}
#ifdef debug_mem_simulat
  cerr << "TSetfile::simulate_timeseries delete[](Kanal) \n"; 
#endif
  delete[](Kanal);	
  if (nr_error > 0)																												//Tobias 2018
  	std::cout<<"==== warning ==== times sidummy out of range: "<<nr_error<<std::endl;  		     
  if ((low_boundary + high_boundary) > 0)
  	std::cout<<"==== warning ==== "<<low_boundary+high_boundary<<" data points out of range < 500, > 65036"<<std::endl;

  return ts;
}

