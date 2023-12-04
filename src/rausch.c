/*******************************************************************
 *  Patch Clamp
 *  2003 Tobias Huth  
 *  CAU KIEL 
 *  FAU Erlangen-Nürnberg
/*******************************************************************/


#include "rausch.h"
#include "bessel.h"		//Tobias 2021 included for 4 pole Bessel LP filter
#include "fft.h"			//2023

TRausch::TRausch()
{
initialized = false;
fitinaction = false;
/*
erzeugt mit matlab:
[b,a] = besself(4,31415.927);
[y,x,t] = step(b,a)  
y = rot90(y)
save bessel_f.mat y
*/
length_noiseseries = 32000;
}

TRausch::~TRausch()
{
if (initialized)
 delete_array();	
}


//*******************************************************************
void TRausch::delete_array()
{
delete[] tnoiseseries;
}


//*******************************************************************
void TRausch::initial_array()
{
 if (tnoiseseries == 0)				//2021 Tobias added this line for fixing memory leak		
 		{
			tnoiseseries = new double[length_noiseseries];
		}	
}


//*******************************************************************
void TRausch::create_distribution()
{
 int64_t i;
 tnoiseseries_padding = new double[length_noiseseries+2*extralength+6]; //create temporary noise series that will be filtered and transfered to noiseries
 for(i=0; i < (length_noiseseries+2*extralength); i++)
  { 
   tnoiseseries_padding[i] = Daten.snd_01(Daten.random_engine);
  }   
   
 for(i=0; i < length_noiseseries; i++) 
 		tnoiseseries[i] = tnoiseseries_padding[i+extralength];
}


//*******************************************************************
void TRausch::create_lowpass(unsigned int SamplingFrequency, unsigned int FilterFrequency)
{		
 int64_t i;
 double mean, sigma;
 double cutoff;
 
//initialize bessel filter 
cutoff = 2.0 *(double) FilterFrequency/(double) SamplingFrequency;			//calculate cutoff
bessel *lpfilter = initialize_filter(4, cutoff, length_noiseseries, extralength-1);		//order of Bessel filter, cutoff frequency, length?

//apply LP Bessel filter to white noise
filter_signal(tnoiseseries, tnoiseseries_padding, lpfilter, length_noiseseries);		//data, padding, filter, length	
 
for(i=0; i < length_noiseseries; i++)
   tnoiseseries[i] = tnoiseseries_padding[i+ extralength];
  
delete[] tnoiseseries_padding; //get rid off this memory hog
free_filter(lpfilter);

//Einstellung der neuen Standardabweichung
mean = 0;
sigma = 0;
for(i=0; i < length_noiseseries; i++)
  mean += tnoiseseries[i];
mean = mean / length_noiseseries;

for(i=0; i < length_noiseseries; i++)
  sigma += (tnoiseseries[i]- mean) * (tnoiseseries[i]-mean);
sigma = sqrt (sigma / (length_noiseseries-1));  


for(i=0; i < length_noiseseries; i++)
   tnoiseseries[i] = tnoiseseries[i]  / sigma;	

initialized = true;
}



//*******************************************************************
void TRausch::load_noisefile(char* Filename)
{
 unsigned long int i;	
 unsigned short  dummy;
 double mean, sigma;
   	
 std::ifstream in(Filename, std::ifstream::binary);	
 if (!in)
  {
   std::cout<<"Cant load noisefile:"<<Filename<<" ,loading aborted!!!!!!!!!!!!!"<<std::endl;
   exit(0); 
  }
 else
  {
 
   in.seekg(0, std::ios::end);                                // set to end of time series
   	

   length_noiseseries = ((long int)in.tellg()) / 2; //Tobias
   std::cout<<"###################################### length_noiseseries: "<<length_noiseseries<<" Daten.Anz_Samples: "<<Daten.Anz_Samples<<std::endl;
 if (length_noiseseries < Daten.Anz_Samples)
 	{
 	std::cout<<"Cant load noisefile:"<<Filename<<" ,too short! loading aborted!"<<std::endl;
   exit(0); 
  }	
  
  length_noiseseries = Daten.Anz_Samples; 
  in.seekg(0, std::ios::beg);                                  // streampos auf ersten wert

  tnoiseseries = new double[length_noiseseries];          // pointer auf neue daten

  for (i=0; i < length_noiseseries; i++)
    {
     in.read((char *)&dummy,2);           // Werte sind 2byte gross 
     tnoiseseries[i] = (double) dummy + Daten.uni_11(Daten.random_engine);  //Bitrauschenn durch antialiasing vermeiden: VNI = random(-1..+1) 	
    }   
  std::cout<<"loaded noisefile:"<<Filename<<" ,samples: "<<length_noiseseries<<std::endl;   
//Einstellung der neuen Standardabweichung, normieren auf 1
  mean = 0;
  sigma = 0;
  for(i=0; i < length_noiseseries; i++)
    mean += tnoiseseries[i];
  mean = mean / length_noiseseries;

  for(i=0; i < length_noiseseries; i++)
    sigma += (tnoiseseries[i]- mean) * (tnoiseseries[i]-mean);
  sigma = sqrt (sigma / (length_noiseseries-1));  
  std::cout<<"mean: "<<mean<<" sigma: "<<sigma<<std::endl; 
  std::cout<<"loaded data:"<<length_noiseseries<<" from:"<<Filename<<" with mean:"<<mean<<" and sigma:"<<sigma<<std::endl;	
 
  for(i=0; i < length_noiseseries; i++)
   {
    tnoiseseries[i] = (tnoiseseries[i] - mean)  / sigma;  
   } 
  initialized = true;
  }	

}
//*******************************************************************
void TRausch::load_noise_spectrum(char* Filename)		     				//2023	
	{
		double dummy;
	  std::ifstream in(Filename, std::ifstream::binary);	
 		if (!in)
  		{
   		 std::cout<<"Cant load noise spectrum:"<<Filename<<" ,loading aborted!!!!!!!!!!!!!"<<std::endl;
   		 exit(0);
  		}
 		else	
 			{
 			 initial_array();	
 			 in.seekg(0, std::ios::end); 	
 			 length_spectrum = in.tellg() / sizeof(double);
 			 if (2*length_spectrum < Daten.Anz_Samples)     //Efthymios 2023 old: if (length_spectrum < Daten.Anz_Samples) //the algorithm expects the spectrum up to f_sampling/2 for the real and imaginary part, including the zero frequency
 			 		{
                     std::cout<<"length_spectrum: "<< length_spectrum<<std::endl;
 			 		 std::cout<<"lenght of spectrum data is shorter than the time series, aborting-------------------->"<<std::endl;
 			 		 exit(0);
 			 		} 	 	
 			 noise_spectrum.resize(length_spectrum, 0);	
 			 in.seekg(0, std::ios::beg);   
 			 for( long int i=0; i < length_spectrum; i++)	
 			   {
 			   	in.read((char *) &dummy,sizeof(double)); 	
 			   	noise_spectrum[i] = {dummy,0};
 			   	//std::cout<<i<<":"<<noise_spectrum[i]<<std::endl;
 			   }	
 			 	
 			 std::cout<<"loaded spectral file: "<<Filename<<" with "<<length_spectrum<<" data points"<<std::endl;	
 		  }	
	}


//*******************************************************************
void TRausch::create_spectral_noise()
	{
	 double dummy = 1;
	 std::complex <double> dummy1;
	 long int i;
	 double mean;
	 double sigma;	
	 const double Pi = 3.14159265358979323846;
	 	
	noise_spectrum_calc.resize((length_spectrum - 1) * 2, 0);

	 for(i=0; i < length_spectrum; i++)
 		{
 	    dummy = Daten.uni_01(Daten.random_engine) * 2 * Pi;
 		dummy1 = {cos(dummy), sin(dummy)};
		noise_spectrum_calc[i] = dummy1 * noise_spectrum[i];
 		}
			   
	 for(i=0; i < length_spectrum - 1; i++)
	     {
		 noise_spectrum_calc[(length_spectrum - 1) * 2 - 1 - i] = {noise_spectrum_calc[i + 1].real() , noise_spectrum_calc[i + 1].imag() * (-1)};
		 }
       
	 noise_spectrum_calc = fft::ifft(noise_spectrum_calc); 

  for(i=0; i < length_noiseseries; i++)
  	{
  	 tnoiseseries[i] = noise_spectrum_calc[i].real();
     mean += tnoiseseries[i];    
    } 
  mean = mean / length_noiseseries;

  for(i=0; i < length_noiseseries; i++)
    sigma += (tnoiseseries[i]- mean) * (tnoiseseries[i]-mean);
    sigma = sqrt (sigma / (length_noiseseries-1));  
  //std::cout<<"mean: "<<mean<<" sigma: "<<sigma<<std::endl; 	
  for(i=0; i < length_noiseseries; i++)
   {
    tnoiseseries[i] = (tnoiseseries[i] - mean)  / sigma;  
   } 	
 			   	
  }
//*******************************************************************

void TRausch::save_noisefile(char Filename[255])
{
int	i;
int  dummy;
	
std::ofstream out(Filename, std::ios::out|std::ios::binary);

out.write((char *)&length_noiseseries,4);

for(i=0; i < length_noiseseries; i++)
  {
   dummy = int(tnoiseseries[i]);
   out.write((char *) &dummy,4);
  } 

out.close();	
}	

//*******************************************************************


void TRausch::makenoise(unsigned long int length,						//TObias 2022
  			unsigned int sampling_frequency,
  			unsigned int filter_frequency, 
  			bool resample,
  			char* filename)
{	

extralength = (unsigned long int)sampling_frequency / filter_frequency * 4;
length_noiseseries = length;



if (Daten.Geneticfitparameter.noise_generation == 2)																						//load recorded noise
   {			
    if (Daten.Geneticfitparameter.file_noise_series_loaded == false) //Tobias 2022
      {
       load_noisefile(filename);
       Daten.Geneticfitparameter.file_noise_series_loaded = true;
      } 
   }     
   
else if (Daten.Geneticfitparameter.noise_generation == 1)																				//generate noise with spectrum
   {			
    if (noise_spectrum.size() == 0) //Tobias 2022
      {
       load_noise_spectrum(filename);
      } 
    create_spectral_noise();  
   }      
   
   
else if (Daten.Geneticfitparameter.noise_generation == 0)																				//generate simple simulated noise
	 {
	 	if (resample)
	 		{
	 		 initial_array();
       create_distribution();
       create_lowpass(sampling_frequency, filter_frequency);
       resample = false;          	
      }
    else
    	{
    	 if (!fitinaction)
    	 {
    	 	initial_array();
        create_distribution();
        create_lowpass(sampling_frequency, filter_frequency);
        fitinaction = true;
    	 }		
      }	    
	 }
else	 
	{
	 std::cout<<"no valid method for noise generation, aborting"<<std::endl;
	 exit(0);
  }
}






