//*******************************************************************
// *  Patch Clamp
// *  2002 Tobias Huth  
// *  CAU KIEL 
//*******************************************************************/
#include <cmath>
#include <complex>
#include "gtkfunc.h"
#include "datarray.h"

class TRausch
{
	
 private:
//schon mal initialisiert? 
 bool initialized;
 bool fitinaction;

//extralaenge um eine zirkul�re noiseserie zu erzeugen, damit der Filter am Anfang der serie schon "eingeschwungen" ist
 unsigned long int extralength; 		//Tobias 2022 long int
 // char storefilename[255];				//Tobias 2022
 


 void delete_array();
 void initial_array();
 void create_distribution();
 void load_noisefile(char name[255]);
 void save_noisefile(char name[255]);
 void create_lowpass(unsigned int SamplingFrequency, 
 		     						 unsigned int FilterFrequency);
 void create_spectral_noise();											//2023
 void load_noise_spectrum(char name[255]);		     	//2023		
 
 
 public:
//array für die Rauschreihe
 double *tnoiseseries;
 double *tnoiseseries_padding; //Tobias 2021
 std::vector<std::complex<double>> noise_spectrum;				//2023
 std::vector<std::complex<double>> noise_spectrum_calc;		//2023
 long int length_spectrum; 																//2023
 

//Länge der Rauschreihe 
 unsigned long int length_noiseseries; 					//Tobias 2022 long int
  //unsigned long int length_noiseseries_load; 		//Tobias 2022 long int
 

 TRausch();
 ~TRausch();
 

 
 void makenoise(unsigned long int length,					//Tobias 2022 long int
  		unsigned int sampling_frequency,
  		unsigned int filter_frequency, 
  		bool resample,
  		char filename[255]);

};









