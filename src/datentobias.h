/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  2000 Erweiterung von Tobias Huth
 *  2011 Erweiterung von Tobias Huth
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

// fuegt weitere globale Parameter in die Struktur Daten(in daten.h) ein
//Daten.lauf: Iterator fuer die verkettete Liste Markerbits
//Daten.markernumber fuer die Position des Markers in der Liste
//Daten.markerstart fuer die Markierung des ersten Aufrufs der Menues
//Daten.iwvisible fuer die Sichtbarkeit des infowindows
//Daten.nrselected fuer die Anzahl der insgesamt ausgewaehlten Abschnitte
//Daten.paintinaction ist dafuer, dass die paintprozedure nicht nochmal waehrend des Malens aufgerufen wird
//Daten.TwoD_fitinaction dient zum Abfangen von zeitaufwendigen Ausgaben



std::vector<TMarkBits>::iterator lauf;
int markernumber;  
bool markerstart; 
bool iwvisible;
int nrselected;
bool paintinaction;        
bool TwoD_fitinaction;
char logfilename[255];
char logdwellfilename[255];
int loeschanzahl;
int   loeschvalue; 
bool geneticfit; 
unsigned short calculate_error;
bool firstsimul;
bool simulate_memory;
char errormessage[255];
int dwellmatrix[11][51];
double logx;
int cut_marked_save;
int cut_start;
double original_ampl_histogram[65536];

TSimSettings SimSettings;

TImportSettings ImportSettings;

struct TSaveTimeseries		// Tobias 2020 data imported from setfile to control automated saving
 {
 	int automated;			// Tobias 2020 programm starts saving timeseries if in multiprocessor mode
 	std::string	 folder;		// Tobias 2020 in this folder
 	int  number;				// Tobias 2020 number times
 	bool random_rate;		// Tobias 2020 if true randomises kij
 	bool random_lvl;		// Tobias 2020 if true randomises lvl
 	double min_noise;
 	double max_noise;
 	bool random_noise;
 };

TSaveTimeseries SaveTimeseries; //Tobias 2020

TFitparameter Fitparameter;
TGeneticfitparameter Geneticfitparameter;
 
bool dependency_plot_draw_diff;
bool matlab_file_loaded;						//Tobias 2022 flag for preventing multiple loading of stepresponse file
double *stepresponse_loaded;								//Tobias 2022 storage for stepresponses
long unsigned int Stuetzpunkte_loaded;				//Tobias 2022 storage for stepresponses

//new random generator with distributions
std::random_device os_seed; 		//for os dependent seed, does now work with mingw
std::mt19937 random_engine;		//Mersenne twister 64 bit
std::uniform_real_distribution<double> uni_01{0,1};
std::uniform_real_distribution<double> uni_11{-1,1}; 
std::normal_distribution<double> snd_01{0.0,1.0};

