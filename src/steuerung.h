/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/


#ifndef _STEUER_H_
#define _STEUER_H_

//#include "glib.h"

#include "declare.h"

// funktion zum steuern der erstellung der dwelltime-dateien mit vogegebenen werten (aus simulations-set-datei)
void  
Steuerung(char* setfilename, short histopts, double koop, bool test);

// zum simulieren von 2d dwell time histogrammen (mit simulat!)
void
Steuerung(char* setfilename, char* dwellfilename, short histopts, double minevents);


// zum testen des expfits mit einem c-o modell
void Steuerung(char* setfilename, short histopts);

struct TChannel_type 
{
  short  n_Channels;
  short  i_null;
  short  i_channel;
  short  n_states;
  TRate_matrix  Simulation_Matrix;
  TRate_matrix  StartValues;

  TChannel_type(short an_channels, 
                TRate_matrix aSimulation_Matrix, 
				TRate_matrix aStartValues,
                short ai_null, 
				short ai_channel, 
                short an_states);
};

#endif /* _STEUER_H_ */

