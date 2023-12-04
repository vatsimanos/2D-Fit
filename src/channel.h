#ifndef _channel_h
#define _channel_h

#include<stdlib.h>
#include<iostream>	//.h weg
#include<math.h>
#include<sys/types.h>
#include<time.h>
#include"datarray.h"
#include"matrix2.h"
#include "daten.h"

#if !defined(M_PI)
#define M_PI 3.141592654
#endif

struct TTSprung     //Tobias TSprung->TTSprung
{
  double Zeit, Differenz;
  int alter_Zustand,neuer_Zustand;
  double altes_Zustandsniveau, neues_Zustandsniveau;
  TTSprung* naechster;  //Tobias Name ge�ndert
};

struct TKanal
{
  Matrix Ratenmatrix, Uebergangsmatrix;
  Vektor Stromvektor, Bes_Wahrsch, initial_state_value;
  // **N_Matrix ist die Matrix der Ratenkonstanten
  // **N_Stromvektor ist die Matrix die den Zustaenden einen Strom zuordnet
  // sie hat nur eine Reihe
  // N_Zustaende ist die Anzahl der Zustaende des Kanals
  // N_Sigma ist das Sigma des Rauschens
  TKanal(int N_Zustaende, double **N_Stromvektor, double **N_Matrix,
         //double N_Sigma = 0, int Startzustand = 1);
	 double  **N_initial_state_value, double N_Sigma = 0); 
  ~TKanal();
  void naechster_Zustand();
  unsigned short int Strom(double Zeit);
  double Bes_W(short int Niveau);
  double Bes_W(int Zustand);
  int Zust() { return Zustand;};
  Matrix Uebergaenge();
  TTSprung *jumps;   //Tobias Name ge�ndert
  protected :
  int Zustaende, Zustand, alter_Zustand, Ziel, k;
  // in t wird von naechster_Zustand der Zeitpunkt des naechsten
  // Ereignisses (Sprunges) gespeichert
  // Ereignis ist in diesem Fall ein Zustandswechsel
  // 
  double t, Verweildauer, Sigma;
  // Im Ratenvektor ist die Summe aller vom jeweiligen Zustand wegfuehrenden 
  // Ratenkonstanten gespeichert
  // also Ratenvektor(i) liefert die Summe aller vom i-ten Zustand 
  // wegfuehrenden Ratenkonstanten
  Vektor Ratenvektor;
  int Startzustand(double  **N_initial_state_value);
};

struct TKanal_gefiltert : TKanal
{	
  double tau, Ausgang;
  TKanal_gefiltert(int N_Zustaende, double **N_Stromvektor, double **N_Matrix,
                   //double f_3_dB, double N_Sigma = 0, int Startzustand = 1);
		   double **N_initial_state_value, double f_3_dB, 
		   double N_Sigma = 0);
  unsigned short int Strom(double Zeit);
};




struct TKanal_digfi : TKanal
{	
  TKanal_digfi(int N_Zustaende, double **N_Stromvektor, double **N_Matrix,
               double N_dT_Sprung, TDoublefeld& N_Sprungantwort,
	       double  **N_initial_state_value, double N_Sigma = 0);
  ~TKanal_digfi();
  unsigned short int Strom(double Zeit);
  void clear();
  Matrix Differenzenmatrix;
  protected :
  TDoublefeld Sprungantwort;
  TTSprung *Spruenge, *letzter_Sprung;   //Tobias
  short int No_Zeiten;
  double dT_Sprung, Summe;
  int alter_Zustand;
};

struct TKanal_mehrfach_gefiltert : public TKanal
{
  unsigned short int N_Taus, *Ergebnis;
  Vektor tau, Ausgang;
  TKanal_mehrfach_gefiltert(int N_Zustaende, double **N_Stromvektor,
                            double **N_Matrix, Vektor& f_3_dB,
                            //double N_Sigma = 0, int Startzustand = 1);
			    double  **N_initial_state_value,
			    double N_Sigma = 0);
  ~TKanal_mehrfach_gefiltert();
  unsigned short int * Strom(double Zeit);
};

#endif /* !defined _channel_h */
