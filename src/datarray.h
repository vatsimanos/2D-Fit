/* Diese Includedatei enthaelt das Basisobjekt aller Objekte, die auf ein
   Datenfeld zur�ckgreifen m�ssen. Das Datenfeld kann sich auf dem Bildschirm
   darstellen, sein Maximum bestimmen, sich ein eine Datei schreiben, Daten
   aus einer Datei einlesen und alle St�tzpunkte auf einen bestimmten Wert
   setzen
*/

#if !defined _datarray_h
#define _datarray_h
#include<string.h>
#include<iostream>
#include<iomanip>  //.h weg
#include<fstream>	//.h weg
#include "daten.h"

#ifdef unix
//  #define ios::bin 0   // dirty hack for getting rid of the bin problem when
                     // compiling on ibm and unix platforms
#endif

//     Zum einlesen von Datenfiles aus Matlab
typedef struct {
     unsigned int type;        //2019 converted from long
     unsigned int mrows;
     unsigned int ncols;
     unsigned int imagf;
     unsigned int namlen;
} Fmatrix;

   typedef struct {
     int LX,  RX, OY, UY;
   } TBildschirmausschnitt;
//   TBildschirmausschnitt ganzer_Bildschirm;

   struct TDatenfeld                      // abstraktes Objekt
   {

   public:
     unsigned long int Stuetzpunkte, Max_Daten, Min_Daten;
     char *DADaten;  //Tobias: NAME ge�ndert alt:Daten
/*
   um Matlabfiles zu laden jetzt die explizite Methode
              TDoublefeld::load_matlab_file(char *FName);
    benutzen !!!
     TDatenfeld(int Handle, char *n_name="");
*/
     TDatenfeld(TDatenfeld& N_Datenfeld, char *n_name=(char *)"");
     TDatenfeld(unsigned long int N_Stuetzpunkte, int Datengroesse,
                char *n_name=(char *)"");
     virtual ~TDatenfeld();
     void * nr(unsigned long int Pos = 1);
     void neue_Dimension(unsigned long int N_Stuetzpunkte);
     void Set_Max_Daten(unsigned long int Startpos = 0,
                        unsigned long int Endpos = 0);
     void Set_Min_Daten(unsigned long int Startpos = 0,
                        unsigned long int Endpos = 0);
     virtual double D_Wert(unsigned long int Pos) = 0;
     virtual int Int_Wert(unsigned long int Pos) = 0;
     double Max() { return D_Wert(Max_Daten); };
     double Min() { return D_Wert(Min_Daten); };
     void Darstellen(TBildschirmausschnitt Fenster,
                     unsigned long int Startpos,
                     unsigned long int Endpos,
                     int Farbe,
                     float Faktor = 1);
     int ASCII_file(const char* Filename);
//     TDatenfeld(FILE *S);
//     void Store(FILE *S);
    friend std::ostream& operator << (std::ostream& str, TDatenfeld& Datenfeld);
    friend std::istream& operator >> (std::istream& str, TDatenfeld& Datenfeld);
   protected :
    int eigene_Daten;
    int Datengroesse;
    char *name;
    void free_Daten();
   };
   typedef TDatenfeld *PDatenfeld;

   struct TDoublefeld : TDatenfeld
   {
	
   public:
/*
     TDoublefeld(int Handle, char *n_name="") :
       TDatenfeld(Handle, n_name) {D=(double*)Daten;};
*/
     TDoublefeld(unsigned long int N_Stuetzpunkte, char *n_name=(char *)"") :
       TDatenfeld(N_Stuetzpunkte, sizeof(double), n_name) {D=(double*)DADaten;};  //TOBIAS
     void neue_Dimension(unsigned long int N_Stuetzpunkte)
       { TDatenfeld::neue_Dimension(N_Stuetzpunkte); D=(double*)DADaten; };       //TOBIAS
     int ASCII_file_lesen(const char* Filename);
     int load_matlab_file(char *FName);
     virtual double D_Wert(unsigned long int Pos);
     virtual int Int_Wert(unsigned long int Pos);
     double& operator [] (unsigned long int Pos);
     TDoublefeld& operator = (TDoublefeld& Feld);
     TDoublefeld& operator += (double Summand);
     TDoublefeld& operator *= (double Faktor);
     double Summieren(int anfang = 0,int ende = 0);
     void Funktion(double (*F)(double));
     void Reset(double Wert = 0);
     double* nr(unsigned long int Pos=1){return (double*)TDatenfeld::nr(Pos);};
   protected :
     double *D;
   };
   typedef TDoublefeld *PDoublefeld;

   struct TFloatfeld : TDatenfeld
   {
   public:
/*
     TFloatfeld(int Handle, char *n_name="") :
       TDatenfeld(Handle, n_name) {F=(float*)Daten;};
*/
     TFloatfeld(unsigned long int N_Stuetzpunkte, char *n_name=(char *)"") :
       TDatenfeld(N_Stuetzpunkte, sizeof(float), n_name) {F=(float*)DADaten;};  //TOBIAS
     void neue_Dimension(unsigned long int N_Stuetzpunkte)
       { TDatenfeld::neue_Dimension(N_Stuetzpunkte); F=(float*)DADaten; };      //TOBIAS
     int ASCII_file_lesen(const char* Filename);
     int load_matlab_file(char *FName);
     virtual double D_Wert(unsigned long int Pos);
     virtual int Int_Wert(unsigned long int Pos);
     float& operator [] (unsigned long int Pos);
     TFloatfeld& operator = (TFloatfeld& Feld);
     TFloatfeld& operator += (float Summand);
     TFloatfeld& operator *= (float Faktor);
     void Funktion(float (*Funk)(float));
     void Reset(float Wert = 0);
     float* nr(unsigned long int Pos=1){return (float*)TDatenfeld::nr(Pos);};
   protected :
     float *F;
   };
   typedef TFloatfeld *PFloatfeld;

   struct TIntegerfeld : TDatenfeld
   {	
   public:
     TIntegerfeld(unsigned long int N_Stuetzpunkte, char *n_name=(char *)"") :
       TDatenfeld(N_Stuetzpunkte, sizeof(int), n_name) {I=(int*)DADaten;};
     void neue_Dimension(unsigned long int N_Stuetzpunkte) {
        TDatenfeld::neue_Dimension(N_Stuetzpunkte);
        I=(int*)DADaten;
     };
     virtual double D_Wert(unsigned long int Pos);
     virtual int Int_Wert(unsigned long int Pos);
     int& operator [] (unsigned long int Pos);
     void Reset(int Wert = 0);
   protected:
     int *I;
   };
   typedef TIntegerfeld *PIntegerfeld;

   struct TShortIntfeld : TDatenfeld
   {	
   public:
     TShortIntfeld(TShortIntfeld& N_ShortIntfeld, char *n_name) :
       TDatenfeld(N_ShortIntfeld, n_name){};
     TShortIntfeld(unsigned long int N_Stuetzpunkte, char *n_name=(char *)"") :
       TDatenfeld(N_Stuetzpunkte, sizeof(short int), n_name) {
        SI=(short int*)DADaten;
     };
     void neue_Dimension(unsigned long int N_Stuetzpunkte) {
        TDatenfeld::neue_Dimension(N_Stuetzpunkte);
        SI=(short int*)DADaten;
     };
     virtual double D_Wert(unsigned long int Pos);
     virtual int Int_Wert(unsigned long int Pos);
     short int& operator [] (unsigned long int Pos);
     void operator *= (double Faktor);
     void operator += (short int Summand);
     void Reset(int Wert = 0);
   protected:
     short int *SI;
   };
   typedef TShortIntfeld *PShortIntfeld;

#endif  /* defined _datarray_h */
