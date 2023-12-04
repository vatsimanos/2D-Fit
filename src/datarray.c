/* Diese Includedatei enthaelt das Basisobjekt aller Objekte, die auf ein
   Datefeld zur�ckgreifen m�ssen. Das Datenfeld kann sich auf dem Bildschirm
   darstellen, sein Maximum bestimmen, sich ein eine Datei schreiben, Daten
   aus einer Datei einlesen und alle St�tzpunkte auf einen bestimmten Wert
   setzen
*/

#include"datarray.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"patchio.h"
// #include<graph.h>

#define _Range_check

void TDatenfeld::free_Daten()
{
  if (eigene_Daten == 1 && Stuetzpunkte > 0 && DADaten != NULL)
  {
  	//std::cout<<"eigene_daten: "<<eigene_Daten<<" Stuetzpunkte: "<<Stuetzpunkte<<std::endl;
    free((char *)DADaten);
	  //std::cout << "free (char*)DADaten)"<< std::endl; 
  }
}

/*
  Dafuer gibt es jetzt die explizite Methode:
        TDoublefeld::load_matlab_file(char *FName);

TDatenfeld::TDatenfeld(int Handle, char *n_name)
{
  ifstream File(Handle);
  Fmatrix Matrixinfo;
  char *Name[20];

  if (!File) abort();
  File.read((char *) &Matrixinfo,sizeof(Fmatrix));
  if (!File) { std::cerr << "Fehler beim einlesen eines Matlabfiles"; abort(); };
  Stuetzpunkte = Matrixinfo.ncols;
  File.read((char *) Name,Matrixinfo.namlen);
  Datengroesse = sizeof(double);
  eigene_Daten = 1;
  DADaten = (char *)malloc((Stuetzpunkte+1) * Datengroesse);
  File.read(DADaten,Stuetzpunkte * Datengroesse);
  if (!File) { std::cerr << "Fehler nach lesen des Matlabfiles"; abort(); };
  File.close();
  name = (char*)malloc(strlen(n_name)+1);
  strcpy(name,n_name);
}
*/

TDatenfeld::TDatenfeld(TDatenfeld& N_Datenfeld, char *n_name)
{
  char name[255]; //TOBIAS	
  Stuetzpunkte = N_Datenfeld.Stuetzpunkte;
  Datengroesse = N_Datenfeld.Datengroesse;
  DADaten = N_Datenfeld.DADaten;
  eigene_Daten = 0;
//  name = (char*)malloc(strlen(n_name)+1); //TOBIAS 
//#ifdef debug_mem_datarray
//       std::cout << "1:malloc (char*)name)"<<(strlen(n_name)+1)<<std::endl;
//#endif
  strcpy(name,n_name);
}

TDatenfeld::TDatenfeld(unsigned long int N_Stuetzpunkte, int N_Datengroesse,
                       char *n_name)
{	
  char name[255];	//TOBIAS
  Stuetzpunkte = N_Stuetzpunkte;
  Datengroesse = N_Datengroesse;
  if (Stuetzpunkte)
     {
    DADaten = (char *)malloc((Stuetzpunkte+1) * Datengroesse); //TOBIAS
#ifdef debug_mem_datarray
    std::cout << "1:malloc (char*)Daten)"<<(Stuetzpunkte+1) * Datengroesse<<std::endl;
#endif
  }
  else DADaten = 0;
  eigene_Daten = 1;
//  name = (char*)malloc(strlen(n_name)+1);
//#ifdef debug_mem_datarray
//    std::cout << "2:malloc (char*)name)"<<strlen(n_name)+1<<std::endl;
//#endif
  strcpy(name,n_name);
}

TDatenfeld::~TDatenfeld()
{
  free_Daten();
}

void * TDatenfeld::nr(unsigned long int Pos)
{
  return DADaten + Datengroesse * (Pos-1);
}

void TDatenfeld::neue_Dimension(unsigned long int N_Stuetzpunkte)
{
  free_Daten();
  Stuetzpunkte = N_Stuetzpunkte;
  if (Stuetzpunkte)
  {
    DADaten = (char *)malloc((Stuetzpunkte+1) * Datengroesse);
#ifdef debug_mem_datarray
    std::cout << "2:malloc (char*)DADaten)"<<(Stuetzpunkte+1) * Datengroesse<<std::endl;
#endif
  }
  else DADaten = 0;
}
/*
 TDatenfeld::TDatenfeld(FILE *S)
{

  S.Read(Stuetzpunkte,sizeof(Stuetzpunkte));
  S.Read(Max_Daten,sizeof(Max_Daten));
  Feldgroesse = (Stuetzpunkte + 1) * Datengroesse;
  getmem(DADaten,Feldgroesse);
  S.Read(*Daten,Feldgroesse);

};

void TDatenfeld::Store(FILE *S)
{

  S.Write(Stuetzpunkte,sizeof(Stuetzpunkte));
  S.Write(Max_Daten,sizeof(Max_Daten));
  S.Write(*DADaten,Feldgroesse);

};
*/
/*
TDatenfeld::File_lesen(file Quelle)
{

  Initialisieren(filesize(Quelle));
  blockread(Quelle,*DADaten,filesize(Quelle))

};
*/
/*
void TDatenfeld::in_File_Schreiben(file Quelle)
{

 unsigned long int Vergleich, Groesse;

  blockwrite(Ziel,*DADaten,Stuetzpunkte,Vergleich);
  if Vergleich <> Stuetzpunkte then write(#7);

};
*/
/*
void TDatenfeld::Darstellen(TBildschirmausschnitt Fenster,
                            unsigned long int Startpos,
                            unsigned long int Endpos,
                            int Farbe,
                            float Faktor)

{

double X_Faktor, Y_Faktor, Mitteln;
int X, Y, Alt_X, Alt_Y;
unsigned long int R;

//                     procedure Zeichnen

Set_Max_Daten(Startpos, Endpos);
Set_Min_Daten(Startpos, Endpos);
if (Max() - Min() != 0)
  {
    Y_Faktor = Faktor * ( Fenster.UY -  Fenster.OY) /
               (Max() - Min());
    if (Endpos == Startpos) Endpos = Stuetzpunkte;
    X_Faktor = (Endpos - Startpos + 1.0) / ( Fenster.RX -  Fenster.LX);
    Y = Fenster.UY - round(Y_Faktor * D_Wert(Startpos)- Min());
    Alt_X = 0;
    if ((Y >= Fenster.OY) && (Y <= Fenster.UY))
      Alt_Y = Y;
     else if (Y < Fenster.UY)
      Alt_Y = Fenster.UY;
      else Alt_Y = Fenster.OY;
    for (X = 0;X <= (Fenster.RX - Fenster.LX - 1);X++)
    {
      Mitteln = D_Wert(Startpos + (unsigned long int)(X * X_Faktor + 0.5));
      if (X_Faktor > 1)
      {
        for (R = Startpos + (unsigned long int)(X * X_Faktor + 0.5) + 1;
           R <= Startpos + (unsigned long int)(X * X_Faktor + 0.5) +
            (unsigned long int)X_Faktor - 1; R++)
              Mitteln += D_Wert(R);
        Mitteln = Mitteln / (unsigned long int)X_Faktor;
      }
      Y = Fenster.UY - round(Y_Faktor * (Mitteln - Min()));
      if ((Y >= Fenster.OY) && (Y <= Fenster.UY))
        g_line(Alt_X + Fenster.LX, Alt_Y, X + Fenster.LX, Y, Farbe);
      Alt_X = X;
      Alt_Y = Y;
    }
  }

}
*/
void TDatenfeld::Set_Max_Daten(unsigned long int Startpos,
                               unsigned long int Endpos)
{
  unsigned long int i;

  if (Endpos == 0) Endpos = Stuetzpunkte;
  Max_Daten = Startpos;
  if (Startpos < Endpos) {
    for (i = Startpos + 1; i <= Endpos; i++)
      if (D_Wert(i) >= D_Wert(Max_Daten))
      {
        Max_Daten = i;
      }
  } else {
    for (i = Startpos - 1; i >= Endpos; i--)
      if (D_Wert(i) >= D_Wert(Max_Daten))
      {
        Max_Daten = i;
      }
  } /* endif */
}

void TDatenfeld::Set_Min_Daten(unsigned long int Startpos,
                               unsigned long int Endpos)
{
  unsigned long int i;

  if (Endpos == 0) Endpos = Stuetzpunkte;
  Min_Daten = Startpos;
  if (Startpos < Endpos) {
    for (i = Startpos + 1; i <= Endpos; i++)
      if (D_Wert(i) <= D_Wert(Min_Daten))
      {
        Min_Daten = i;
      }
  } else {
    for (i = Startpos - 1; i >= Endpos; i--)
      if (D_Wert(i) <= D_Wert(Min_Daten))
      {
        Min_Daten = i;
      }
  } /* endif */
}

int TDatenfeld::ASCII_file(const char* Filename)
{
  unsigned long int i;
  std::ofstream ASCII_file(Filename);

  if (!ASCII_file) return ASCII_file.rdstate();
  for (i = 0; i <= Stuetzpunkte;i++)
  {
//    std::cout << D_Wert(i);
    ASCII_file << D_Wert(i) << "\n";
  }
  ASCII_file.close();
  return ASCII_file.rdstate();
}

std::ostream& operator <<(std::ostream& str, TDatenfeld& Datenfeld)
{
	 std::cout << "vorher\n";
  int SP;

  SP = (int)Datenfeld.Stuetzpunkte;
  str.write((char *) &SP,4);
  str.write((char *) Datenfeld.DADaten + Datenfeld.Datengroesse,
            (int)(Datenfeld.Stuetzpunkte * Datenfeld.Datengroesse));
            std::cout << "nacher\n";
  return str;
}

std::istream& operator >>(std::istream& str, TDatenfeld& Datenfeld)
{
  long int SP, Zahl = 1;

  str.read((char *) &SP,4);
  if (SP < (int)Datenfeld.Stuetzpunkte)
  {
    std::cout <<
"\nWarnung : Es wurde versucht ein Datenfeld mit einem File zu \n" <<
"initialisieren, das kleiner als das Datenfeld ist.\n" ;
    abort();
  }
  if (Datenfeld.Stuetzpunkte == 0)
    Datenfeld.neue_Dimension((unsigned long int)SP);
  str.read(Datenfeld.DADaten + Datenfeld.Datengroesse,
            (int)(Datenfeld.Stuetzpunkte * Datenfeld.Datengroesse));
  Zahl = str.rdstate();
  std::cout << Zahl;
  if (str.gcount() != Datenfeld.Stuetzpunkte * Datenfeld.Datengroesse)
  {
    std::cerr << "Fehler beim lesen eines Datenfeldes\n" <<
            str.gcount() << " bytes gelesen\n";
    abort();
  }
  return str;
}

double& TDoublefeld::operator [] (unsigned long int Pos)
{

#if defined (_Range_check)
  if (Pos > Stuetzpunkte)
  {
    printf("Zugriff auf Element auserhalb des g�ltigen Bereiches\n");
    printf("Verlangte Position : %lu \n", Pos);
    if (name[0]) {
      printf("Name der Varibale: %s", name);
    } /* endif */
    abort();
  }
  else
#endif
  return D[Pos];
  
}

double TDoublefeld::D_Wert(unsigned long int Pos)
{
  double Test;

  Test = operator [] (Pos);
  return operator [] (Pos);
}

int TDoublefeld::Int_Wert(unsigned long int Pos)
{
  return (int)operator [] (Pos);
}

void TDoublefeld::Reset(double Wert)
{
  unsigned long int i;

  for(i = 0; i <= Stuetzpunkte; operator [] (i) = Wert, i++);
}

TDoublefeld& TDoublefeld::operator = (TDoublefeld& Feld)
{
  unsigned long int i;

  if (Feld.Stuetzpunkte != Stuetzpunkte)
  {
    std::cerr << "Doublefeld 'operator =' mit ungleich langen Feldern aufgerufen" << std::endl;
    abort();
  }
  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = Feld.operator [] (i);
  return *this;
}

TDoublefeld& TDoublefeld::operator += (double Summand)
{
  unsigned long int i;

  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = operator [] (i) + Summand;
  return *this;
}

TDoublefeld& TDoublefeld::operator *= (double Faktor)
{
  unsigned long int i;

  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = operator [] (i) * Faktor;
  return *this;
}

double TDoublefeld::Summieren(int anfang,int ende)
{
  int i;
  double S = 0;

  if (ende == 0) {
     ende = Stuetzpunkte;
  } /* endif */
  for (i=anfang; i<=ende; i++) {
     S += operator [](i);
  } /* endfor */
  return S;
}

void TDoublefeld::Funktion(double F(double))
{
  unsigned long int i;

  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = F(operator [] (i));
}


int TDoublefeld::ASCII_file_lesen(const char* Filename)
{
  unsigned long int i = 0;
  std::ifstream ASCII_file(Filename);
  double dummy;

  if (!ASCII_file) return ASCII_file.rdstate();
  while((ASCII_file >> dummy) && (++i <= Stuetzpunkte)) operator[] (i) = dummy;
#ifdef DEBUG
  if (i > Stuetzpunkte) {std::cerr << "Zu grosses ASCII File eingelesen\n";  }
#endif
  // Tobias 2017 auskommentiert ASCII_file.clear((!std::ios::eofbit) & ASCII_file.rdstate());
  ASCII_file.close();
  return ASCII_file.rdstate();
}

int TDoublefeld::load_matlab_file(char *FName)
{
  char name[256]; //TOBIAS	
  short i;	
  std::ifstream File;
  Fmatrix Matrixinfo;
  char n_name[256];
  long unsigned int ii;
if (Daten.matlab_file_loaded == false)						//Tobias 2022 to prevent multiple load from file, stepresponses are globally stored and passed to this routine at subsequent loads
   {
   	Daten.matlab_file_loaded = true;
   }
else   											//Tobias 2022 to prevent multiple loading
	 {
	  Stuetzpunkte = Daten.Stuetzpunkte_loaded;
	  //std::cout<<"Stuetzpunkte: "<<Stuetzpunkte<<std::endl;
	  DADaten = (char *)malloc((Stuetzpunkte+1) * sizeof (double));  
	  //std::cout<<"Malloc"<<std::endl;
	  D = (double *) DADaten;	
	  for (int ii = 0; ii < Stuetzpunkte; ii++)
	     {
  	   D[ii] = Daten.stepresponse_loaded[ii];
  	   //std::cout<<"i:"<<ii<<" D:"<<D[ii]<<" memory:"<<Daten.stepresponse_loaded[ii]<<std::endl;
  	   }		 
    return 0;	 
   }		 

/* #ifdef unix 
  if (FName[strcspn(FName, ".")+1] != 's')
    {
      std::cerr << "\nConvert the stepresponse file bevor using it on the " 
	   << "sun workstation with 'conv -s setfilename' .\nThe stepresponsefile "
	   << "given in the setfile must be in the PC-Format.\n"
	   << "conv will distribute a *.sun file. This file can be used as"
	   << " the stepresponsefile on the sun" << std::endl;
      exit(1);
    }   
#endif
*/  	
  find_file(FName, "DATEN", 0, 0);   
  File.open(FName, std::ios::in | std::ios::binary);


  if (!File) {
     sprintf(Daten.errormessage,"Fehler beim oeffnen des Matlabfiles1");	
     std::cout << "Fehler beim oeffnen des Matlabfiles" << std::endl;
     return 1;
  }

  File.read((char *) &Matrixinfo,sizeof(Fmatrix));
	
  if (!File) {
     sprintf(Daten.errormessage,"Fehler beim einlesen eines Matlabfiles2");		
     //std::cerr << "Fehler beim einlesen eines Matlabfiles";
     return 1;
  };
  File.read((char *) n_name,Matrixinfo.namlen);  

		
  
#ifdef debug_datarray 
std::cout << Matrixinfo.type << "\t" << Matrixinfo.mrows << "\t" << Matrixinfo.ncols << "\t" << Matrixinfo.imagf << "\t" << Matrixinfo.namlen << "\t" << n_name << std::endl;
#endif
  Stuetzpunkte = Matrixinfo.ncols;
    	 //std::cout<<"Stuetzpunkte:"<<Stuetzpunkte<<std::endl;	
  if (Stuetzpunkte == 0 ){
    sprintf(Daten.errormessage,"stepresponse file has zero response points in it");
    //std::cerr << "stepresponse file has zero response points in it" << std::endl;
    return 1;
  }
  if (Stuetzpunkte <= 5 ){
    std::cerr << "stepresponse file has " << Stuetzpunkte 
	 << " response points in it" << std::endl;
  }
#ifdef debug_datarray  
  std::cout << "Stuetzpunkte : " << Stuetzpunkte << std::endl;
  std::cout << "Datengroesse : " << Datengroesse << std::endl;
#endif
  //eigene_Daten = 1; 
  //free_Daten();
  DADaten = (char *)malloc((Stuetzpunkte+1) * Datengroesse);    
  //std::cout<<"Malloc"<<std::endl;
  //std::cout<<"FNAME: "<<FName<<" stuetz chars: "<<(Stuetzpunkte+1) * Datengroesse<<std::endl;
#ifdef debug_mem_datarray
  std::cout << "3:malloc (char*)DADaten)"<<(Stuetzpunkte+1) * Datengroesse<<std::endl;
#endif
  D = (double *) DADaten;
  File.read((char *)D,Stuetzpunkte * Datengroesse); //pinpoint Tobias2017 File.read((char *)D+1,Stuetzpunkte * Datengroesse);
  
  Daten.stepresponse_loaded = new double[(Stuetzpunkte+1)];
  Daten.Stuetzpunkte_loaded = Stuetzpunkte;
  for (int ii = 0; ii < Stuetzpunkte; ii++)
  	{
  	 Daten.stepresponse_loaded[ii] = D[ii];	
     //std::cout<<"i: "<<ii<<":"<<D[ii]<<":"<<Daten.stepresponse_loaded[ii]<<"  stp: "<<Stuetzpunkte<<" FName: "<<FName<<std::endl;
    } 
  //for (i=0; i<Stuetzpunkte; i++)
  //  std::cout<<"D"<<i<<":"<<D[i]<<std::endl;

  if (!File) {
    sprintf(Daten.errormessage,"Fehler nach lesen des Matlabfiles");	
    std::cerr << "Fehler nach lesen des Matlabfiles" << std::endl;
    return 1;
  };
  File.close();
  

   
#ifdef debug_mem_datarray
  for (unsigned int n=0; n<= Stuetzpunkte; n++) 
  		std::cout << "Stuetzpunkt " << n << " : " << D[n] << std::endl;
#endif

  if (name[0] == 0) { 
//    name = (char*)malloc(strlen(n_name)+1); //TOBIAS
#ifdef debug_mem_datarray
//    std::cout << "3:malloc (char*)name)"<<strlen(n_name)+1<<std::endl;
#endif
    strcpy(name,n_name);
  }  //endif 
  return 0;
}

float& TFloatfeld::operator [] (unsigned long int Pos)
{
#if defined (_Range_check)
  if (Pos > Stuetzpunkte)
  {
    printf("Zugriff auf Element auserhalb des g�ltigen Bereiches\n");
    printf("Verlangte Position : %lu \n", Pos);
    if (name[0]) {
      printf("Name der Varibale: %s", name);
    } /* endif */
    abort();
  }
  else
#endif
  return F[Pos];
}

double TFloatfeld::D_Wert(unsigned long int Pos)
{
  double Test;

  Test = operator [] (Pos);
  return Test;
}

int TFloatfeld::Int_Wert(unsigned long int Pos)
{
  return (int)operator [] (Pos);
}

void TFloatfeld::Reset(float Wert)
{
  unsigned long int i;

  for(i = 0; i <= Stuetzpunkte; operator [] (i) = Wert, i++);
}

TFloatfeld& TFloatfeld::operator = (TFloatfeld& Feld)
{
  unsigned long int i;

  if (Feld.Stuetzpunkte != Stuetzpunkte)
  {
    std::cerr << "Doublefeld 'operator =' mit ungleich langen Feldern aufgerufen" << std::endl;
    abort();
  }
  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = Feld.operator [] (i);
  return *this;
}

TFloatfeld& TFloatfeld::operator += (float Summand)
{
  unsigned long int i;

  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = operator [] (i) + Summand;
  return *this;
}

TFloatfeld& TFloatfeld::operator *= (float Faktor)
{
  unsigned long int i;

  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = operator [] (i) * Faktor;
  return *this;
}

void TFloatfeld::Funktion(float F(float))
{
  unsigned long int i;

  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = F(operator [] (i));
}


int TFloatfeld::ASCII_file_lesen(const char* Filename)
{
  unsigned long int i = 0;
  std::ifstream ASCII_file(Filename);
  float dummy;

  if (!ASCII_file) return ASCII_file.rdstate();
  while((ASCII_file >> dummy) && (++i <= Stuetzpunkte)) operator[] (i) = dummy;
  if (i > Stuetzpunkte) {std::cerr << "Zu grosses ASCII File eingelesen\n";  }
  //Tobias 2017 auskommentiert ASCII_file.clear((!ios::eofbit) & ASCII_file.rdstate());
  ASCII_file.close();
  return ASCII_file.rdstate();
}

int& TIntegerfeld::operator [] (unsigned long int Pos)
{
#if defined (_Range_check)
  if (Pos > Stuetzpunkte)
  {
    printf("Zugriff auf Element auserhalb des g�ltigen Bereiches \n");
    printf("Verlangte Position : %lu \n", Pos);
    if (name[0]) {
      printf("Name der Varibale: %s", name);
    } /* endif */
    abort();
  }
  else
#endif
  return I[Pos];
  return I[Pos];
}

double TIntegerfeld::D_Wert(unsigned long int Pos)
{
  return (double)operator [] (Pos);
}

int TIntegerfeld::Int_Wert(unsigned long int Pos)
{
  return operator [] (Pos);
}

void TIntegerfeld::Reset(int Wert)
{
  unsigned long int i;

  for(i = 0; i <= Stuetzpunkte; operator [] (i) = Wert, i++);
}

short int& TShortIntfeld::operator [] (unsigned long int Pos)
{
#if defined (_Range_check)
  if (Pos > Stuetzpunkte)
  {
    printf("Zugriff auf Element auserhalb des g�ltigen Bereiches \n");
    printf("Verlangte Position : %lu \n", Pos);
    if (name[0]) {
      printf("Name der Varibale: %s", name);
    } /* endif */
    abort();
  }
  else
#endif
  return SI[Pos];
  return SI[Pos];
}

double TShortIntfeld::D_Wert(unsigned long int Pos)
{
  return (double)operator [] (Pos);
}

int TShortIntfeld::Int_Wert(unsigned long int Pos)
{
  return (int)operator [] (Pos);
}

void TShortIntfeld::Reset(int Wert)
{
  unsigned long int i;

  for(i = 0; i <= Stuetzpunkte; operator [] (i) = (short int)Wert, i++);
}

void TShortIntfeld::operator *= (double Faktor)
{
  unsigned long int i;
  short int dummy;

  for (i = 0; i <= Stuetzpunkte; i++)
    if ((dummy=operator[] (i)) > 0) {
        operator [] (i) = (short int)((operator [] (i) * Faktor) + .5);
    } else {
        operator [] (i) = (short int)((operator [] (i) * Faktor) - .5);
    } /* endif */
}

void TShortIntfeld::operator += (short int Summand)
{
  unsigned long int i;

  for (i = 0; i <= Stuetzpunkte; i++)
    operator [] (i) = operator [] (i) + Summand;
}

