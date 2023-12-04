#include "channel.h"

TKanal::TKanal(int N_Zustaende,double  **N_Stromvektor, double **N_Matrix,
               //double N_Sigma, int Startzustand) felixw 17.11.99
	        double  **N_initial_state_value, double N_Sigma)
   : Ratenmatrix(N_Matrix, N_Zustaende, N_Zustaende),
     Uebergangsmatrix(N_Zustaende, N_Zustaende),
     Stromvektor(N_Stromvektor, N_Zustaende),
     Bes_Wahrsch(N_Zustaende),
     // eingefuegt fw 17.11.
     initial_state_value(N_initial_state_value, N_Zustaende),
     Ratenvektor(N_Zustaende)
{
  int i, j;

  Zustaende = N_Zustaende;
  Zustand = alter_Zustand = Ziel = Startzustand(N_initial_state_value);
  Verweildauer = 0;
  t = 0;
  // Initialisierung des in chanel.h beschriebenen Ratenvektors
  for (i=1; i<=Zustaende; i++)
    for (j=1; j<=Zustaende; j++)
      if (i!=j) Ratenvektor(i) += Ratenmatrix(i,j);
  Uebergangsmatrix = 0;
  // der erste punkt wird auch als sprung benutzt (hoehe 0),
  // soll aber in der Matrix nicht auftauchen
  // warum weiss ich nicht
  Uebergangsmatrix(Zustand,Zustand) = -1;
  Sigma = N_Sigma;
  k = 1;
  jumps=new TTSprung;  //Tobias
  jumps->naechster = NULL;
  //tmp_jump=jumps;
}

TKanal::~TKanal() 
{  	
  TTSprung *akt,*next;    //Tobias
  int i=0;
if (jumps !=NULL)
  {  
   akt=jumps;
   do{
      i++;
     next=akt->naechster;
     delete(akt);
     akt=next;
   }while(next!=NULL);
  } 
}

int TKanal::Startzustand(double  **N_initial_state_value)
{
  double Wuerfel,counter;
  int i,Ende=0; 

  counter=0;

  Wuerfel = Daten.uni_01(Daten.random_engine); // Tobias 
  i=0;
  do
    {
      i++;
      counter+=N_initial_state_value[1][i];
      if (counter>Wuerfel)
	Ende=1;
      if (i>=Zustaende)
	Ende=1;
    }
  while(!Ende);
  return 1; //Unterbinden des zuf�lligen Anfangszustandes
}


void TKanal::naechster_Zustand()
{
  double Wuerfel;

  alter_Zustand = Zustand;
  Zustand = Ziel;
  Ziel = 0;
  Verweildauer = 0;


  Wuerfel = - Ratenvektor(Zustand) * Daten.uni_01(Daten.random_engine); //TOBIAS 
  do {
    Ziel += 1;
    if (Zustand!=Ziel) Wuerfel += Ratenmatrix(Zustand,Ziel);
  } while (Wuerfel < 0);
  // alle im folgenden veraenderten Variablen gehoeren zu der speziellen
  // Instanz von TKanal und ind deshalb beim naechsten aufruf von z.B.
  // der Fkt Strom noch erhalten insbesondere gilt dies fuer
  // Ziel und t
  ++Uebergangsmatrix(alter_Zustand, Zustand);
  // log ist der nat. Logarithmus also ln, d.h.
  // 0 < -1 * log(ran1(&k)) < infinity
  // die Wahrscheinlichkeit fuer eine Verweildauer t in einem Zustand ist 
  // Proportional zu exp(-t) 
  // also liefert ran1(&k) die Wahrscheinlichkeit der Verweildauer
  Verweildauer = -1 * log(Daten.uni_01(Daten.random_engine)+1e-17) / Ratenvektor(Zustand);  //Tobias Verweildauer = -1 * log(ran1(&k)) / Ratenvektor(Zustand);//1e-17 to avoid zero
  // 
  Bes_Wahrsch(Zustand) += Verweildauer;

  t += Verweildauer;  
}

unsigned short int TKanal::Strom(double Zeit)
{	
  // waehle den naechsten Zustand aus und bestimme wann der Uebergang ist,
  // wenn dies noch nicht (beim letztenmal) bestimmt worden ist
  if (Zeit<t)
    ++Uebergangsmatrix(Zustand, Zustand);
  while (Zeit > t)
    naechster_Zustand();
  return (unsigned short int) floor (0.5 + Stromvektor(Zustand)); //Tobias rint durch floor + 0.5 ersetzt und if weg

}

double TKanal::Bes_W(short int Niveau)
{	
  int i;
  double Bes_Wahr = 0;

  for (i = 1; i <= Zustaende; i ++)
   if ((short int)Stromvektor(i) == Niveau) Bes_Wahr += Bes_Wahrsch(i);
  return Bes_Wahr / t;
}

double TKanal::Bes_W(int Zustand)  { return Bes_Wahrsch(Zustand) / t;  }

Matrix TKanal::Uebergaenge() {return Uebergangsmatrix; };

TKanal_gefiltert::TKanal_gefiltert(int N_Zustaende, double **N_Stromvektor,
        double **N_Matrix, double  **N_initial_state_value, double f_3_dB, 
				   double N_Sigma) :
 TKanal(N_Zustaende, N_Stromvektor, N_Matrix, N_initial_state_value, N_Sigma)
{
  tau = 1 / (f_3_dB * 2 * M_PI);
  Ausgang = Stromvektor(Zustand);
}

unsigned short int TKanal_gefiltert::Strom(double Zeit)
{	
  static double Ergebnis;

  if (Zeit<t)
    ++Uebergangsmatrix(Zustand, Zustand);
  while (Zeit > t)
  {
    Ausgang = Stromvektor(Zustand) -
             (Stromvektor(Zustand) - Ausgang) * exp(-Verweildauer / tau);
    naechster_Zustand();
  }
  Ergebnis = Stromvektor(Zustand) - (Stromvektor(Zustand) - Ausgang) *
             exp((t - Verweildauer - Zeit) / tau);
  
  return (unsigned short int) floor (0.5 + Ergebnis);  //Tobias rint durch floor + 0.5 ersetzt und if weg

}

TKanal_digfi::TKanal_digfi(int N_Zustaende, double **N_Stromvektor,
        double **N_Matrix, double N_dT_Sprung,
        TDoublefeld& N_Sprungantwort, double  **N_initial_state_value,
			   double N_Sigma) :
 TKanal(N_Zustaende, N_Stromvektor, N_Matrix, N_initial_state_value, N_Sigma),
 Sprungantwort(N_Sprungantwort),Differenzenmatrix(N_Zustaende, N_Zustaende)	
{	
  int i, j;

  Spruenge = new TTSprung;   //Tobias
  Spruenge->Zeit = 0;
  Spruenge->Differenz = 0;
  Spruenge->naechster = NULL;
  letzter_Sprung = Spruenge;
  No_Zeiten = 0;
  Summe = Stromvektor(Zustand);
  // Sprungantwort *= 0.1; // was generated by matlab in a step with hight 10
  dT_Sprung = N_dT_Sprung;
  for(i = 1; i <= Zustaende; i++)
   for(j = 1; j <= Zustaende; j++)
    Differenzenmatrix(i,j) = (Stromvektor(j) - Stromvektor(i)); 
};

TKanal_digfi::~TKanal_digfi()
{ 	
 
  TTSprung *akt,*next;   //Tobias
  
  akt=Spruenge;
  do{
    next=akt->naechster;
    delete(akt);
    akt=next;
    //cerr << "hier werden Spruenge geplaettet\n";
  }while(next!=NULL); 
   
}

void TKanal_digfi::clear()
{			
  TTSprung *akt,*next;   //Tobias
  
  akt=Spruenge;
  do{
    next=akt->naechster;
    delete(akt);
    akt=next;
    //cerr << "hier werden Spruenge geplaettet\n";
  }while(next!=NULL);
}

unsigned short int TKanal_digfi::Strom(double Zeit)
{	
  static unsigned long int lDt;
  static double Ergebnis, Beitrag, dDt;
  static TTSprung *neuer_Sprung, *aktueller_Sprung;   //Tobias
  static char Ende;

  if (Zeit<t)
    ++Uebergangsmatrix(Zustand, Zustand);
  while (Zeit > t)
  {
    alter_Zustand = Zustand;
    naechster_Zustand();

    Summe += Differenzenmatrix(alter_Zustand,Zustand);
    /*cout << "Summe " << Summe << "\tZeit " << Zeit 
	 << "\talter_Zustand " << alter_Zustand
	 << "\tZustand " << Zustand << endl;*/
    if (letzter_Sprung->naechster == NULL)
     {
      neuer_Sprung = new TTSprung;   //Tobias
     } 
    else
    {
      neuer_Sprung = letzter_Sprung->naechster;
      letzter_Sprung->naechster = neuer_Sprung->naechster;
    }
    neuer_Sprung->Zeit = t;
    neuer_Sprung->Differenz = Differenzenmatrix(Zustand, Ziel);
    neuer_Sprung->naechster = Spruenge;
    Spruenge = neuer_Sprung;
  }
    
  aktueller_Sprung = Spruenge->naechster;

  Ergebnis = Summe;
  Ende = 0;
  do
  {
    Beitrag = 0;
    dDt = (Zeit - aktueller_Sprung->Zeit) / dT_Sprung + 1;
    lDt = (unsigned long int) (dDt);
    	
    if ((lDt + 1) < Sprungantwort.Stuetzpunkte)
    {
        //std::cout<<"Sprungantwort.Stuetzpunkte:"<< Sprungantwort.Stuetzpunkte <<std::endl;

        Beitrag = 1 - Sprungantwort[lDt]* 0.1; //Efthymios 2023 old: Beitrag = 1 - (Sprungantwort[lDt+1]*(dDt-lDt)+ Sprungantwort[lDt] * (lDt + 1 - dDt)) * 0.1; 
        //adapted since  recorded step response is sampled at 100 kHz and there is no need to scale it since simulations are also at 100 kHz. 
        //Therefore: needs to be adapted if sampling frequency deviates from 100 kHz             
            
      Ergebnis -= Beitrag * aktueller_Sprung->Differenz;


      if (aktueller_Sprung->naechster != NULL)
        aktueller_Sprung = aktueller_Sprung->naechster;
      else Ende = 1;
    }
    else Ende = 1;
  } while (!Ende);
  letzter_Sprung = aktueller_Sprung;

  return (unsigned short int) floor (0.5 + Ergebnis);  //Tobias rint durch floor + 0.5 ersetzt und if weg

}

TKanal_mehrfach_gefiltert::TKanal_mehrfach_gefiltert(int N_Zustaende,
   double **N_Stromvektor, double **N_Matrix,
     Vektor& f_3_dB, double  **N_initial_state_value, double N_Sigma)
 : TKanal(N_Zustaende, N_Stromvektor, N_Matrix, N_initial_state_value,N_Sigma),
   tau(f_3_dB.E()), Ausgang(f_3_dB.E())
{	
  int i;

  N_Taus = f_3_dB.E();
  Ergebnis = (unsigned short int *)malloc(N_Taus * sizeof(short int));
  for (i = 1; i <= N_Taus; i++)
  {
    tau(i) = 1 / (f_3_dB(i) * 2 * M_PI);
    Ausgang(i) = Stromvektor(Zustand);
  }
}

TKanal_mehrfach_gefiltert::~TKanal_mehrfach_gefiltert()
{
  free(Ergebnis);
}

unsigned short int * TKanal_mehrfach_gefiltert::Strom(double Zeit)
{	
  int I;
  static int i;
  static double Rauschen;

  while (Zeit > t)
  {
    for (i = 1; i <= N_Taus; i++)
      Ausgang(i) = Stromvektor(Zustand) -
             (Stromvektor(Zustand) - Ausgang(i)) * exp(-Verweildauer / tau(i));
    naechster_Zustand();
  }
  if (Sigma==0) Rauschen = 0;
  else Rauschen = Daten.snd_01(Daten.random_engine)*Sigma;		//gasdev (nr) ersetzt
  for (i = 1; i <= N_Taus; i++)
    Ergebnis[i] = (unsigned short int) floor (0.5 + Stromvektor(Zustand)-(Stromvektor(Zustand)- //tobias rint durch floor + 0.5 ersetzt
                  Ausgang(i))* exp((t-Verweildauer-Zeit)/tau(i))+Rauschen); 
  return Ergebnis;
}
