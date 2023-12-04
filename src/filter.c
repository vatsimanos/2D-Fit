#include "filter.h"

TFilter::TFilter(int Ordnung, double* ZP, double* NP, double Startwert)
  : Zaehlerp(ZP, Ordnung,0), Nennerp(NP, Ordnung,0),
    Eingang(Ordnung,Startwert,0), Ausgang(Ordnung,Startwert,0)
{}

double TFilter::Out(double In)
{
  double Ergebnis;

  Eingang >>= 1;
  Ausgang >>= 1;
  Eingang(0) = In;
//  Ergebnis = Eingang * Zaehlerp;
  Ergebnis = (Eingang * Zaehlerp) - (Ausgang * Nennerp);
  if (Nennerp(0) != 1) Ergebnis = Ergebnis / Nennerp(0);
  Ausgang(0) = Ergebnis;
  return Ergebnis;
}
