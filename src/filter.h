#include "matrix2.h"
#include "daten.h"

struct TFilter
{	
  Vektor Zaehlerp, Nennerp, Eingang, Ausgang;
  TFilter(int Ordnung, double* ZP, double* NP, double Startwert = 0);
  double Out(double In);
};
