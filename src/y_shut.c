/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <iostream>
#include <stdio.h>
#include <cmath>

#include "declare.h"
#include "steady.h"
#include "y_shut.h"
#include "matu.h"

row_vec row;

void 
ShowMatrix(Tmatrix A, char* text)
{
short  i,j;

  std::cout << text << std::endl;;
  for (i=0;i <= TargetFitMaxChannels;i++)
  {
    for (j=0;j <= TargetFitMaxChannels;j++)
	{
      std::cout << A[i][j] << "   ";
	}
    std::cout << std::endl;
  }
}

void
ShowVector(Vector v, char* text)
{
short  i;

  std::cout << text << std::endl;
  for (i=0;i <= TargetFitMaxChannels;i++)
  {
    std::cout << v[i] << "   ";
  }
  std::cout << std::endl;
}

void
ShowDimension(short dim)
{
short dummy;

  std::cout << "Dim: " << dim << std::endl;
  std::cin >> dummy;
}

void 
calc_k_weg(short level,Tmama mama, Vector k_weg)
{
//k_weg ist die summe ueber alle Ratenkonstanten die von einem Zustand weg-fuehren
short i,j;

  for (i=0;i <= level;i++)
  {
    k_weg[i] = 0;
    for (j=0;j <= level-1;j++)
	{
      k_weg[i] += mama[row[level]+i][row[level-1]+j] + mama[row[level]+i][row[level+1]+j];
	}
	for (j=level;j <= level+1;j++)
	{
      k_weg[i] += mama[row[level]+i][row[level+1]+j];
	}
  }
}

double 
calc_f_trace(short level, short channels, Tmama mama)
{
short  i,j;
double f_trace;

  f_trace = 0;
  for (j=0;j <= level;j++)
  {
    if (level > 0)
	{
      for (i=0;i <= level-1;i++)
	  {
        f_trace += p_steady[level][j] * mama[row[level]+j][row[level-1]+i];
	  }
	}
    if (level < channels)
	{
      for (i=0;i <= level+1;i++)
	  {
        f_trace += p_steady[level][j] * mama[row[level]+j][row[level+1]+i];
	  }
	}
  }
  f_trace = 2 * f_trace;
  return f_trace;
}

double 
calc_c_amplitude(Vector e_vec, short level, double f_trace, Vector k_weg)
{
short  i;
double dum_amp;

  dum_amp = 0;
  for (i=0;i <= level;i++)
  {
    dum_amp += p_steady[level][i] * e_vec[i] * k_weg[i];
  }
  return dum_amp * dum_amp * 2 * c_strich(e_vec, level) / f_trace;
}

void 
calc_c_histo(Vector lam_vec, Vector amp_vec, short level, short channels, Tcmatrix cmatrix, Tmama mama)
{
short  k;
Vector  e_vektor[6];
Vector  k_weg;
double f_trace;

  eigensystem(level+1, cmatrix, e_vektor, lam_vec);
  calc_k_weg(level, mama, k_weg);
  f_trace = calc_f_trace(level, channels, mama);

  for (k=0;k <= level;k++)
  {
    amp_vec[k] = calc_c_amplitude(e_vektor[k], level, f_trace, k_weg);
  }
}

void check (short dim, Tmatrix A, Vector EV[], Vector EW)
{
short i;
char  buffer[MaxTextLen];

  ShowMatrix(A, (char *)"Matrix: ");
  ShowVector(EW, (char *)"Eigenwerte: ");
  for (i=0;i <= dim-1;i++)
  {
      sprintf(buffer,(char *)"%i-te Komponente der Eigenvektoren",i);
      ShowVector(EV[i], buffer);
  }
  ShowDimension(dim);
}


void 
eigensystem(short rdim, Tmatrix matti, Vector eig_vec[], Vector e_val_r)
{
double  esum;
short   a,b;
MatD     mutti, evecs;
VecID    selm;
VECn     evalr, evali;
short   ierr;

  for (a=1;a <= rdim;a++) //Matrix bereitstellen
  {
    for (b=1;b <= rdim;b++)
	{
      mutti[a][b] = matti[a-1][b-1];
	}
  }
  selmhes(rdim, 1, rdim, mutti, selm); // Transformation auf obere Hessenbergform

  //Mit Einheitsmatrix initialisieren
  for (a=1;a <= rdim;a++)
  {
    for (b=1;b <= rdim;b++)
	{
      evecs[a][b] = 0;
	}
    evecs[a][a] = 1;
  }

  shqr2(rdim, 1, rdim, mutti, evalr, evali, evecs, ierr); // Berechnung der EW + EV
  selmbak(1, rdim, rdim, mutti, selm, evecs); // R�cktransformation
  for (a=1;a <= rdim;a++) // Ergebnis sichern
  {
	esum=0;
    e_val_r[a-1] = evalr[a];
    for (b=1;b <= rdim;b++)
    {
      eig_vec[a-1][b-1] = evecs[b][a];
      esum += evecs[b][a] * evecs[b][a];
	}
    for (b=1;b <= rdim;b++)
    {
      eig_vec[a-1][b-1] = eig_vec[a-1][b-1] / sqrt(esum);
	}
  }
  // Check (rdim,matti,eig_vec,e_val_r); //PASCAL kommentar
}

double 
c_strich(Vector dum_vec, short level)
{
short  i;
double dum_c;

  dum_c = 0;
  for (i=0;i <= level;i++)
  {
    dum_c += dum_vec[i] * dum_vec[i] * p_steady[level][i];
  }
  return 1 / dum_c;
}
