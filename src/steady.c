/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

//#include "glib.h"

#include "steady.h"

double l_steady[7];
double p_steady[7][7];

void 
single_steady(double k12, double k21, double k23, double k32, double &c_steady, double &g_steady, double &o_steady)
{
double  nenner = k12 * k32 + k12 * k23 + k32 * k21;

  c_steady = k32 * k21 / nenner;
  g_steady = k12 * k23 / nenner;
  o_steady = k12 * k32 / nenner;
}

void 
multi_steady(double k12, double k21, double k23, double k32, short ch)
{
short  k,l;
double c_steady, g_steady, o_steady;

  single_steady(k12, k21, k23, k32, c_steady, g_steady, o_steady);
  //writeln('c_steady: ', c_steady);
  //writeln('g_steady: ', g_steady);
  //writeln('o_steady: ', o_steady);
  for (k=0;k <= ch;k++)
  {
    for (l=0;l <= k;l++)
	{
      p_steady[k][l] = bin(ch,k) * bin(k,l) * pot(c_steady,l) * pot(g_steady,k-l) * pot(o_steady,ch-k);
      // writeln('multi_steady: ', p_steady[k,l]);
	}
    l_steady[k] = 0;
    for (l=0;l <= k;l++)
	{
      l_steady[k] = l_steady[k] + p_steady[k][l];
	}
  }
}

double 
pot(double x, short k)
{
  if (k == 0) 
  { 
	return 1;
  } else {
    if (k < 0) 
	{
	  return 1 / pot(x,-k);
	} else {
      return pot(x,k-1) * x;
	}
  }
}

int 
bin(short n, short k)
{
double  prod=1;
short i;

  if ((k == 0) || (k==n)) 
  {
	return 1;
  } else {
    for (i=0;i <= k-1;i++)
    {
      prod = prod * (n-i) / (i+1);
    }
    return (int)prod;
  }
}
