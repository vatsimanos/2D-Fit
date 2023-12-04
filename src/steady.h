/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#ifndef	_STEADY_STATE_H_
#define	_STEADY_STATE_H_

// #include "glib.h"

extern double l_steady[7];
extern double p_steady[7][7];

void single_steady(double k12, double k21, double k23, double k32, double &c_steady, double &g_steady, double &o_steady);
void multi_steady(double k12, double k21, double k23, double k32, short ch);
double pot(double x, short k);
int bin(short n, short k);

#endif /* _STEADY_STATE_H_ */
