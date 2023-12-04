/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

//hier werden aus den steady-state bedingungen sowie Eigenvektoren und Eigenwerten
//der Ueberganfgsmatrix des absorbierenden Modells die Dwell-Time Histogramme berechnet

#ifndef	_Y_SHUT_H_
#define	_Y_SHUT_H_

#include "declare.h" 

extern row_vec row;

void calc_k_weg(short level,Tmama mama, Vector k_weg);
double calc_f_trace(short level, short channels, Tmama mama);
double calc_c_amplitude(Vector e_vec, short level, double f_trace, Vector k_weg);
void calc_c_histo(Vector lam_vec, Vector amp_vec, short level, short channels, Tcmatrix cmatrix, Tmama mama);
void check (short dim, Tmatrix A, Vector EV[], Vector EW);
void eigensystem(short rdim, Tmatrix matti, Vector eig_vec[], Vector e_val_r);
double c_strich(Vector dum_vec, short level);

#endif /* _Y_SHUT_H_ */