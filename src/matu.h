/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/
#ifndef	_MATU_H_
#define	_MATU_H_

#include "declare.h"

const short matu_u_dim = max_dim_unter_matrix; // dim auch in fit.pas

typedef double MatD[matu_u_dim+1][matu_u_dim+1];
typedef double VecD[matu_u_dim+1];
typedef short  VecID[matu_u_dim+1];
typedef VecD    VECn;

void selmhes(short n, short low, short igh, MatD a, VecID Int);
void selmbak(short low, short igh, short m, MatD a, VecID Int, MatD z);
void sbalanc(short n, short &low, short &igh, MatD a, VecD scale);
void sbalbak(short n, short low, short igh, short m, VecID scale, MatD z);
void shqr2(short n, short low, short igh, MatD h, VECn wr, VECn wi, MatD z, short &ierr);

#endif /* _MATU_H_ */

