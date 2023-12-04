/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#ifndef	_ROUND_H_
#define	_ROUND_H_

//funktion zum runden eines doubles auf den naechsten integer
inline int my_round(double zahl)		//Tobias
{
  if (zahl<0) 
  {
     return (int)(zahl - 0.5);
  } else {
     return (int)(zahl + 0.5);
  }
}

#endif /* _ROUND_H_ */
