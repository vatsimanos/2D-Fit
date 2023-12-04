/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <stdio.h>
#include "gtkfunc.h"
#include "gtkapp.h"
#include "daten.h"
#include "error.h"


inline void 
store_path(char* FileName)
{
  // nur dir, ohne filename sichern
  strcpy(Daten.cdir, FileName);
  char* p = strrchr(Daten.cdir,'\\');
  if (p != NULL)
  {
    p[1] = 0; 
  } else {
	Daten.cdir[0]=0; 
  }
}


