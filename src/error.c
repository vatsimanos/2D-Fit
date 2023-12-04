/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <iostream>

//#include "glib.h"
//#include "gdk/gdk.h"

void
warning (char *text)
{
  if (0)  // mit GUI?
  {

  } else {
	//cerr
	std::cout << "*** Warning *** " << text << std::endl;
//2019    gdk_beep(); 
    return;
  }
}

short
request(const char *text)
{
  if (0)  // mit GUI?
  {

  } else {
    //g_warning("%s\n",text);
//2019    gdk_beep(); 
    return 0;
  }
}

