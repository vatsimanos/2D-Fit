#include <stdio.h> 

#include "gtkapp.h"
#include "declare.h"
#include "daten.h"
#include "gtkfunc.h"


// global damit auf diesen style auch noch zugeriffen werden kann,
// wenn window schon erzeugt ist.
// wird nur gebraucht, um statusbar ohne rahmen zu zeichnen!


TApplication::TApplication()
{

  autom_steuer = false;

  BtnDownFirst = false;

}

void 
TApplication::cleardrawingarea()
{
  Daten.graphtype = nothing;
  Daten.paint();
}

void 
TApplication::statustext(const char *text)
{
}


void
TApplication::set_menu_state(bool unblock) //Tobias weg: = true
{
}

void
TApplication::set_gcscolors()
{
}

void 
TApplication::init_gcs()
{ 
}
