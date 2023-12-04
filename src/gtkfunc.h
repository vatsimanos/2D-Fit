/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#ifndef	_GTKFUNC_H_
#define	_GTKFUNC_H_

#include "daten.h"

// erstellt gewuenschten fileselektor
//void gtk_file_sel(char* titel, char* file, gpointer OK_func);		//gpointer OK_func

// diese Funktionen reichen die filenamen vom fileselector weiter
//void openjob_OK (GtkWidget *widget);
//void savejob_OK (GtkWidget *widget);
//void opentime_OK (GtkWidget *widget);
//BEGIN KARSTENASCI
//void import_OK (GtkWidget *widget);
//END KARSTENASCI
//void opensettings_OK (GtkWidget *widget);
//void savesettings_OK (GtkWidget *widget);
//void savetmatrix_OK (GtkWidget *widget);
//void savedlevel_OK (GtkWidget *widget);
//void savedthisto_OK (GtkWidget *widget);
//void exportbmp_OK (GtkWidget *widget);
//void savegnu_OK (GtkWidget *widget);
//void opendwelllb_OK(GtkWidget* widget);
//void savedwelllb_OK(GtkWidget* widget);

// macht aus RGB werten GdkColor 
/*inline GdkColor 
RGB2Gdk(unsigned short R, unsigned short G, unsigned short B)
{
GdkColor color;

  color.pixel = 0;
  color.red   = R * 256;
  color.green = G * 256;
  color.blue  = B * 256;

  return color;
}*/

// erzwingt bildaufbau
/*inline void
refresh()
{
  if (!Application.autom_steuer)
  {
    while (gtk_events_pending()) gtk_main_iteration();
  }
} */

// style draws (z.b. draw_shadow) koennen auf diese leere funktion umgebogen werden
inline void empty(...){};

// aendert style fuer alle widgets eines containers
//void 
//SetStyleRecursively (GtkWidget *widget, gpointer data);

// kopiert styleclass (in style ex. pointer auf diese struct)
//GtkStyleClass*	
//GtkStyleClassCopy (GtkStyleClass* SCTo, const GtkStyleClass* SCFrom);

#endif /* _GTKFUNC_H_ */
