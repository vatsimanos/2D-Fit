/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#ifndef _GTKAPP_H_
#define _GTKAPP_H_


class TApplication
{
private:
  

  int      onoise_id, showmark_id, bw_id, ls_id; 

public:

  short     xMouse1;
  short     xMouse2;
  short     xMouse2old;
  short     yMouse1;
  short     yMouse2;
  short     yMouse2old;
  short     i_window, j_window;
  bool   BtnDownFirst;

  bool autom_steuer; // automatische steuerung
 
  TApplication();

  void set_gcscolors();
  void init_gcs();
  void set_menu_state(bool unblock = true);

  void cleardrawingarea();

  void statustext(const char *text);
};

extern TApplication Application;

#endif /* _GTKAPP_H_ */
