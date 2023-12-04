/*******************************************************************
 *  Patch Clamp
 *  1999 Oliver Radomski
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#ifndef	_ERROR_H_
#define	_ERROR_H_

// diese funktion entscheidet, ob ein fenster offen ist
// dann fehlermeldung in app, sonst g_warning
void warning (char* text);

// return val fuer moegliche user entscheidung in app
// (OK:0 oder Cancel:-1)
short request(char* text);

#endif /* _ERROR_H_ */
