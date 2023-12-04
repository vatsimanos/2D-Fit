#include "err.h"

Terror::Terror()
{
}

Terror::~Terror()
{}


char *Terror::get_errorstring(int ernr, char *location)
{
  switch (ernr){
  case OUT_OF_MEMORY :
    sprintf(errorstring,"%s : out of memory",location); break;
  case FILE_OPEN_ERROR :
    sprintf(errorstring,"%s : file open error",location); break;
  case DATA_FILE_LENGTH_WARNING :
    sprintf(errorstring,"%s : unexpected filelength",location); break;
  case EMPTY_VARIABLE :
    sprintf(errorstring,"%s : empty variable",location); break;
  case ALREADY_PRESENT :
    sprintf(errorstring,"%s : data already loaded",location); break;    
  case NO_FILENAME :
    sprintf(errorstring,"%s : no datafile specified",location); break;
  case FILE_READ_ERROR :
    sprintf(errorstring,"%s : error while reading",location); break;
  case ENVVAR_NOT_SET :
    sprintf(errorstring,"%s : environment variable not set",location); break;
  case NOT_IMPLEMENTED :
    sprintf(errorstring,"%s : not implemented yet",location); break;
  case ZBRENT:
    sprintf(errorstring,"%s : numerical recipies error",location); break;
  case NOT_CALLED_WITH_NULL :
    sprintf(errorstring,"%s : function 'char *get_ (char *)' must be called "
	    "wih  NULL Pointer",location); break;
  case NOT_PRESENT :
    sprintf(errorstring,"%s : empty variable ",location); break;
  case SYNTAX_ERROR :
    sprintf(errorstring,"%s : syntax error while reading ",location); break;
  case WRONG_CONTENT :
    sprintf(errorstring,"%s : tried to set variable with illegal content",
	    location); break;
  case OUT_OF_RANGE :
    sprintf(errorstring,"%s : illegal index",location); break;
  case READ_ERROR :
    sprintf(errorstring,"%s : error while reading File",location); break;
  case NO_VALUES :
    sprintf(errorstring,"%s : array contains no values",location); break;
  case STUD_WRONG :
    sprintf(errorstring,"%s :",location); break;
  case VALUE_NOT_MASKED :
    sprintf(errorstring,"%s : wrong content of timeseries if flags are present",location); break;
    /*
  case  :
    sprintf(errorstring,"%s :",location); break;
    */
  default:
    sprintf(errorstring,"%s : unknown error code %i",location,ernr); break;
  }
  strcat(errorstring,"\n");
  return errorstring;
}


char *Terror::get_errorstring(int ernr, char *location,char *arg)
{
  strcpy(errorstring,get_errorstring(ernr,location));
  strcat(errorstring,"\t");
  switch (ernr){
  case OUT_OF_MEMORY : 
    strcat(errorstring," while allocating memory for "); break;
  case FILE_OPEN_ERROR : strcat(errorstring," while opening "); break;
  case DATA_FILE_LENGTH_WARNING : strcat(errorstring," of "); break;
  case EMPTY_VARIABLE : strcat(errorstring," "); break;
  case ALREADY_PRESENT : strcat(errorstring,"  "); break;
  case NO_FILENAME : 
    strcat(errorstring," try loading a setfile first "); break;
  case FILE_READ_ERROR: strcat(errorstring,"  "); break;
  case ENVVAR_NOT_SET: strcat(errorstring,"  "); break;
  case NOT_IMPLEMENTED : strcat(errorstring,"  "); break;
  case ZBRENT : strcat(errorstring,"  "); break;
  case NOT_CALLED_WITH_NULL : strcat(errorstring,"  "); break;
  case NOT_PRESENT : strcat(errorstring,"  "); break;
  case SYNTAX_ERROR : strcat(errorstring," "); break;
  case WRONG_CONTENT : strcat(errorstring,"  "); break;
  case OUT_OF_RANGE : strcat(errorstring,"  "); break;
  case READ_ERROR : strcat(errorstring,"  "); break;  
  case VALUE_NOT_MASKED : strcat(errorstring,"  "); break;
  case NO_VALUES : strcat(errorstring,"  "); break;
  case STUD_WRONG : strcat(errorstring,"  "); break;
   /*  
  case  : strcat(errorstring,"  "); break;
    */
  default: strcat(errorstring," unknown error code "); break;
    }
  strcat(errorstring,arg);
  strcat(errorstring,"\n");
  return errorstring;
}

void Terror::error(int ernr, char *location)
{
#ifndef __GTK__
  std::cerr << get_errorstring(ernr,location);
  exit (ernr);
#else
#endif
}

void Terror::error(int ernr, char *location,char *arg)
{
#ifndef __GTK__
  std::cerr << get_errorstring(ernr,location,arg);
  exit (ernr);
#else
#endif
}

void Terror::warning(int ernr, char *location)
{
#ifndef __GTK__
  std::cerr << get_errorstring(ernr,location);
#else
#endif

}

void Terror::warning(int ernr, char *location,char *arg)
{
#ifndef __GTK__
  std::cerr << get_errorstring(ernr,location,arg);
#else
#endif 
}










