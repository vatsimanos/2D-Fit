#ifndef __ERR_H__
#define __ERR_H__

#ifndef MAX_ERROR_STRING_LENGTH
#define MAX_ERROR_STRING_LENGTH 1024
#endif

#ifndef MAX_STRING_LENGTH
#define MAX_STRING_LENGTH 256
#endif

#define OK 0
//define errors
#define OUT_OF_MEMORY 10
#define FILE_OPEN_ERROR 11
#define NO_VALUES 12
#define EMPTY_VARIABLE 13
#define ALREADY_PRESENT 14
#define NO_FILENAME 15
#define FILE_READ_ERROR 16
#define ZBRENT 17
#define SYNTAX_ERROR 18
#define OUT_OF_RANGE 19
#define READ_ERROR 20
//define warnings
#define DATA_FILE_LENGTH_WARNING 100
#define ENVVAR_NOT_SET 101
#define NOT_IMPLEMENTED 102
#define NOT_CALLED_WITH_NULL 103
#define NOT_PRESENT 104
#define WRONG_CONTENT 105
#define VALUE_NOT_MASKED 106
#define STUD_WRONG 107

#include <stdio.h>
#include <iostream>	//.h weg
#include <string.h>
#include <sstream>	//.h weg & sstream instead of strstream

class Terror
{
 private:
  std::ostringstream os;
  
  char errorstring[MAX_ERROR_STRING_LENGTH];
  char *get_errorstring(int ernr, char *location);
  char *get_errorstring(int ernr, char *location,char *arg);
  
 public:
  Terror();

  ~Terror();

  void error(int ernr, char *location);

  void error(int ernr, char *location,char *arg);

  void warning(int ernr, char *location);

  void warning(int ernr, char *location,char *arg);

  //ostream& errorstream();

  //void errorstream_show();

};


#endif /* __ERROR_H__ */
























