#include "timeseries.h"
#include "daten.h"
 
TTimeseries::TTimeseries()
{
  Value=NULL;
  FileName=NULL;
  Num_Samples=0;
  delta_t=0,0;
  f_3dB=0;
}

TTimeseries::~TTimeseries()
{
#ifdef debug_mem_timeseries
  std::cerr << "TTimeseries::~TTimeseries delete[] FileName  " << std::endl;
#endif
  delete[] FileName;
  FileName=NULL;
#ifdef debug_mem_timeseries
  std::cerr << "TTimeseries::~TTimeseries delete[] Value  " << std::endl;
#endif
  delete[] Value;
  Value=NULL;
}

char *TTimeseries::get_filename(char *Filename)
{
  if (Filename!=NULL){
    err.warning(NOT_CALLED_WITH_NULL,(char *)"TTimeseries::get_filename",Filename);
#ifdef debug_mem_timeseries
    std::cerr << "TTimeseries::get_filename delete[] (Filename)  " << std::endl;
#endif
    delete[] (Filename);
    Filename=NULL;
  }
  else if (FileName==NULL)
    err.error(EMPTY_VARIABLE,(char *)"TTimeseries::get_filename",(char *)"FileName");
#ifdef debug_mem_timeseries
  std::cerr << "TTimeseries::get_filename new(Filename)  " << strlen(FileName)+1 
       <<std::endl;
#endif
  Filename=new char[strlen(FileName)+1];
  if (Filename!=NULL){
      strcpy(Filename,FileName);
      return Filename;}
  else{
    err.error(OUT_OF_MEMORY,(char *)"TTimeseries::FName",Filename);
    return NULL;}
};

unsigned long int TTimeseries::get_num_samples()
{
  if (Num_Samples!=0)
    return (Num_Samples);
  else{
    err.error(NO_VALUES,(char *)"TTimeseries::get_num_samples()",(char *)"Num_Samples");
    return 0;}
}

void TTimeseries::set_num_samples(unsigned long int num_samples)
{
  if (num_samples>=0)
    Num_Samples=num_samples;
  else
    err.warning(WRONG_CONTENT,(char *)"TTimeseries::set_num_samples",(char *)"Num_Samples");
}

void TTimeseries::set_value(unsigned short int *value,unsigned long int length)
{
  if (value==NULL)
    err.error(WRONG_CONTENT,(char *)"TTimeseries::set_value",(char *)"value");
  else if (Value!=NULL){
#ifdef debug_mem_timeseries
    err.warning(ALREADY_PRESENT,"TTimeseries::set_value","Value");
    std::cerr << "TTimeseries::set_value delete[](Value)\n";
#endif
    delete[](Value);
    Value=NULL;
  }
  else if (length==0)
    err.error(WRONG_CONTENT,(char *)"TTimeseries::set_value",(char *)"length");
#ifdef debug_mem_timeseries
  std::cerr << "TTimeseries::set_value new(Value)  " << length+1 << std::endl;
#endif
  Value=new unsigned short int[length+1];
  memcpy(Value+1,value+1,length*sizeof(unsigned short int));
}

void TTimeseries::set_delta_t(double deltat)
{
  if (deltat>0)
    delta_t=deltat;
  else
    err.warning(WRONG_CONTENT,(char *)"TTimeseries::set_delta_t",(char *)"delta_t");
}

void TTimeseries::set_filename(char *name)
{
  if (name!=NULL){
    if (FileName!=NULL){
#ifdef debug_mem_timeseries
    std::cerr << "TTimeseries::set_filename delete[](FileName)  " << std::endl;
#endif
      delete[](FileName);
      FileName=NULL;
    }
#ifdef debug_mem_timeseries
    std::cerr << "TTimeseries::set_filename new(FileName)  " 
	 << strlen(name)+1 << std::endl;
#endif
    FileName=new char[strlen(name)+1];
    strcpy(FileName,name);
  }
  else
    err.warning(WRONG_CONTENT,(char *)"TTimeseries::set_filename",(char *)"FileName");
}

void TTimeseries::set_f_3dB(double f3dB)
{
  if ((f3dB>0)&&(f3dB<1/delta_t))
    f_3dB=f3dB;
  else
    err.warning(WRONG_CONTENT,(char *)"TTimeseries::set_filename",(char *)"f_3dB , try to set delta_t first");
}

int TTimeseries::load(char *FName)
{
  //if (FileName==NULL){ // nur laden falls noch nicht im Speicher
    if (0==0){
    if (FName!=NULL){
      if (FileName!=NULL){
#ifdef debug_mem_timeseries
	std::cerr << "TTimeseries::load delete[](FileName)  " << std::endl;
#endif
	delete[](FileName);
	FileName=NULL;
}
#ifdef debug_mem_timeseries
      std::cerr << "TTimeseries::load new(FileName)  " 
	   << strlen(FName)+1 << std::endl;
#endif
      FileName=new char[strlen(FName)+1];
      if (FileName==NULL){
	err.error(OUT_OF_MEMORY,(char *)"TTimeseries::load",FName);
	return -1;}
      else
	strcpy(FileName,FName);}
#ifdef debug_timeseries
    std::cerr << FName << std::endl;
#endif
    std::ifstream in(FName,std::ios::in|std::ios::binary);
    if  (!in)
      {
	char buffer[MAX_STRING_LENGTH], daten[]="DATEN";
	const char *Pfad;
	Pfad = getenv(daten);
	if (Pfad)
	  strcpy(buffer, Pfad);
	else err.warning(ENVVAR_NOT_SET,(char *)"TTimeseries::load",daten);
      if  (! (buffer[strlen(buffer)-1] == '\\' ||
	      buffer[strlen(buffer)-1] == '/'))  
	strcat(buffer,"/");
      strcat(buffer,FName);
      in.open(buffer,std::ios::in);
      if (in)
	{
#ifdef debug_mem_timeseries
	std::cerr << "TTimeseries::load delete[](FileName)  " << std::endl;
#endif
	  delete[](FileName);
	  FileName=NULL;
#ifdef debug_mem_timeseries
	  std::cerr << "TTimeseries::load new(FileName)  " 
	       << strlen(FName)+1 << std::endl;
#endif
	  FileName=new char[strlen(buffer)+1];
	  strcpy(FileName,buffer);
	}
      std::cerr << buffer << std::endl;
      }

    if (!in)
      {
	err.error(FILE_OPEN_ERROR,(char *)"TTimeseries::load",FName);   
	return -1;
      }
    in.read((char *)&(Num_Samples), 4);  
    in.seekg(0, std::ios::end);                                // Filelaenge in bytes
    if (Num_Samples != (unsigned long int)(((int)in.tellg() - 4) / 2))
      {
	Num_Samples = (((int)in.tellg() - 4) / 2);
	err.warning(DATA_FILE_LENGTH_WARNING,(char *)"TTimeseries::load",FileName);
	in.seekg(0, std::ios::beg);
      }
    else
      in.seekg(4, std::ios::beg);
#ifdef debug_timeseries
    // std::cerr << Num_Samples << "sizeof short int" << sizeof(short int) << std::endl;    
#endif
    if (Value!=NULL){
#ifdef debug_mem_timeseries
    std::cerr << "TTimeseries::load delete[](Value)\n";
#endif
      delete[](Value);
      Value=NULL;
    }
#ifdef debug_mem_timeseries
    std::cerr << "TTimeseries::load new(Value)  " 
	 << Num_Samples+1 << std::endl;
#endif
    Value = new unsigned short int[Num_Samples+1];  // pointer auf neue daten
    if (Value==NULL)
      err.error(OUT_OF_MEMORY,(char *)"TTimeseries::load",(char *)"Value");
    // Value[1] is first Value of timeseries 
    // not very computer science like but thilo wants it this way
    in.read((char *)Value+1, Num_Samples * 2);   // werte sind 2byte gross 
#ifdef debug_timeseries  
    std::cerr << FileName << " : " << std::endl;
    //for (int ui=0;ui<=20;ui++)
    //  std::cerr << "befor Value [" << ui << "] = " << Value[ui] << std::endl;
#endif
    for (unsigned int i=0; i<=Num_Samples; i++) 
	Value[i] = convert_short_int(Value[i]);
#ifdef debug_timeseries  
    // std::cerr << FileName << " : " << std::endl;
    // for (int ui=0;ui<=20;ui++)
    //  std::cerr << "Value [" << ui << "] = " << Value[ui] << std::endl;
#endif
    in.close();
    return 0;
  }
  else
    err.error(ALREADY_PRESENT,(char *)"TLeveldet::load",FileName);
}

void TTimeseries::save(char *FName)
{
  int i;	
  std::ofstream out(FName,std::ios::out|std::ios::binary|std::ios::trunc);
  if (!out)
    err.error(FILE_OPEN_ERROR,(char *)"TTimeseries::save",FName);
  //Tobias:anders; funktioniert nur f�r die Zeitreihe mit Daten.Wert
  
  out.write((char *)&Daten.Anz_Samples,sizeof(Daten.Anz_Samples));
  std::cout<<"ANZ_SAMPLES:"<<Daten.Anz_Samples<<std::endl;
  for (i=1; i <= Daten.Anz_Samples; i++)
    out.write((char *)&Daten.Wert[i],sizeof(Daten.Wert[i]));
  
  out.close();
}

//Tobias: zum Ausgeben der Amplitudenhistogramme f�r Phillips direkten Zeitreihenfit!
void TTimeseries::save_amplhis(char *FName)
{
	
}	


std::ostream& operator<< (std::ostream& os, const TTimeseries& ts)
{
  os << "name :" << ts.FileName << std::endl;
  os << ts.Num_Samples << " samples with " << (1/ts.delta_t) << "Hz" << std::endl;
  for (unsigned int i=1;i<=20;i++)
    os << i << "\t" << ts.Value[i] << std::endl;
  for (unsigned int i=1;i<=20;i++)
    os << ts.Num_Samples-i << "\t" << ts.Value[ts.Num_Samples-i] << std::endl; 
  return os;
}

void TTimeseries::clear()
{
  //std::cerr << "TTimeseries::clear()\n";
  if (Value!=NULL){
    delete[](Value);
    Value=NULL;
  }
  else
    err.warning(NO_VALUES,(char *)"TTimeseries::clear_mem",(char *)"Value");
  if (FileName!=NULL){
    delete[](FileName);
    FileName=NULL;
  }
  else
    err.warning(NO_VALUES,(char *)"TTimeseries::clear_mem",(char *)"FileName");
  //std::cerr << "TTimeseries::clear() end\n";
}




