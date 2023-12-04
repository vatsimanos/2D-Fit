#include "setfile.h"
 
TSetfile::TSetfile()
{
  line_counter=0;
}

void TSetfile::init(double deltat,double sigm, double zerocurr, double f3dB, int nochanneltypes, int *nochanneltype, int nosamples, char *dataf, char *rauschfile, char *stepresp, Tchannel_type **channeltype)
{
  if (no_channels_type!=NULL){
    delete[](no_channels_type);
    no_channels_type=NULL;
#ifdef debug_mem_setfile
    cerr << "TSetfile::init delete[](no_channels_type) " << endl;
#endif
  }
  if (channel_type!=NULL){
    for (int i=1;i<=no_channel_types;i++)
      if (channel_type[i]!=NULL){
#ifdef debug_mem_setfile
	cerr << "TSetfile::init delete(channel_type[" << i << "])\n";
#endif
	delete(channel_type[i]);
	channel_type[i]=NULL;
      }
#ifdef debug_mem_setfile
    cerr << "TSetfile::init delete[](channel_type)\n";
#endif
    delete[](channel_type);
    channel_type=NULL;
  }
  delta_t=deltat;
  sigma=sigm;
  zero_current=zerocurr;
  f_3dB_digfi=f3dB;
  no_channel_types=nochanneltypes;
#ifdef debug_mem_setfile
    cerr << "TSetfile::init new(no_channels_type)\n";
#endif
  no_channels_type=new int[no_channel_types+1];
  no_channels_total=0;
  for (int i=1;i<=no_channel_types;i++){
    no_channels_type[i]=nochanneltype[i];
    no_channels_total+=no_channels_type[i];}
  line_counter=0;
  no_samples=nosamples;
  strcpy(datafile,dataf);
  strcpy(rausch_file,rauschfile);
  strcpy(stepresponse_file,stepresp);
#ifdef debug_mem_setfile
    cerr << "TSetfile::init new(channel_type)\n";
#endif
  channel_type=new Tchannel_type*[no_channel_types+1];
  for (int i=1;i<=no_channel_types;i++){
    channel_type[i]=NULL;
    channel_type[i]=copy(channel_type[i],*channeltype[i]);
#ifdef debug_mem_setfile
    //cerr << *channel_type[i] << endl;
#endif
  }
}

TSetfile::~TSetfile()
{
  
  // do not delete arrays !!!!!!!!!!!
  //delete(datafile);
  //delete(rausch_file);
  //delete(stepresponse_file);


#ifdef debug_mem_setfile
    cerr << "TSetfile::~TSetfile delete[](no_channels_type)\n";
#endif
  delete[](no_channels_type);
  no_channels_type=NULL;
  /*
  for (int i=1;i<=no_channel_types ;i++)
    if (channel_type[i]!=NULL){
#ifdef debug_mem_setfile
      cerr << "TSetfile::~TSetfile delete(channel_type[" << i << "])\n";
#endif      
      delete(channel_type[i]);
      channel_type[i]=NULL;
}
#ifdef debug_mem_setfile
  cerr << "TSetfile::~TSetfile delete[](channel_type)\n";
#endif  
  delete[](channel_type);
  channel_type=NULL;
  */
}

long int TSetfile::get_no_samples()
{ return no_samples; }

double TSetfile::get_sigma()
{ return sigma; }

double TSetfile::get_zero_current()   
{ return zero_current;  }

double TSetfile::get_single_channel_current()  //by Tobias
 {
  unsigned short int i;
  double current;
  
  current = 0;
  for (i=1; i <= channel_type[1]->no_states; i++)
   if (fabs(channel_type[1]->level[i]) > fabs(current))
     current = channel_type[1]->level[i];
  return current;  
 }  
 
int TSetfile::get_no_rconst()   // by Tobias
 {
  int i,j;
  int nr=0;
   	
  for (i=1; i<=channel_type[1]->no_states; i++)
    for (j=1; j<=channel_type[1]->no_states; j++)
      {
        if  (channel_type[1]->k[i][j] > 0 )
            nr = nr+1;   
	//cout <<i<<":"<<j<<":"<<channel_type[1]->k[i][j]<<endl; 	                                    
      }	
  return nr;	
 }  	
  
double TSetfile::get_rconst(unsigned int nr_rconst)
 {
  int i,j;
  unsigned int nr=0;
  double rconst=-1;
    	
  for (i=1; i<=channel_type[1]->no_states; i++)
    for (j=1; j<=channel_type[1]->no_states; j++)
        if  (channel_type[1]->k[i][j] > 0 )
          {
            nr = nr+1;   
            if (nr == nr_rconst)
               rconst = channel_type[1]->k[i][j];
 	  }
 	  	         
  return rconst;    
 }       

 void TSetfile::set_rconst(unsigned int nr_rconst, double rconst)
  {
   int i,j;
   unsigned int nr=0;
   double diag;
     	
   for (i=1; i<=channel_type[1]->no_states; i++)
     for (j=1; j<=channel_type[1]->no_states; j++)
        if  (channel_type[1]->k[i][j] > 0 )
          {       
             nr = nr+1;   
             if (nr == nr_rconst)
                channel_type[1]->k[i][j] = rconst;
          }	
   //Setzen der Diagonale       
   for (i=1; i<=channel_type[1]->no_states; i++)
     {
      diag=0;	
      for (j=1; j<=channel_type[1]->no_states; j++)              
         if (i != j)
            diag = diag + channel_type[1]->k[i][j];
      channel_type[1]->k[i][i] = -diag;      
     } 
  }

void TSetfile::set_sigma (double noise)				//Tobias 2021
 	{
 	 sigma = noise;
 	}  

double TSetfile::get_rconstxy (int x, int y)
 { return channel_type[1]->k[x][y];} 

int TSetfile::get_matrix_dimension()
 { return channel_type[1]->no_states; }

double TSetfile::get_delta_t()
{ return delta_t; }

double TSetfile::get_f3dB()
{ return f_3dB_digfi; }

int TSetfile::get_no_channels_total()
{ return no_channels_total; }
 
int TSetfile::get_no_channel_types()
{ return no_channel_types; }

char *TSetfile::get_timeseries_filename(char *Filename)
{
  if (datafile!=NULL){
    if (Filename!=NULL){
      err.warning(NOT_CALLED_WITH_NULL,(char *)"TSetfile::get_timeseries_filename",Filename);
      delete[](Filename);
      Filename=NULL;
    }
    Filename=new char[strlen(datafile)+1];
    strcpy(Filename,datafile);
    return Filename;
  }
  else{
    err.warning(NOT_PRESENT,(char *)"TSetfile::get_timeseries_filename",(char *)"datafile");
    return NULL;
  }
}

char *TSetfile::get_noiseseries_filename(char *Filename)
{
  if (datafile!=NULL){
    if (Filename!=NULL){
      err.warning(NOT_CALLED_WITH_NULL,(char *)"TSetfile::get_timeseries_filename",
		  Filename);
      delete[](Filename);
      Filename=NULL;
    }
    Filename=new char[strlen(datafile)+1+6];
    if (strrchr(strcpy(Filename,datafile),'.')!=NULL)
      *strrchr(Filename,'.')='\0';
    strcat(Filename,".noise");
    //cerr << Filename << endl;
    return Filename;
  }
  else{
    err.warning(NOT_PRESENT,(char *)"TSetfile::get_timeseries_filename",(char *)"datafile");
    return NULL;}
}


char *TSetfile::get_rausch_filename(char *Filename)
{
  if (rausch_file!=NULL){
    if (Filename!=NULL){
      err.warning(NOT_CALLED_WITH_NULL,(char *)"TSetfile::get_rausch_filename",
		  Filename);
      delete[](Filename);
      Filename=NULL;
    }
    Filename=new char[strlen(rausch_file)+1];
    strcpy(Filename,rausch_file);
    return Filename;
  }
  else{
    err.warning(NOT_PRESENT,(char *)"TSetfile::get_rausch_filename",(char *)"rausch_file");
    return NULL;
  }
}

char *TSetfile::get_stepresponse_filename(char *Filename)
{
  if (stepresponse_file!=NULL){
    if (Filename!=NULL){
      err.warning(NOT_CALLED_WITH_NULL,(char *)"TSetfile::get_rausch_filename",
		  Filename);
      delete[](Filename);
      Filename=NULL;
    }
    Filename=new char[strlen(stepresponse_file)+1];
    strcpy(Filename,stepresponse_file);
    return Filename;
  }
  else{
    err.warning(NOT_PRESENT,(char *)"TSetfile::get_rausch_filename",
		(char *)"stepresponse_file");
    return NULL;
  }
}


char *TSetfile::next_line(FILE *Fsetfile)
{
   int i, len;
   bool booldummy;  //2019
   static char line[256], sdummy[256];

   do {
      booldummy = fgets(sdummy, 256, Fsetfile);   //2019
      line_counter++;
   } while (sdummy[0] == '#' || sdummy[0] == '\n' || sdummy[0] == '\r'); /* enddo */
        len = strcspn(sdummy,"#");
   strncpy(line, sdummy, len);
   for (i=0; i<len; i++) {
      line[i] = tolower(line[i]);
   } /* endfor */
   line[len] = '\0';
   return line;
}

double **TSetfile::read_dMatrix(FILE *Fsetfile, int *Zeilen, int *Spalten)
{
   double **M;
   int i,j, idummy;
   char *line, *line_end;
   char buffer[256]; 
   line = next_line(Fsetfile);
   idummy = sscanf(line, "%d %d", Zeilen, Spalten);		//tobias 2017 
   if (idummy!=2) {
      sprintf(buffer,"error reading matrix dimensions in line %d", 
	      line_counter);
      err.error(SYNTAX_ERROR,(char *)"TSetfile::read_dMatrix",buffer);
   }
   M = real_matrix(1, *Zeilen, 1, *Spalten);

   for (i=1; i<=*Zeilen; i++) {	 
      line = next_line(Fsetfile);
      line_end = line + strlen(line);
      for (j=1; j<=*Spalten; j++) {
         line += strspn(line, " \t");
         if (line < line_end) {
	       //std::cout << "\t" << M[i][j];
	   idummy = sscanf(line,"%lf", &(M[i][j]));
	   // if (i==j) 
	   // cout << " " << M[i][j];
         }
         else {
	   sprintf(buffer,"exeded line end while reading matrix in line %d"
		   "\nline pointer = %p, line end = %p\n",line_counter,line
		   ,line_end);
	   err.error(SYNTAX_ERROR,(char *)"TSetfile::read_dMatrix",buffer);
	 }
         line += strcspn(line, " \t");
         if (idummy != 1) {
            sprintf(buffer,"error reading setfile in line %d\n"
                    "while reading double matrix row %d column %d\n",
                    line_counter, i, j);
            err.error(SYNTAX_ERROR,(char *)"TSetfile::read_dMatrix",buffer);
         } /* endif */
      } /* endfor */
   } /* endfor */ 
   return M;
}

int **TSetfile::read_iMatrix(FILE *Fsetfile, int* Zeilen, int *Spalten)
{
   int i,j, idummy;
   int **Matrix;
   char *line;
   char buffer[256];

   line = next_line(Fsetfile);
   idummy = sscanf(line, "%d %d", Zeilen, Spalten);  //Tobias 2017 
   if (idummy!=2) {
      sprintf(buffer,"error reading matrix dimensions in line %d", 
	      line_counter);
      err.error(SYNTAX_ERROR,(char *)"TSetfile::read_iMatrix",buffer);
   }
   Matrix = integer_matrix(1, *Zeilen, 1, *Spalten);
   for (i=1; i<=*Zeilen; i++) {
      line = next_line(Fsetfile);
      for (j=1; j<=*Spalten; j++) {
         line += strspn(line, " \t");
         idummy = sscanf(line,"%d",&Matrix[i][j]);
         line += strcspn(line, " \t");
         if (idummy != 1) {
            sprintf(buffer,"error reading setfile in line %d\n"
                    "while reading integer matrix row %d column %d\n",
                    line_counter, i, j);
            err.error(SYNTAX_ERROR,(char *)"TSetfile::read_iMatrix",buffer);
         } /* endif */
      } /* endfor */
   } /* endfor */
   return Matrix;
}

double TSetfile::read_setfile(char *setfile_name)
{  	
	
   int i, j, k, idummy, Z, S, Z1, S1, no_states, no_channel_type_1, **ppidummy,
       no_fit, **fitdummy, n_lines;
   FILE *FHsetfile;
   double version, **ppddummy, **kdummy, **k_startdummy, *leveldummy=NULL, *initial_state_dummy,sum=0;
   char filename[256], *line, sdummy[256];
   char buffer[366];
   bool booldummy;  //2019

   strcpy(filename, setfile_name);
   if (strstr(filename,".set") == 0) strcat(filename,".set");
   FHsetfile = find_file(filename, "SETFILES", 0, 1);
   if (FHsetfile == 0) {	
     //err.error(FILE_OPEN_ERROR,"TSetfile::read_setfile",setfile_name);
     sprintf (Daten.errormessage,"Error opening Setfile %s",setfile_name);
     return 0;
   } else {
   } /* endif */
   rewind(FHsetfile);
   idummy = fscanf(FHsetfile, "# SV %lf", &version);
   if (idummy != 1) {
      version = 1;
   } /* endif */
   booldummy = fgets(sdummy, 256, FHsetfile); //2019
   line_counter = 1;
   if (version >= 3.0) n_lines = 9;
   else
     if (version >= 2.1) n_lines = 8;
     else n_lines = 7;
   for (i=1; i<=n_lines; i++) {
      line = next_line(FHsetfile);
      idummy = 0;
      idummy += sscanf(line,"delta%*c%*c%*c%lf", &delta_t);
      idummy += sscanf(line,"datafile%*c%s", datafile);
      idummy += sscanf(line,"no%*cchannel%*ctypes%*c%d", &no_channel_types);
      idummy += sscanf(line,"channel%*cno%*c%d", &no_channel_type_1);
      idummy += sscanf(line,"nchan%*c%d", &no_channel_type_1);
      idummy += sscanf(line,"no%*csamples%*c%lu", &no_samples);
      idummy += sscanf(line,"noise%*cfile%*c%s", rausch_file);
      idummy += sscanf(line,"rauschen%*c%s", rausch_file);
      idummy += sscanf(line,"sigma%*c%lf", &sigma);
      idummy += sscanf(line,"sigma%*cnoise%*c%lf", &sigma);
      idummy += sscanf(line,"stepresponse%*c%s", stepresponse_file);
      idummy += sscanf(line,"sprungantwort%*c%s", stepresponse_file);
      idummy += sscanf(line,"zero%*ccurrent%*c%lf", &zero_current);
      idummy += sscanf(line,"f%*c3db%*cdigfi%*c%lf", &f_3dB_digfi);
      if (idummy != 1) {
         sprintf(buffer,"error while reading header of setfile in line "
		 "%d\n%s\n",line_counter,line);	 
         //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
         strcpy(Daten.errormessage,buffer);
         return 0;
      } /* endif */
   } /* endfor i n_lines */
   if (version == 1) {
      no_channel_types = 1;
   } /* endif */
   if (version <3.0 )
     f_3dB_digfi=1/(4*delta_t);
   no_channels_total = 0;
   channel_type =
      (Tchannel_type **) malloc(sizeof(Tchannel_type*)*(no_channel_types+1));
   if (channel_type == NULL)
     err.error(OUT_OF_MEMORY,(char *)"TSetfile::read_setfile",(char *)"channel_type");
   no_channels_type = (int *) malloc(sizeof(int) *
				     no_channel_types+1); 
   for (i=1; i<=no_channel_types; i++) {
      if (version >= 2) {
         line = next_line(FHsetfile);
         line = strstr(line,"=");
         if (line == NULL) {
            sprintf(buffer,"no \'=\' character in line %d\n", line_counter);
            strcpy(Daten.errormessage,buffer);
            return 0;
            //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
         } /* endif */
         line += 1;
         // einlesen der Anzahl der Kanaele des Typs i
	 idummy = sscanf(line, "%d", &(no_channels_type[i]));
         if (idummy != 1) {
            sprintf(buffer, "error reading number of equal channels of channel"
                            " type %d in line %d\n", i, line_counter);
            strcpy(Daten.errormessage,buffer);
            return 0;                            
            //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
      } /* endif */
      } else {
         no_channels_type[i] = no_channel_type_1;
         line = next_line(FHsetfile);
      } /* endif */
      // zaehlen wieviele Kanaele insgesamzt existieren
      no_channels_total += no_channels_type[i];
      // Einlesed der simulation matrix
      kdummy = read_dMatrix(FHsetfile, &Z, &S);
      if (Z != S) {
         sprintf(buffer,"matrix of rate constants of channel type %d "
		 "is not quadratic (line %d)\n",i, line_counter);
         strcpy(Daten.errormessage,buffer);
         return 0;		 
	 //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
      } /* endif */
      no_states = Z;
      if (version >= 2) {
	// Einlesen der start value matrix
	k_startdummy = read_dMatrix(FHsetfile, &Z, &S);
	if (Z != S) {
	  sprintf(buffer,"matrix of start values of channel type %d "
		  "is not quadratic (line %d)\n",i, line_counter);
          strcpy(Daten.errormessage,buffer);
          return 0;		  
	  //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
	} /* endif */
	if (Z != no_states) {
	  sprintf(buffer,"matrix of start values of channel type %d "
		  "has a wrong dimension (line %d)\n",i, line_counter);
          strcpy(Daten.errormessage,buffer);
          return 0;		  
	  //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
	} /* endif */
      } else {
	k_startdummy = real_matrix(1, no_states, 1, no_states);
	// start value matrix selbst generieren, falls version < 2
	for (j=1; j<=no_states; j++) {
	  for (k=1; k<=no_states; k++) {
	    k_startdummy[j][k] =
	      kdummy[j][k];
	  } /* endfor */
	} /* endfor */
      } /* endif */    
      if (version <2) line = next_line(FHsetfile);
      // Einlesen der Indizes die gefittet werden sollen  
      fitdummy = read_iMatrix(FHsetfile, &Z, &S);		      //fitdummy = read_iMatrix(FHsetfile, &Z, &S);      
      if (Z > no_states*(no_states-1)) {
         sprintf(buffer,"too much elements to fit (%d) are choosen for "
		 "channel type %d above line %d\n(it don't makes sense to "
		 "fit the same parameter twice or to fit the diagonal\n"
		 "elements)\n", Z, i, line_counter);
	 strcpy(Daten.errormessage,buffer);
         return 0;	 
         //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
      } /* endif */
      if (S != 2) {
         sprintf(buffer,"a fit matrix has exactly 2 columns, not %d (above "
		 "line %d)\n",S, line_counter);
	 strcpy(Daten.errormessage,buffer);
         return 0;	 
         //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
      } /* endif */
      no_fit = Z;
      // Einlesen initial state distribution
      if (version < 2.1)
	line = next_line(FHsetfile);
      if ((version < 2.1)||(version>=3.0)) {
	//line = next_line(FHsetfile);

	ppddummy = read_dMatrix(FHsetfile, &Z, &S); 
	if (Z>1){
	  sprintf(buffer,"initial state matrix has exactly 1 line, not %d "
		  "(above line %d)\n",Z, line_counter);
          strcpy(Daten.errormessage,buffer);
          return 0;		  
	  //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
	}
	sum=0.0;
	for (j=1; j<=no_states; j++)
	  sum += ppddummy[1][j];
	double one=1.0;
	if ((sum-one)>0.0000009){
	  sprintf(buffer,"sum of the initial state distribution values has "
		  "to be exactly 1.0, not %f (above line %d)\n",sum, 
		  line_counter);
	  strcpy(Daten.errormessage,buffer);
          return 0;	  
	  //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
	}
	initial_state_dummy=(double *)malloc(sizeof(double)*no_states)-1;  
	for (j=1; j<=no_states; j++) {
	  // cout << "no_states " << no_states << "\tj " << j << endl;
	  initial_state_dummy[j] = ppddummy[1][j];
	}
	free_real_matrix(ppddummy,1,Z,1,S);			//TObias 2019 free_dmatrix(ppddummy,1,Z,1,S);
      } /* endif */
      else // genereate initial state value if not given in setfile 
	{
	  initial_state_dummy=(double *)malloc(sizeof(double)*no_states)-1;
	  for (j=1; j<=no_states; j++) {
	    initial_state_dummy[j] = 1 / no_states;
	  }
	}
      // Einlesen der levels 
      if (version < 2.1) {
         leveldummy=(double *)malloc(sizeof(double)*no_states)-1;  
         if (version <2) line = next_line(FHsetfile);
         ppidummy = read_iMatrix(FHsetfile, &Z, &S);
         if (version <2) line = next_line(FHsetfile);
         ppddummy = read_dMatrix(FHsetfile, &Z1, &S1);
         for (j=1; j<=no_states; j++) {
            leveldummy[j] = ppddummy[1][ ppidummy[1][j] ] - ppddummy[1][1];
         } /* endfor */
         if (i == 1) zero_current = ppddummy[1][1];
         free_integer_matrix(ppidummy, 1, Z, 1, S);
         free_real_matrix(ppddummy, 1, Z1, 1, S1);
      } else {
         ppddummy = read_dMatrix(FHsetfile, &Z, &S);
	 if (Z==1){
	   leveldummy=(double *)malloc(sizeof(double)*no_states)-1;
	   for (j=1; j<=no_states; j++) 
	     leveldummy[j] = ppddummy[1][j];
	   free_real_matrix(ppddummy,1,Z,1,S);
	   // for (j=1;j<=S;j++)
	   //   cout << "leveldummy " << leveldummy[j] << endl;
	 }
	 else{
	   sprintf(buffer,"output level matrix has exactly 1 row, not %d "
		   "(above line %d)\n",Z,line_counter);
           strcpy(Daten.errormessage,buffer);
           return 0;		   
	   //err.error(SYNTAX_ERROR,"TSetfile::read_setfile",buffer);
	   }
      } /* endif */      
#ifdef debug_mem_setfile
      cerr << "TSetfile::read_setfile new(channel_type[" << i << "])\n";
#endif
      channel_type[i] = new Tchannel_type(no_states, no_fit);
      channel_type[i]->k         = kdummy;
      channel_type[i]->k_start   = k_startdummy;
      channel_type[i]->initial_state_value = initial_state_dummy;
      channel_type[i]->fit_index = fitdummy;
      channel_type[i]->level     = leveldummy;
      channel_type[i]->init();
   } /* endfor */   
   return version;
}

// writes Version 3.0 setfile to name
void TSetfile::write_setfile(char *name)
{
  if (name==NULL){
    err.warning(EMPTY_VARIABLE,(char *)"Setfile::write_setfile",(char *)"name empty, setfile not written");}
  else{
#ifndef debug_setfile
    std::ofstream setfile(name,std::ios::out);
    //ofstream setfile(name,ios::out|ios::noreplace);
    if (setfile.fail())
      err.error(FILE_OPEN_ERROR,(char *)"TSetfile::write_setfile",name);
    else{
#endif
#ifdef debug_setfile
      cout
#else
      setfile 
#endif
	<< "# SV 3.0  (version of the setfile structure, has to be in "
	<< "the first line)\n" 
	<< "# lines starting with a '#' are treated as a comment and are "
	<< "ignored by the\n# program\n#"
	<< "\ndelta t\t\t\t" << delta_t
	<< "\ndatafile\t\t" << datafile
	<< "\nno channel types\t" << no_channel_types
	<< "\nno samples\t\t" << no_samples
	<< "\nnoise file\t\t" << rausch_file
	<< "\nsigma\t\t\t" << sigma
	<< "\nstepresponse\t\t" << stepresponse_file
	<< "\nzero current\t\t" << zero_current
	<< "\nf 3dB digfi\t\t" << f_3dB_digfi << std::endl << std::endl;
	for (int i=1;i<=no_channel_types;i++)
	  {
#ifdef debug_setfile
	    cout
#else
	    setfile
#endif
	      << "\nno type " << i <<" channels\t=" << no_channels_type[i]
	      << " #number of equal channels has to be behind a '='"
	      << *channel_type[i];
	  } // end for no_channel_types
#ifndef debug_setfile
    } // end if setfile.fail()
    setfile.close();
#endif
  } // end if name==NULL
} // end write_setfile



void TSetfile::clear()
{
  //cerr << "TSetfile::clear()\n";
  if (no_channels_type!=NULL){
    delete[](no_channels_type);
    no_channels_type=NULL;
#ifdef debug_mem_setfile
    cerr << "TSetfile::init delete[](no_channels_type) " << endl;
#endif
  }
  if (channel_type!=NULL){
    for (int i=1;i<=no_channel_types;i++)
      if (channel_type[i]!=NULL){
#ifdef debug_mem_setfile
	cerr << "TSetfile::init delete(channel_type[" << i << "])\n";
#endif
	delete(channel_type[i]);
	channel_type[i]=NULL;
      }
#ifdef debug_mem_setfile
    cerr << "TSetfile::init delete[](channel_type)\n";
#endif
    delete[](channel_type);
    channel_type=NULL;
}
  //cerr << "TSetfile::clear() end\n";
}

#include "simulat.c"

