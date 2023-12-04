/*library to handle file file-io for patch clamp data                        */

#include "patchio.h"

void Tchannel_type::init_steady_state()
{
	std::cout<<"################ warning: not supported ##############"<<std::endl; //Tobias 2019

}

void Tchannel_type::init_symbols()
{
std::cout<<"################ warning: not supported ##############"<<std::endl; //Tobias 2019

}

double **Tchannel_type::kronecker_add_rate_constants(double **channel1, int N1,
                                 double **channel2, int N2)
{
   int i1, j1, i2, j2, N = N1 * N2;
#ifdef debug_mem_patchio
   nm++;cerr << "Tchannel_type::kronecker_add_rate_constants Result+ " << nm << endl;
#endif

   double **Result = real_matrix(1, N, 1, N);

   for (i1=1; i1<=N; i1++) {
      for (j1=1; j1<=N; j1++) {
         Result[i1][j1] = 0;
      } /* endfor */
   } /* endfor */
   for (i2=0; i2<N2; i2++) {
      for (i1=1; i1<=N1; i1++) {
         for (j1=1; j1<=N1; j1++) {
            Result[i2*N1+i1][i2*N1+j1] += channel1[i1][j1];
         } /* endfor */
      } /* endfor */
   } /* endfor */
   for (i2=0; i2<N2; i2++) {
      for (j2=0; j2<N2; j2++) {
         for (i1=1; i1<=N1; i1++) {
            Result[i2*N1+i1][j2*N1+i1] += channel2[i2+1][j2+1];
         } /* endfor */
      } /* endfor */
   } /* endfor */
   return Result;
}

Tchannel_type *Tchannel_type::reorder_states(Tchannel_type *channel)
{
 std::cout<<"################ warning: not supported ##############"<<std::endl;  //Tobias 2019
return 0;
}

Tchannel_type *Tchannel_type::reduce_states(Tchannel_type *channel)
{
   int i, j, k, no_new_states,
       *index = integer_vector(1, channel->no_states),
       *equal_index = integer_vector(1, channel->no_states);
   double akt_level, summ_i = 0, summ_j = 0,
          **k_new, **k_start_new, *level_new;

   for (i=1; i<=channel->no_states; i++) {
      equal_index[i] = 0;
   } /* endfor */
   for (i=1; i<=channel->no_states; i++) {
      j = i;
      while ((channel->level[++j] == channel->level[i]) && !equal_index[j]) {
         k = 1;
         while (k <= channel->no_states) {
            akt_level = channel->level[k];
            summ_i = 0;
            summ_j = 0;
            while (channel->level[k] == akt_level) {
               summ_i += channel->k[i][k];
               summ_j += channel->k[j][k];
               k++;
            } /* endwhile */
            if (summ_i == summ_j) {
               equal_index[j] = i;
            } else {
               equal_index[j] = 0;
               break;
            } /* endif */
         } /* endwhile */
      } /* endwhile */
   } /* endfor */
   for (i=1; i<=channel->no_states; i++) {
      if (equal_index[i]) {
         for (j=1; j<=channel->no_states; j++) {
            channel->k[j][equal_index[i]] += channel->k[j][i];
            channel->k_start[j][equal_index[i]] += channel->k_start[j][i];
         } /* endfor */
      } /* endif */
   } /* endfor */
   j = 0;
   for (i=1; i<=channel->no_states; i++) {
      if (!equal_index[i]) {
         index[++j] = i;
      } /* endif */
   } /* endfor */
   no_new_states = j;
   std::cout << "number of new states : " << no_new_states << std::endl;
#ifdef debug_mem_patchio
   nm++;cerr << "Tchannel_type::reduce_states k_new+ " << nm << endl;
#endif
   k_new = real_matrix(1, no_new_states, 1, no_new_states);
#ifdef debug_mem_patchio
   nm++;cerr << "Tchannel_type::reduce_states k_start_new+ " << nm << endl;
#endif
   k_start_new = real_matrix(1, no_new_states, 1, no_new_states);
   
   level_new = double_vector(1, no_new_states);
   for (i=1; i<=no_new_states; i++) {
      for (j=1; j<=no_new_states; j++) {
         k_new[i][j] = channel->k[index[i]][index[j]];
         k_start_new[i][j] = channel->k_start[index[i]][index[j]];
      } /* endfor */
      level_new[i] = channel->level[index[i]];
   } /* endfor */
#ifdef debug_mem_patchio
   nm--;cerr << "Tchannel_type::reduce_states channel->k- " << nm << endl;
#endif
   free_real_matrix(channel->k, 1, channel->no_states, 1, channel->no_states);
#ifdef debug_mem_patchio
   nm--;cerr << "Tchannel_type::reduce_states channel->k_start- " << nm << endl;
#endif
   free_real_matrix(channel->k_start, 1, channel->no_states, 1, channel->no_states);
   free_double_vector(channel->level, 1, channel->no_states);
   channel->k = k_new;
   channel->k_start = k_start_new;
   channel->level = level_new;
   channel->no_states = no_new_states;
   return channel;
}

Tchannel_type *Tchannel_type::merge_channels(Tchannel_type *channel1, Tchannel_type *channel2)
{
   int i, j;
   Tchannel_type *Result = (Tchannel_type *) malloc(sizeof(Tchannel_type));

   Result->no_states = channel1->no_states * channel2->no_states;
   Result->level = double_vector(1, Result->no_states);
   Result->k = kronecker_add_rate_constants(channel1->k, channel1->no_states,
                                            channel2->k, channel2->no_states);
   Result->k_start = kronecker_add_rate_constants(
                                       channel1->k_start, channel1->no_states,
                                       channel2->k_start, channel2->no_states);
   for (i=1; i<=channel1->no_states; i++) {
      for (j=1; j<=channel2->no_states; j++) {
         Result->level[(j-1)*channel1->no_states+i] =
                  channel1->level[i] + channel2->level[j];
      } /* endfor */
   } /* endfor */
   Result = reorder_states(Result);
//   Result = reduce_states(Result);
   return Result;
}

Tchannel_type::Tchannel_type(int No_states, int No_fit, char *Name)
{
  nm=0;
  nm=0;
   no_states = No_states;
   no_fit = No_fit;
   no_symbols = 0;
   if (Name!="") {
      name = Name;
   } else {
      name = NULL;
   } /* endif */
#ifdef debug_mem_patchio
   nm++;cerr << "Tchannel_type::Tchannel_type k+ " << nm << endl;
#endif
   k = real_matrix(1, no_states, 1, no_states);
#ifdef debug_mem_patchio
   nm++;cerr << "Tchannel_type::Tchannel_type k_T+ " << nm << endl;
#endif
   k_T = real_matrix(1, no_states, 1, no_states);
#ifdef debug_mem_patchio
   nm++;cerr << "Tchannel_type::Tchannel_type k_start+ " << nm << endl;
#endif
   k_start = real_matrix(1, no_states, 1, no_states);
   initial_state_value = double_vector(1, no_states);
   q = double_vector(1, no_states);
   level = double_vector(1, no_states);
   symbol = double_vector(1, no_states);
   symbol_q = double_vector(1, no_states);
   symbol_index = integer_vector(1, no_states);
   if (no_fit)
      fit_index = integer_matrix(1, no_fit, 1, 2);
   error_index = NULL;
}

Tchannel_type::~Tchannel_type()
{
  if (k!=NULL){
#ifdef debug_mem_patchio
    nm--;cerr << "Tchannel_type::~Tchannel_type k- " << nm << endl;
#endif
    free_real_matrix(k, 1, no_states, 1, no_states);}
  if (k_T!=NULL){
#ifdef debug_mem_patchio
    nm--;cerr << "Tchannel_type::~Tchannel_type k_T- " << nm << endl;
#endif
    free_real_matrix(k_T, 1, no_states, 1, no_states);}
  if (k_start!=NULL){
#ifdef debug_mem_patchio
    nm--;cerr << "Tchannel_type::~Tchannel_type k_start- " << nm << endl;
#endif
    free_real_matrix(k_start, 1, no_states, 1, no_states);}
  if (q!=NULL)
    free_double_vector(q, 1, no_states);
  if (level!=NULL)
    free_double_vector(level, 1, no_states);
  if (symbol!=NULL)
    free_double_vector(symbol, 1, no_states);
  if (symbol_q!=NULL)
    free_double_vector(symbol_q, 1, no_states);
  if (symbol_index!=NULL)
    free_integer_vector(symbol_index, 1, no_states);
  if (no_fit)
    if (fit_index!=NULL)
      free_integer_matrix(fit_index, 1, no_fit, 1, 2);
  if (error_index)
    if (error_index!=NULL)
      free_integer_matrix(error_index, 1, no_symbols, 1, no_states);
  //falsch gibt segfault??
  // HIER delete[](name);
  
}

void Tchannel_type::init()
{
   int i, j, ei = 0;
   //init_steady_state(); Tobias 2019
   //init_symbols();			Tobias 2019
   error_index = integer_matrix(1, no_symbols, 1, no_states);
   for(i = 1; i <= no_symbols; i++)
      for(j = 1; j <= no_states; j++)
         if (symbol_index[j] != i)
            error_index[i][j] = ei++;
         else
            error_index[i][j] = 0;
}

Tchannel_type * Tchannel_type::operator + (Tchannel_type &channel2)
{
   int i,j;

   Tchannel_type *Result = new Tchannel_type(no_states * channel2.no_states);

   Result->no_states = no_states * channel2.no_states;
   Result->level = double_vector(1, Result->no_states);
   Result->k = kronecker_add_rate_constants(k, no_states,
                                            channel2.k, channel2.no_states);
   Result->k_start = kronecker_add_rate_constants( k_start, no_states,
                                       channel2.k_start, channel2.no_states);
   for (i=1; i<=no_states; i++) {
      for (j=1; j<=channel2.no_states; j++) {
         Result->level[(j-1)*no_states+i] =
                  level[i] + channel2.level[j];
      } /* endfor */
   } /* endfor */
   Result = reorder_states(Result);
//   Result = reduce_states(Result);
   return Result;
}

Tchannel_type * Tchannel_type::operator * (int factor)
{
   std::cerr << "Function '*' of class Tchannel_type is not yet implemented"
	<< factor << std::endl;
   exit(1);
   return NULL;
}

Tchannel_type *copy(Tchannel_type *d,const Tchannel_type& r)
{
  if (d!=NULL)
    delete[](d);
  d=new Tchannel_type(r.no_states,r.no_fit);
  d->no_states=r.no_states;
  d->no_symbols=r.no_symbols;
  d->no_fit=r.no_fit;
  //cout << "r.name " << r.name << endl;
  if (strlen(r.name)!=0){
    d->name=new char[strlen(r.name)+1];
    d->name=strcpy(d->name,r.name);
  }
  else
    d->name=NULL;

  for (int x=1;x<=r.no_states;x++){
    for (int y=1;y<=r.no_states;y++){
    d->k[x][y]=r.k[x][y];
    //cerr << "\t" << d->k[x][y];
    d->k_T[x][y]=r.k_T[x][y];
    d->k_start[x][y]=r.k_start[x][y];
    }
    d->initial_state_value[x]=r.initial_state_value[x];
    d->q[x]=r.q[x];
    d->level[x]=r.level[x];
    d->symbol[x]=r.symbol[x];
    d->symbol_q[x]=r.symbol_q[x];
    d->symbol_index[x]=r.symbol_index[x];
  } 
  //cerr << endl;
  if (r.no_fit){
    d->fit_index = integer_matrix(1, r.no_fit, 1, 2);
    for (int x=1;x<=r.no_fit;x++){
      for (int y=1;y<=2;y++){
	d->fit_index[x][y]=r.fit_index[x][y];}
    }
  }
  d->error_index = NULL;
  return d;
}

std::ostream& operator << (std::ostream& str, const Tchannel_type &channel)
{
   int i, j;

   str << "\n#rate (simulation) matrix K:\n";
   str << channel.no_states << " " << channel.no_states << std::endl;
   for (i=1; i<=channel.no_states; i++) {
      for (j=1; j<=channel.no_states; j++) {
	// auf der Diagonalen werden Hilfswerte gespeichert
	if (i==j)
	  str << "\t0";
	else
	  str << "\t" << channel.k[i][j];
      } /* endfor */
      str << std::endl;
   } /* endfor */
   str << std::endl;
   str << "#start values of rate matrix for fitting\n";
   str << channel.no_states << " " << channel.no_states << std::endl;
   for (i=1; i<=channel.no_states; i++) {
      for (j=1; j<=channel.no_states; j++) {
	str << "\t" << channel.k_start[i][j];
      } /* endfor */
      str << std::endl;
   } /* endfor */
   str << std::endl;
   str << "#indices of the parameters in K to be fitted:\n";
   str << channel.no_fit << " 2" << std::endl;
   for (i=1; i<=channel.no_fit; i++) {
      for (j=1; j<=2; j++) {
	str << "\t" << channel.fit_index[i][j];
      } /* endfor */
      str << std::endl;
   } /* endfor */
   str << std::endl;
   str << "#initial state distribution (t=0):\n";
   str << "1 " << channel.no_states << std::endl;
   for (i=1; i<=channel.no_states;i++)
     str << "\t" << channel.initial_state_value[i];
   str << std::endl << std::endl;
   
   str << "# output level in DAC units, depending on symbol:\n";
   str << "1 " << channel.no_states << std::endl;
   for (i=1; i<=channel.no_states; i++) {
     str << "\t" << channel.level[i];
   } /* endfor */
   str << std::endl << std::endl;
   return str;
}

FILE* find_file(char *name, const char *envvar, int write, int open_file)
{
  char *Pfad, Gefunden[256] = "", *mode, wb[] = "wb", rb[] = "rb";
  FILE *Datei = 0;
  if (write)   mode = wb;
    else       mode = rb;
  if ((Datei = fopen(name,mode)) != 0) {
    if (!open_file) fclose(Datei);
    return Datei;
  } /* endif */
  else
  {
    Pfad = getenv(envvar);
    if (Pfad)
      strcpy(Gefunden, Pfad);
      else std::cout << "Die Environmentvariable " << envvar << " ist nicht gesetzt" << std::endl;  //TOBIAS 2013 auskommentiert
    if (Gefunden[strlen(Gefunden)-1] == '\\' ||
        Gefunden[strlen(Gefunden)-1] == '/')  {
    } else {
      strcat(Gefunden,"/");
    } /* endif */
    strcat(Gefunden,name);
  }
  if ((Datei = fopen(Gefunden,mode)) != 0)
  {
    strcpy(name, Gefunden);
    if (!open_file) fclose(Datei);
    return Datei;
  }
  return 0;
}

unsigned short int * load_data(const char *fname, unsigned long int &length)
{
  const int len = 256;
  unsigned long int i;
  unsigned short int *time_series;
  char Datenfile[len];
  FILE *Datafile;

  strcpy(Datenfile, fname);
  if ((Datafile = find_file(Datenfile,"DATEN",READ,OPEN)) == 0)
  { std::cerr << "error finding data file" << std::endl; exit(1); }
  else {
    if ((fread(&length,4,1,Datafile) != 1) || ferror(Datafile))
    {
      std::cerr << "error reading length of time series" << std::endl;
      exit(1);
    }
    else {
      length = convert_long_int(length);
      //cout << "length of time series: " << length << endl;
      time_series = (unsigned short int *)
                      malloc(length*sizeof(unsigned short int));
      if ((i=fread(time_series,2,length,Datafile)) != length || ferror(Datafile))
      {
        std::cerr << "error reading time series at position: " << i << std::endl;
        if (ferror(Datafile)) perror("Daten_laden");
        exit(1);
      }
      else
      std::cout << length << " samples read" << std::endl;
      for (i=0; i<length; i++) 
	time_series[i] = convert_short_int(time_series[i]);
    }
    if (fclose(Datafile))
    {
      std::cout << "error closing data file (error no: " << errno << ")" << std::endl;
    }
  }
  return (time_series-1);
}

int save_data(const char* n_fname, short int *time_series, 
	      unsigned long int length)
{
  int i;
  char fname[256];
  FILE *outfile;

  strcpy(fname, n_fname);
  outfile = find_file(fname, "DATEN", WRITE, OPEN);
  if (outfile == NULL) {
    std::cerr << "error opening file: " << fname << " for writing\n";
    exit(1);
  }
  if (fwrite(inv_convert_li(length), 4, 1, outfile) != 1) {
    std::cerr << "error writing length in file: " << fname << std::endl;
    exit(1);
  }
  for (i=1; (unsigned)i<=length; i++)
    if (fwrite(inv_convert_si(time_series[i]), 2, 1, outfile) != 1) {
      std::cerr << "error writing data in file: " << fname << " at position " << i 
	   << std::endl;
      exit(1);
    }
  return (0);
}








