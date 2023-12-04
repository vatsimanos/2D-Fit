/*******************************************************************
 *  Kiel-Patch
 *  1999 Oliver Radomski 
 *  based on uwe kirst's daytest 
 *  CAU KIEL
 *******************************************************************/

#include <math.h>

#include "declare.h"
#include "daten.h"
#include "fit.h"
#include "round.h"
#include <chrono>	//2019
#include <thread>	//2019

TSettings::TSettings()
{
  Ampl = Dagan;
  VFactor = VDagan;       

  Zoom = false;
  
  Gain = 50.0;
  Samplininterval = 5.0;
  PositiveMembranPotential = true;

  Zoomed.StackedWindows = 1;
  Zoomed.AvarageFilter = normal;
  Zoomed.LengthOfAvarageFilter = 4;

  UnZoomed.StackedWindows = 4;
  UnZoomed.AvarageFilter = normal;
  UnZoomed.LengthOfAvarageFilter = 20;

  enlarge_current = false;

  SamplesToSkip = 0;
  SamplesToProcess = 0;

//BEGIN KARSTEN (von 0,0 auf 2000,2100 ge�ndert)
  CurrentFrom = 2000;
  CurrentTo = 2100;
//END KARSTEN

  multi_samples = 54000;

  x = y = w = h = 0;
}

short   
TSettings::StackedWindows()
{
  if (Zoom)
  {
    return Zoomed.StackedWindows;
  } else {
    return UnZoomed.StackedWindows;
  }
}

Mittelwertbildung     
TSettings::AvarageFilter()
{
  if (Zoom)
  {
    return Zoomed.AvarageFilter;
  } else {
    return UnZoomed.AvarageFilter;
  }
}

short    
TSettings::LengthOfAvarageFilter()
{
  if (Zoom)
  {
    if (Zoomed.AvarageFilter == keine) return 1;
    return Zoomed.LengthOfAvarageFilter;
  } else {
    if (UnZoomed.AvarageFilter == keine) return 1;
    return UnZoomed.LengthOfAvarageFilter;
  }
}

void 
TSettings::load(std::istream& in)	//std::istream& in 
{
// zum laden aus setfile
  in.read((char *)&Ampl, sizeof(Ampl));
  in.read((char *)&VFactor, sizeof(VFactor));  
  in.read((char *)&Gain, sizeof(Gain));
  in.read((char *)&Samplininterval, sizeof(Samplininterval));
  in.read((char *)&PositiveMembranPotential, sizeof(PositiveMembranPotential));

  in.read((char *&)Zoomed.StackedWindows, sizeof(Zoomed.StackedWindows));
  in.read((char *&)Zoomed.AvarageFilter, sizeof(Zoomed.AvarageFilter));
  in.read((char *&)Zoomed.LengthOfAvarageFilter, sizeof(Zoomed.LengthOfAvarageFilter));

  in.read((char *)&UnZoomed.StackedWindows, sizeof(UnZoomed.StackedWindows));
  in.read((char *)&UnZoomed.AvarageFilter, sizeof(UnZoomed.AvarageFilter));
  in.read((char *)&UnZoomed.LengthOfAvarageFilter, sizeof(UnZoomed.LengthOfAvarageFilter));

  in.read((char *)&enlarge_current, sizeof(enlarge_current));

  in.read((char *)&Zoom, sizeof(Zoom));

  in.read((char *)&SamplesToSkip, sizeof(SamplesToSkip));
  in.read((char *)&SamplesToProcess, sizeof(SamplesToProcess));
  in.read((char *)&CurrentFrom, sizeof(CurrentFrom));
  in.read((char *)&CurrentTo, sizeof(CurrentTo));

  int l;
  in.read((char *)&l, sizeof(l));
  in.read((char *)&Daten.cdir, l * sizeof(char));

  in.read((char *)&x, sizeof(x));
  in.read((char *)&y, sizeof(y));
  in.read((char *)&w, sizeof(w)); 
  in.read((char *)&h, sizeof(h)); 
}

void 
TSettings::store(std::ostream& out)
{
/*2019	
// abspeichern ins jobfile
  out.write((char *)&Ampl, sizeof(Ampl));
  out.write((char *)&VFactor, sizeof(VFactor));  
  
  out.write((char *)&Gain, sizeof(Gain));
  out.write((char *)&Samplininterval, sizeof(Samplininterval));
  out.write((char *)&PositiveMembranPotential, sizeof(PositiveMembranPotential));

  out.write((char *)&Zoomed.StackedWindows, sizeof(Zoomed.StackedWindows));
  out.write((char *)&Zoomed.AvarageFilter, sizeof(Zoomed.AvarageFilter));
  out.write((char *)&Zoomed.LengthOfAvarageFilter, sizeof(Zoomed.LengthOfAvarageFilter));

  out.write((char *)&UnZoomed.StackedWindows, sizeof(UnZoomed.StackedWindows));
  out.write((char *)&UnZoomed.AvarageFilter, sizeof(UnZoomed.AvarageFilter));
  out.write((char *)&UnZoomed.LengthOfAvarageFilter, sizeof(UnZoomed.LengthOfAvarageFilter));

  out.write((char *)&enlarge_current, sizeof(enlarge_current));

  out.write((char *)&Zoom, sizeof(Zoom));

  out.write((char *)&SamplesToSkip, sizeof(SamplesToSkip));
  out.write((char *)&SamplesToProcess, sizeof(SamplesToProcess));
  out.write((char *)&CurrentFrom, sizeof(CurrentFrom));
  out.write((char *)&CurrentTo, sizeof(CurrentTo));

  int l = strlen(Daten.cdir)+1; // mit term 0
  out.write((char *)&l, sizeof(l));
  out.write((char *)&Daten.cdir, l * sizeof(char));

  int x, y, w, h;
  gdk_window_get_root_origin(Application.DLG_Main->window, &x, &y);
  gdk_window_get_size(Application.DLG_Main->window, &w, &h);

  out.write((char *)&x, sizeof(x));
  out.write((char *)&y, sizeof(y));
  out.write((char *)&w, sizeof(w));
  out.write((char *)&h, sizeof(h));
2019*/  
}

std::ostream& operator<<(std::ostream& out, Verstaerk ver)
{
// Achtung hier kein << Verstaerk, dann rekursion!

  return out << (int)ver << "   # Ampl: "
	         << "Dagan=" << (int)Dagan
		     << "; EPC7=" << (int)EPC7;
}

std::ostream& operator<<(std::ostream& out, Mittelwertbildung mit)
{
// Achtung hier kein << Mittelwertbildung, dann rekursion!

  return out << (int)mit << "   # Mittelwertbildung: "
	         << "keine=" << (int)keine
		     << "; normal=" << (int)normal
			 << "; gleitend=" << (int)gleitend;
}


std::ostream& operator<<(std::ostream& out, TSettings& Settings)
{
	
int x, y, w, h;

//2019  gdk_window_get_root_origin(Application.DLG_Main->window, &x, &y);
//2019  gdk_window_get_size(Application.DLG_Main->window, &w, &h);

  // abspeichern ascii //z.b. ini file
  return out << "#settings" << std::endl
	         << Settings.Ampl << std::endl
	         << Settings.VFactor << "   # VFactor" << std::endl
             << Settings.Gain    << "   # Gain" << std::endl
             << Settings.Samplininterval         << "   # Samplininterval" << std::endl
             << Settings.PositiveMembranPotential << "   # PositiveMembranPotential" << std::endl
             << Settings.Zoomed.StackedWindows    << "   # StackedWindows (Zoomed)" << std::endl
             << Settings.Zoomed.AvarageFilter     << " (Zoomed)" << std::endl
             << Settings.Zoomed.LengthOfAvarageFilter   << "   # LengthOfAvarageFilter (Zoomed)" << std::endl
             << Settings.UnZoomed.StackedWindows        << "   # StackedWindows (UnZoomed)" << std::endl
             << Settings.UnZoomed.AvarageFilter         << " (UnZoomed)" << std::endl
             << Settings.UnZoomed.LengthOfAvarageFilter << "   # LengthOfAvarageFilter (UnZoomed)" << std::endl
             << Settings.enlarge_current  << "   # enlarge_current" << std::endl
             << Settings.Zoom             << "   # Zoom" << std::endl
             << Settings.SamplesToSkip    << "   # SamplesToSkip" << std::endl
             << Settings.SamplesToProcess << "   # SamplesToProcess" << std::endl
             << Settings.CurrentFrom      << "   # CurrentFrom" << std::endl
             << Settings.CurrentTo        << "   # CurrentTo" << std::endl
             << Settings.multi_samples << "   # Multi_Samples" << std::endl
			 << x << "   # x" << std::endl
		     << y << "   # y" << std::endl
			 << w << "   # width" << std::endl
		     << h << "   # height" << std::endl
			 // bins per log erstma fuer alle 1d gleich 
			 // sonst geh�ren diese werte auch nicht in die settings!
			 << Daten.Dwell_1d_A.bins_per_log() << "   # Dwell_1d bins_per_log" << std::endl
			 // auch fuer alle 2d gleich!
			 << Daten.Dwell_2d_A.bins_per_log() << "   # Dwell_2d bins_per_log" << std::endl
			 << Daten.cdir << "   # current dir" << std::endl;
		 
return out<<"false"<<std::endl; //2019
}

std::istream& operator>>(std::istream& in, TSettings& Settings)
{
char buffer[MaxTextLen];
char c;
int  i;
short buf_i16;

  in.get(buffer, MaxTextLen, '\n');
  if (!strncmp(buffer, "#settings", 9))	//2019 strlwr removed
  {
	while (in.get(c)) if (c == '\n') break;
	in >> i;
	Settings.Ampl = (Verstaerk)i;
	while (in.get(c)) if (c == '\n') break;
	in >> Settings.VFactor;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.Gain;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.Samplininterval;
	while (in.get(c)) if (c == '\n') break;
	in >> Settings.PositiveMembranPotential;
	while (in.get(c)) if (c == '\n') break;
	in >> Settings.Zoomed.StackedWindows;
	while (in.get(c)) if (c == '\n') break;
    in >> i;
	Settings.Zoomed.AvarageFilter = (Mittelwertbildung)i;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.Zoomed.LengthOfAvarageFilter;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.UnZoomed.StackedWindows;
	while (in.get(c)) if (c == '\n') break;
    in >> i;
	Settings.UnZoomed.AvarageFilter = (Mittelwertbildung)i;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.UnZoomed.LengthOfAvarageFilter;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.enlarge_current;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.Zoom;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.SamplesToSkip;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.SamplesToProcess;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.CurrentFrom;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.CurrentTo;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.multi_samples;
	while (in.get(c)) if (c == '\n') break;
	in >> Settings.x;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.y;
	while (in.get(c)) if (c == '\n') break;
    in >> Settings.w;
	while (in.get(c)) if (c == '\n') break;
	in >> Settings.h;
	while (in.get(c)) if (c == '\n') break;
	in >> buf_i16;
	Daten.Dwell_1d_A.bins_per_log(buf_i16);
	while (in.get(c)) if (c == '\n') break;
	in >> buf_i16;
	Daten.Dwell_2d_A.bins_per_log(buf_i16);
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.cdir;
    while (in.get(c)) if (c == '\n') break;
    	
    	
  in >> Daten.Geneticfitparameter.init_random;		//2021 Tobias
  	while (in.get(c)) if (c == '\n') break; 
  in >> Daten.Fitparameter.log_min_close;					//2021 Tobias
  	while (in.get(c)) if (c == '\n') break; 
  in >> Daten.Fitparameter.log_max_close;					//2021 Tobias
  	while (in.get(c)) if (c == '\n') break; 
  in >> Daten.Fitparameter.log_min_open;					//2021 Tobias
  	while (in.get(c)) if (c == '\n') break; 
  in >> Daten.Fitparameter.log_max_open;					//2021 Tobias
  	while (in.get(c)) if (c == '\n') break; 
  in >> Daten.Fitparameter.bins_per_log;					//2021 Tobias
  	while (in.get(c)) if (c == '\n') break;   		  		  		  		  		
  		 	
	in >> Daten.Geneticfitparameter.popsize;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.ngen;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.elitism;	
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.pmut;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.pcross;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.bits;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.nConvergence;
    while (in.get(c)) if (c == '\n') break; 	
	in >> Daten.Geneticfitparameter.methode;
    while (in.get(c)) if (c == '\n') break;
  in >> Daten.calculate_error;
    while (in.get(c)) if (c == '\n') break;   	
	in >> Daten.Geneticfitparameter.fit_methode;
    while (in.get(c)) if (c == '\n') break;
     in >>Daten.Fitparameter.multiplikator;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.proportionalfit;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.logscale;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.SimSettings.show_2D_Diff;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.SimSettings.dosbox;	
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.nr_calc_mean_result;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[1];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[1];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[2];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[2];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[3];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[3];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[4];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[4];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[5];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[5];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[6];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[6];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[7];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[7];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[8];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[8];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[9];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[9];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmin[10];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.startwertmax[10];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.levelfit;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.levelfit_min[1];	
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.levelfit_max[1];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.levelfit_min[2];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.levelfit_max[2];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[1];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[2];				
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[3];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[4];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[5];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[6];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[7];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[8];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[9];
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.state_conductance[10];
		while (in.get(c)) if (c == '\n') break;
	in >> Daten.SaveTimeseries.min_noise;
	  while (in.get(c)) if (c == '\n') break;
	in >> Daten.SaveTimeseries.max_noise;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.SimSettings.sampling_frequency;	
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.SimSettings.filter_frequency;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.automated_loading;   
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.mpi_support;
    while (in.get(c)) if (c == '\n') break;	
  in >> Daten.SimSettings.result_file_name;
    while (in.get(c)) if (c == '\n') break;		
	in >> Daten.Geneticfitparameter.files_directory;
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.file_time_series;
    while (in.get(c)) if (c == '\n') break;
  in >> Daten.Geneticfitparameter.noise_generation;
    while (in.get(c)) if (c == '\n') break;	
	in >> Daten.Geneticfitparameter.file_noise_series;		//or spectral file
    while (in.get(c)) if (c == '\n') break;
	in >> Daten.Geneticfitparameter.file_setfile;  
    while (in.get(c)) if (c == '\n') break;
    	
  in >> Daten.SaveTimeseries.automated;  
    while (in.get(c)) if (c == '\n') break;
  in >> Daten.SaveTimeseries.folder;  
    while (in.get(c)) if (c == '\n') break;
  in >> Daten.SaveTimeseries.number;  
    while (in.get(c)) if (c == '\n') break;
  in >> Daten.SaveTimeseries.random_rate;  
    while (in.get(c)) if (c == '\n') break;
  in >> Daten.SaveTimeseries.random_lvl;  
    while (in.get(c)) if (c == '\n') break;
  in >> Daten.SaveTimeseries.random_noise;  
    while (in.get(c)) if (c == '\n') break;	  	
    	
	in >> Daten.SimSettings.check;
    if (!Daten.Geneticfitparameter.mpi_support)	
       std::cout<<"settings last read: "<<Daten.SimSettings.check<<std::endl;	
  } else {
  	std::cout<<"setfile.ini corrupted"<<std::endl;
  	std::this_thread::sleep_for(std::chrono::milliseconds(5000));	
    in.clear(std::ios::badbit); // status setzen
  }

  // kein kommentar als pfad!
  if (Daten.cdir[0] == '#') Daten.cdir[0] = 0; 

  // bins per log fuer alle 1d gleich und fuer alle 2d gleich!
  Daten.Dwell_1d_A.init((char *)"A (measured)");
  Daten.Dwell_2d_A.init((char *)"A (measured)");
  Daten.Dwell_2d_A.addjumps(&Daten.Dwell_1d_A);

  Daten.Dwell_1d_B.init((char *)"B (Simulated)");
  Daten.Dwell_2d_B.init((char *)"B (Simulated)");
  Daten.Dwell_2d_B.addjumps(&Daten.Dwell_1d_B);

  Daten.Dwell_2d_diff.init((char *)" ");

  return in;
}


TColorAndOrderRec::TColorAndOrderRec()
{
  HOHD_R = HOHD_default_Color_R;
  HOHD_G = HOHD_default_Color_G;
  HOHD_B = HOHD_default_Color_B;
  SHD_R = SHD_default_Color_R;
  SHD_G = SHD_default_Color_G;
  SHD_B = SHD_default_Color_B;
  HOSHD_R = HOSHD_default_Color_R;
  HOSHD_G = HOSHD_default_Color_G;
  HOSHD_B = HOSHD_default_Color_B;
  DrawingOrder = 0;
  SHD_Showflag = true;
  HOHD_Showflag = true;
  HOSHD_Showflag = true;
}


TMarkBits::TMarkBits(unsigned short aMark, int aStartPos, bool achosen)
{
  mark=aMark;
  startpos=aStartPos;
  //TOBIAS begin
  chosen=achosen;
  //Tobias end
  next=NULL;
}

// Klassenvariablen definieren
short TSprung::i_oldwindow = -1;
short TSprung::x_pixel_old = 0;
short TSprung::y_pixel_old = 0;
short TSprung::y_pixel_min = 0;
short TSprung::y_pixel_max = 0;

int TSprung::p = 0;
short TSprung::j = 0;
int TSprung::pos1 = 0;
int TSprung::pos2 = 0;
bool TSprung::YetStarted = true;
double TSprung::SD = 0.0;

TSprung::TSprung()
{ 
  level_1(0, 0);
  level_2(0, 0);
  position(0);
  setnvs(0, 0, 0);
}

void
TSprung::load(std::istream& in)
{
  in.read((char *)&sprung_i_null, sizeof(sprung_i_null));
  in.read((char *)&sprung_i_channel, sizeof(sprung_i_channel));
  in.read((char *)&sprung_i_extra_channel, sizeof(sprung_i_extra_channel));
  in.read((char *)&sprung_C_1, sizeof(sprung_C_1));
  in.read((char *)&sprung_U_1, sizeof(sprung_U_1));
  in.read((char *)&sprung_C_2, sizeof(sprung_C_2));
  in.read((char *)&sprung_U_2, sizeof(sprung_U_2));
  in.read((char *)&sprung_position, sizeof(sprung_position));  
}

void
TSprung::store(std::ostream& out)
{
  out.write((char *)&sprung_i_null, sizeof(sprung_i_null));
  out.write((char *)&sprung_i_channel, sizeof(sprung_i_channel));
  out.write((char *)&sprung_i_extra_channel, sizeof(sprung_i_extra_channel));
  out.write((char *)&sprung_C_1, sizeof(sprung_C_1));
  out.write((char *)&sprung_U_1, sizeof(sprung_U_1));
  out.write((char *)&sprung_C_2, sizeof(sprung_C_2));
  out.write((char *)&sprung_U_2, sizeof(sprung_U_2));
  out.write((char *)&sprung_position, sizeof(sprung_position));
}

void
TSprung::plot(short strom_min, short strom_max)
{
/*2019	
// aus bigcolle.pas
short   i,i_window = 0;
double  rel_height, rel_pos;
short   y_pixel,y_pixel_2,x_pixel,window_height;
double  a,b;//,c;

  window_height = Daten.g_maxy / Daten.Settings.StackedWindows();

  a = position();
  b = Daten.Settings.SamplesToSkip+1;
  a = a-b;
  b = Daten.Settings.StackedWindows();
  a = a*b;
  b = Daten.Settings.SamplesToProcess;
  a = a/b;
  rel_pos = a;
  // Dies fkt. aus Genauigkeitsgr�nden nicht:   !!!!
  // rel_pos = (Sprung.Position-Daten.Settings.SamplesToSkip+1)*Daten.Settings.StackedWindows()/Daten.Settings.SamplesToProcess;

  if ((rel_pos > 0) && (rel_pos < Daten.Settings.StackedWindows()))
  {
    for (i=0;i <= Daten.Settings.StackedWindows()-1;i++)
    {
      if ((i<=rel_pos) && (rel_pos<=(i+1))) i_window = i;
    }
    rel_pos = rel_pos - i_window;
    rel_height = (level_1() - strom_min)/(strom_max - strom_min);
    y_pixel = Daten.g_maxy - (short)(rel_height * window_height) + (i_window - Daten.Settings.StackedWindows() + 1) * window_height;
    x_pixel = round(rel_pos * Daten.g_maxx);
    rel_height = (level_2() - strom_min)/(strom_max - strom_min);
    y_pixel_2 = Daten.g_maxy - (short)(rel_height * window_height) + (i_window - Daten.Settings.StackedWindows() + 1) * window_height;

	// muss der sprung noch gezeichnet werden, 
	// oder ist an dieser stelle schon einer?
    bool draw = (x_pixel_old != x_pixel);
	short miny = MIN(y_pixel,y_pixel_2);
	short maxy = MAX(y_pixel,y_pixel_2);
	if (x_pixel_old != x_pixel)
	{
	  draw = true;
      y_pixel_min = miny;
	  y_pixel_max = maxy;
	} else {
      if (maxy > y_pixel_max)
	  {
		y_pixel_max = maxy;
		draw = true;
	  }
      if (miny < y_pixel_min) 
	  {
	    y_pixel_min = miny;
		draw = true;
	  }
	}	

	if (i_window != i_oldwindow)
    {
      // sprung muss gezeichnet werden, auch wenn x_Pixel_old == x_pixel
	  draw = true;
      if (true) // nur wenn wirklich alle spruenge erkannt
      {
        if (i_oldwindow != -1)
        {
      	  // niveau nach letztem sprung im vorherigen stackedwindow zeichnen
          gdk_draw_line(Daten.DataPix, Application.draw_gc,
                        x_pixel_old, y_pixel_old,
                        Daten.g_maxx, y_pixel_old);
        }
        // niveau vor erstem sprung im stackedwindow zeichnen 	
        gdk_draw_line(Daten.DataPix, Application.draw_gc,
                      0, y_pixel,
                      x_pixel, y_pixel);
        x_pixel_old = x_pixel;
        y_pixel_old = y_pixel;
      }
    } else {
      // niveau zwischen zwei spruengen zeichnen 
      // LineTo(PaintDC,x_pixel,y_pixel);
	  if (draw)
	  {
        gdk_draw_line(Daten.DataPix, Application.draw_gc,
                      x_pixel_old, y_pixel_old,
                      x_pixel, y_pixel_old);
        x_pixel_old = x_pixel; // wenn !draw, dann schon gleich 
	  }
    }
    // sprung zeichnen
	if (draw)
	{
      gdk_draw_line(Daten.DataPix, Application.draw_gc,
                      x_pixel_old, y_pixel_old,
                      x_pixel_old, y_pixel_2);
      x_pixel_old = x_pixel;  // wenn !draw, dann schon gleich
	}
    i_oldwindow = i_window;
    y_pixel_old = y_pixel_2;
  }
2019*/  
}

void
TSprung::reset(bool draw) // Tobias: = true weg
{
/*2019	
  if (draw && (x_pixel_old != 0))
  {
    // hier wird niveau nach letztem sprung gezeichnet
    gdk_draw_line(Daten.DataPix, Application.draw_gc,
                  x_pixel_old, y_pixel_old,
                  Daten.g_maxx, y_pixel_old);
  }

  i_oldwindow = -1;
  x_pixel_old = 0;
  y_pixel_old = 0;
  y_pixel_min = 0;
  y_pixel_max = 0;
2019*/  
}

void
TSprung::calcStD(double level)
{
  int i;
  double l2;

  if (YetStarted)
  {
    YetStarted = false;
    pos2 = position();
    for (i=pos1;i <= pos2;i++)
    {
      //Hier die kritische PASCAL Typeumwandlung 
		SD = SD + (level - Daten.SimulatWert[i]) * (level - Daten.SimulatWert[i]);
    }
    p = p + (pos2 - pos1);
  } else {
	l2 = level_2();
    if (l2 == level) YetStarted = true;
    pos1 = position();
  }
}

double  
TSprung::calcSD(double level, std::vector<TSprung> &spruenge)
{
std::vector<TSprung>::iterator iter;

  pos1 = 0;
  pos2 = 0;
  YetStarted = false;
  SD = 0;
  p = 0;

  for(iter=spruenge.begin();iter != spruenge.end();iter++) 
  {
    (*iter).calcStD(level);
  }

  if (p != 0)
  {
    return sqrt(fabs(SD / p));
  } else {
    return 0.0;
  }
}


TDwell::TDwell(double ax, double ay)
{
  x=ax;
  y=ay;
}

void
TDwell::load(std::istream& in)
{
  in.read((char *)&x, sizeof(x));
  in.read((char *)&y, sizeof(y));
}

void
TDwell::store(std::ostream& out)
{
  out.write((char *)&x, sizeof(x));
  out.write((char *)&y, sizeof(y));
}

Tniveau::Tniveau(double f, short i, short j)
{
  n=f;
  channel=i;
  extra_channel=j;
}

void
Tniveau::load(std::istream& in)
{
  in.read((char *)&n, sizeof(n));
  in.read((char *)&channel, sizeof(channel));
  in.read((char *)&extra_channel, sizeof(extra_channel));
}

void
Tniveau::store(std::ostream& out)
{
  out.write((char *)&n, sizeof(n));
  out.write((char *)&channel, sizeof(channel));
  out.write((char *)&extra_channel, sizeof(extra_channel));
//  cout << "speichere Niveau: " << n << endl;
}

TRate_matrix::TRate_matrix()
{
short i,j;

  for (i=0;i <= n_states_max-1;i++)
  {
    for (j=0;j <= n_states_max-1;j++)
    {
      r[i][j] = 0;
	}
    open[i] = false;
  }
}
